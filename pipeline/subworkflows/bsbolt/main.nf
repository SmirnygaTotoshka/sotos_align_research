process ALIGN{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(reads, arity: 1..2)
        path reference_dir
    output:
        tuple val(meta), path('*.bam')
    
    script:
    def input_options = meta.reverse.size() == 0 ? "-F1 ${reads[0]}" : "-F1 ${reads[0]} -F2 ${reads[1]}"
    def is_nondirection_option = params.non_directional ? "-UN" : ""       
    """
        python -m bsbolt Align \
        -DB ${reference_dir} \
        ${input_options} \
        -O ${meta.id}_raw \
        -R "@RG ID:${meta.id} SM:${meta.id}" \
        -t ${task.cpus}
        ${is_nondirection_option}
    """
    stub:
    """
    touch ${meta.id}_raw.bam
    """ 
}
process CLEAN{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(raw_alignment)
        path panel_bed
    output:
        tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
    
    script:
    def primers_max_len = getMaxSeqLen(panel_bed)
    """
    
    samtools view -b -h ${params.samtools_flags} -@ ${task.cpus} \
    -e "mapq >= ${params.min_mapq} && sclen < ${params.softclipped_threshold} && hclen < ${params.hardclipped_threshold}" \
    ${raw_alignment} | samtools sort - -O BAM -o sot011_filtered.bam -@ ${task.cpus}

    samtools index ${task.cpus} sot011_filtered.bam

    samtools view -H sot011_filtered.bam | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > tmp.genome

    bedtools flank -i ${panel_bed} -g tmp.genome -l ${primers_max_len} -r ${primers_max_len} | \
    awk 'BEGIN{FS="\\t"; OFS="\\t"} {print(\$1,\$2,\$3,".",0,".")}' | \
    sort -k 1,1 -k2,2n > primers.bed

    samtools ampliconclip -@ ${task.cpus} \
    --hard-clip --both-ends \
    --filter-len ${params.read_min_len} \
    -b primers.bed \
    sot011_filtered.bam > sot011_clipped.bam

    samtools sort -@ ${task.cpus} -n sot011_clipped.bam > sot011_byname.bam
    python ${workflow.projectDir}/scripts/remove_singletons.py sot011_byname.bam sot011_byname_cleaned.bam
    samtools sort -@ ${task.cpus} sot011_byname_cleaned.bam > sot011.bam
    samtools index -@ ${task.cpus} sot011.bam
    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.bam.bai
    """ 
}
process METH_CALL{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(alignment), path(index)
        path reference_dir

    output:
        tuple val(meta), path(alignment), path("${meta.id}.CGmap.gz"), emit: for_startistic
        tuple val(meta), path("${meta.id}.CGmap.gz"), emit: cgmap
        tuple val(meta), path("${meta.id}.bg"), emit: track
    
    script:
    def min_base_quality = Math.max(params.cut_front, params.cut_tail)
    """
        # Call CGmap output
        python -m bsbolt CallMethylation \
        -I ${alignment} \
        -O ${meta.id} \
        -DB ${reference_dir} \
        -t ${task.cpus} \
        -max ${params.max_call_depth} \
        -min ${params.min_call_depth} \
        -BQ ${min_base_quality} \
        -MQ ${params.min_mapq} \
        -IO

        #repeat call for bedGraph output. It cannot generate CGmap and bedGraph for one call
        python -m bsbolt CallMethylation \
        -I ${alignment} \
        -O ${meta.id} \
        -DB ${reference_dir} \
        -t ${task.cpus} \
        -max ${params.max_call_depth} \
        -min ${params.min_call_depth} \
        -BQ ${min_base_quality} \
        -MQ ${params.min_mapq} \
        -IO \
        -BG
    """
    stub:
        """
        touch ${meta.id}.CGmap.gz
        touch ${meta.id}.bg
        """
}
process CALC_STAT{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(alignment), path(cgmap)
        path panel
        path cpg_context
        path conversion_context
    output:
        tuple val(meta), path('*.csv')
    
    script:
    """
        python ${workflow.projectDir}/scripts/calculate_statistics.py \
        --input ${alignment}
        --output ${meta.id}_statistics.csv
        --panel ${panel}
        --cpg ${cpg_context}
        --conversion ${conversion_context}
        --cgmap ${cgmap}
        --depth_threshold ${params.depth_threshold} 
    """
    stub:
    """
    touch ${meta.id}_statistics.csv
    """ 
}

workflow BSBOLT{
    take:
        reads
        reference
        cpg_context
        conversion_context
    main:
        ref_value = reference.first()
        //cpg_value = cpg_context.first()
        //conversion_value = conversion_context.first()

        raw_alignment = ALIGN(reads, ref_value)
        cleaned_alignment = CLEAN(raw_alignment, params.panel_bed)
        calling_result = METH_CALL(cleaned_alignment, ref_value)
        statistics = CALC_STAT(calling_result.for_startistic, params.panel_bed, cpg_context, conversion_context)
    emit:
        alignment = cleaned_alignment
        cgmap = calling_result.cgmap
        track = calling_result.track
        by_locus = statistics
}

def getMaxSeqLen(Path fasta_file){
    def file = new File(fasta_file)
    def sequences = [:]
    def currentHeader = ""
    def currentSeq = ""

    file.eachLine { line ->
        if (line.startsWith('>')) {
            if (currentHeader) sequences[currentHeader] = currentSeq
            currentHeader = line.substring(1).trim()
            currentSeq = ""
        } else {
            currentSeq += line.trim()
        }
    }
    if (currentHeader) sequences[currentHeader] = currentSeq

    return sequences.max { record -> record.value.length() }
}