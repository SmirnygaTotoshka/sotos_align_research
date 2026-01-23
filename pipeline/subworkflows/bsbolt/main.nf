process ALIGN{
    tag "${meta.id}"
    label 'align'
    container "eod-tools.med-gen.ru/sotos-align:1.0"
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
        -t ${task.cpus} \
        ${is_nondirection_option}
    """
    stub:
    """
    touch ${meta.id}_raw.bam
    """ 
}
process CLEAN{
    tag "${meta.id}"
    label 'samtools'
    container "eod-tools.med-gen.ru/sotos-align:1.0"
    input:
        tuple val(meta), path(raw_alignment)
        path panel_bed
    output:
        tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
    
    script:
    """
    samtools view -b -h ${params.samtools_flags} -@ ${task.cpus} \
    -e "mapq >= ${params.min_mapq} && sclen < ${params.softclipped_threshold} && hclen < ${params.hardclipped_threshold}" \
    ${raw_alignment} | samtools sort - -O BAM -o filtered.bam -@ ${task.cpus}

    samtools index -M -@ ${task.cpus} filtered.bam

    samtools view -H filtered.bam | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > tmp.genome

    primers_max_len=\$(awk '/^>/{if(l>m)m=l; l=0; next} {gsub(/[[:space:]]/,""); l+=length(\$0)} END{if(l>m)m=l; print m}' ${params.primers_fasta})
    bedtools flank -i ${panel_bed} -g tmp.genome -l \${primers_max_len} -r \${primers_max_len} | \
    awk 'BEGIN{FS="\\t"; OFS="\\t"} {print(\$1,\$2,\$3,".",0,".")}' | \
    sort -k 1,1 -k2,2n > primers.bed

    samtools ampliconclip -@ ${task.cpus} \
    --hard-clip --both-ends \
    --filter-len ${params.read_min_len} \
    -b primers.bed \
    filtered.bam > clipped.bam

    samtools sort -@ ${task.cpus} -n clipped.bam > byname.bam
    python ${workflow.projectDir}/scripts/remove_singletons.py byname.bam byname_cleaned.bam
    samtools sort -@ ${task.cpus} byname_cleaned.bam > ${meta.id}.bam
    samtools index -M -@ ${task.cpus} ${meta.id}.bam
    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.bam.bai
    """ 
}
process METH_CALL{
    tag "${meta.id}"
    label 'align'
    container "eod-tools.med-gen.ru/sotos-align:1.0"
    input:
        tuple val(meta), path(alignment), path(index)
        path reference_dir

    output:
        tuple val(meta), path(alignment), path(index), path("${meta.id}.CGmap.gz"), emit: for_startistic
        tuple val(meta), path("${meta.id}.CGmap.gz"), emit: cgmap
        tuple val(meta), path("${meta.id}.bg.gz"), emit: track
    
    script:
    def min_base_quality = Math.max(params.cut_front, params.cut_tail)
    """

    num_alignments=\$(samtools view -c ${alignment})

    if [ \$num_alignments -ge ${params.min_call_depth} ]; then
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
    else
        touch ${meta.id}.CGmap.gz
        touch ${meta.id}.bg.gz
    fi
    """
    stub:
    """
    touch ${meta.id}.CGmap.gz
    touch ${meta.id}.bg.gz
    """
}
process CALC_STAT{
    tag "${meta.id}"
    label 'single'
    container "eod-tools.med-gen.ru/sotos-align:1.0"
    input:
        tuple val(meta), path(alignment), path(index), path(cgmap)
        path panel
        path cpg_context
        path conversion_context
    output:
        tuple val(meta), path('*.csv')
    
    script:
    """
        python ${workflow.projectDir}/scripts/calculate_statistics.py \
        --input ${alignment} \
        --output ${meta.id}_statistics.csv \
        --panel ${panel} \
        --cpg ${cpg_context} \
        --conversion ${conversion_context} \
        --cgmap ${cgmap} \
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
