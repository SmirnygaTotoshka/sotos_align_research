process ALIGN{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(reads, arity: 1..2)
        path reference_dir
    output:
        tuple val(meta), path('*.bam')
    
    script:
    def input_options = meta.reverse == "" ? "-F1 ${reads[0]}" : "-F1 ${reads[0]} -F2 ${reads[1]}"
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

    samtools view -H sot011_filtered.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > tmp.genome

    bedtools flank -i ${panel_bed} -g tmp.genome -l ${primers_max_len} -r ${primers_max_len} | \
    awk 'BEGIN{FS="\t"; OFS="\t"} {print($1,$2,$3,".",0,".")}' | \
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
    """ 
}
process METH_CALL{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(reads, arity: 1..2)

    output:
        tuple val(meta), path('*.fastq.gz')
    
    script:
    def io_options = meta.reverse == "" ? "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz" : "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz --in2 ${reads[1]} --out2 ${meta.id}_R2.fastq.gz"      
    """
        bash /bbmap/repair.sh ${io_options} 
    """
    stub:
    if (meta.reverse == ""){
        """
        touch ${meta.id}_R1.fastq.gz
        """ 
    }
    else {
        """
        touch ${meta.id}_R1.fastq.gz
        touch ${meta.id}_R2.fastq.gz
        """
    }
}
process CALC_STAT{
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(reads, arity: 1..2)

    output:
        tuple val(meta), path('*.fastq.gz')
    
    script:
    def io_options = meta.reverse == "" ? "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz" : "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz --in2 ${reads[1]} --out2 ${meta.id}_R2.fastq.gz"      
    """
        bash /bbmap/repair.sh ${io_options} 
    """
    stub:
    if (meta.reverse == ""){
        """
        touch ${meta.id}_R1.fastq.gz
        """ 
    }
    else {
        """
        touch ${meta.id}_R1.fastq.gz
        touch ${meta.id}_R2.fastq.gz
        """
    }
}

workflow BSBOLT{
    take:
        reads
        bisulfite_reference
    main:
        raw_alignment = ALIGN(reads, bisulfite_reference)
        cleaned_alignment = CLEAN(raw_alignment, params.panel_bed)
        calling_result = METH_CALL(cleaned_alignment)
        statistics = CALC_STAT(cgmap)
    emit:
        alignment = cleaned_alignment
        cgmap = calling_result.cgmap
        track = calling_result.track
        by_locus = statistics.by_locus
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

    return sequences.max { it.value.length() }
}