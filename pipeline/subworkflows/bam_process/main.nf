include { samplesheetToList } from 'plugin/nf-schema'


process FILTER_BAM {
    tag "${meta.id}"
    label 'samtools'
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path("*.bam", arity : 1, name: "raw.bam")
    output:
        tuple val(meta), path("${meta.id}_filtered.bam", arity : 1), path("${meta.id}_filtered.bam.bai", arity : 1)
    script:
    def cpus = task.cpus > 1 ? task.cpus.intdiv(2) : 1
    """
    samtools view -b -h ${params.samtools_flags} -@ ${cpus} \
    -e 'mapq >= ${params.min_mapq} && sclen < ${params.softclipped_threshold} && hclen < ${params.hardclipped_threshold}' \
    raw.bam | samtools sort - -O BAM -o ${meta.id}_filtered.bam -@ ${cpus}

    samtools index -@ ${task.cpus} ${meta.id}_filtered.bam
    """
    stub:
    """
    touch ${meta.id}_filtered.bam
    touch ${meta.id}_filtered.bam.bai
    """
}

process SPLIT_BAM {
    tag "${meta.id}"
    label 'samtools'
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(alignment), path(index), val(locus)
    output:
        tuple val(meta), path("*.bam"), val(locus)
    script:
    """
    echo "${locus.chr}\t${locus.start}\t${locus.end}" > regions.bed
    bedtools intersect -a ${alignment} -b regions.bed -wa > ${meta.id}_${locus.name}.bam
    """
    stub:
    """
    touch ${meta.id}_${locus.name}.bam
    """
}

process AMPLICON_CLIP {
    tag "${meta.id}"
    label 'samtools'
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path("*.bam", arity : 1, name: "locus.bam"), val(locus)
    output:
        tuple val(meta), path("${meta.id}_${locus.name}_clipped.bam")
    script:
    """
    echo "${locus.chr}\t${locus.start}\t${locus.end}" > regions.bed
    
    samtools sort -@ ${task.cpus} -O BAM locus.bam > before_filter.bam
    samtools index -@ ${task.cpus} before_filter.bam

    samtools view -H locus.bam | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > tmp.genome

    forward_primer="${locus.forward}"
    reverse_primer="${locus.reverse}"

    bedtools flank -i regions.bed -g tmp.genome -l \${#forward_length} -r 0 | \
    awk 'BEGIN{FS="\t"; OFS="\t"} {print(\$1,\$2,\$3,".",0,"+")}' > left_primers.bed

    bedtools flank -i regions.bed -g tmp.genome -r \${#reverse_length} -l 0 | \
    awk 'BEGIN{FS="\t"; OFS="\t"} {print(\$1,\$2,\$3,".",0,"-")}' > right_primers.bed

    cat left_primers.bed right_primers.bed | sort -k 1,1 -k2,2n > primers.bed
        
    samtools ampliconclip -@ ${task.cpus} \
    --strand --hard-clip --both-ends \
    --filter-len ${params.read_min_len} \
    -b primers.bed \
    before_filter.bam > ${meta.id}_${locus.name}_clipped.bam 
    """
    stub:
    """
    touch ${meta.id}_${locus.name}_clipped.bam
    """
}

process MERGE_CLIPPED_LOCUSES {
    tag "${meta.id}"
    label 'samtools'
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(bam_files)
    output:
        tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
    script:
    """
    samtools merge -o ${meta.id}.bam -@ ${task.cpus} ${bam_files}
    samtools index -@ ${task.cpus} ${meta.id}.bam
    """
    stub:
    """
    touch ${meta.id}.bam
    touch ${meta.id}.bam.bai
    """
}

workflow BAM_PROCESS{
    take:
        alignment
        panel
    main:
        validate_panel = Channel.fromList(samplesheetToList(panel, "assets/schema_panel.json"))
        merged_bam = FILTER_BAM(alignment).combine(validate_panel) | 
                     SPLIT_BAM                                     | 
                     AMPLICON_CLIP                                 |
                     groupTuple                                    |
                     MERGE_CLIPPED_LOCUSES             
    emit:
        merged_bam
}