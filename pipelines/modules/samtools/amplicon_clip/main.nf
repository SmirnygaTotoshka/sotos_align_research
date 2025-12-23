process AMPLICON_CLIP{
    tag "$meta.id"
    label 'validation_filter'
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(alignment)
        path panel
    output:
        tuple val(meta), path("${meta.id}_clip.bam"), path("${meta.id}_clip.bam.bai"), emit: alignments
    script:
    def fix_mate = meta.reverse == "" ? "" : "| samtools fixmate -m - -"
    """
        samtools sort -@ ${task.cpus} -O BAM ${alignment} > ${meta.id}_before_filter.bam
        samtools index -@ ${task.cpus} ${meta.id}_before_filter.bam

        samtools view -H ${alignment} | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > tmp.genome

        bedtools flank -i ${panel} -g tmp.genome -l ${params.primers_maxlen} -r 0 | \
        awk 'BEGIN{FS="\\t"; OFS="\\t"} {print(\$1,\$2,\$3,".",0,"+")}' > left_primers.bed

        bedtools flank -i ${panel} -g tmp.genome -r ${params.primers_maxlen} -l 0 | \
        awk 'BEGIN{FS="\\t"; OFS="\\t"} {print(\$1,\$2,\$3,".",0,"-")}' > right_primers.bed
        cat left_primers.bed right_primers.bed | sort -k 1,1 -k2,2n > primers.bed
        
        samtools ampliconclip -@ ${task.cpus} --strand --hard-clip --both-ends --filter-len ${params.read_min_len} \
        -b primers.bed ${meta.id}_before_filter.bam $fix_mate | samtools sort -@ ${task.cpus} -O bam - > ${meta.id}_clip.bam
        
        samtools index -@ ${task.cpus} ${meta.id}_clip.bam
    """
    stub:
    """
    touch ${meta.id}_clip.bam
    touch ${meta.id}_clip.bam.bai
    """    
}