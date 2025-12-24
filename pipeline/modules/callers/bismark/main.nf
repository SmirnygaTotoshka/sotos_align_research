process METHYLATION_EXTRACTOR{
    tag "${meta.id}"
    label "bismark"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(alignment, arity : 1) 

    output:
        tuple val(meta), path("*.gz"), emit: tracks 
        tuple val(meta), path("*.txt"), emit: statistics

    afterScript "rm ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam"

    script:
    def mode = meta.reverse == "" ? "--single-end" : "--paired-end"
    """
        samtools sort -@ ${task.cpus} -n ${alignment} > ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam
        bismark_methylation_extractor ${task.ext.bismark_param_str} ${mode} ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam 
    """
    stub:
    """
    touch CHG.txt
    touch CHH.txt
    touch CpG.txt
    touch chrom_sizes.txt
    touch mbias.txt
    touch bismark.cov.gz
    touch UCSC.bedGraph.gz
    touch sorted_byname.bedGraph.gz
    touch splitting_report.txt
    """
}