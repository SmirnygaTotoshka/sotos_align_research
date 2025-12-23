process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: reads
    tuple val(meta), path('*.zip'), emit: json
    tuple val(meta), path('*.html'), emit: report

    script:
    def files = reads.collect().join(' ')
    """
        fastqc -t ${task.cpus} ${files}
    """
    stub:
    """
        touch ${meta.id}.html
        touch ${meta.id}.zip
    """
}
