process ALIGN {
    tag "${meta.id}"
    label "samtools"

    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(reads)
        path reference_file
    output:
        tuple val(meta), path("${meta.id}_preprocessing.bam")

    script:
    """
    /new_meth/bwa-meth-master/bwameth.py --reference ${reference_file} ${reads} | samtools sort - -O BAM -@ ${task.cpus} -o ${meta.id}_preprocessing.bam
    """
    stub:
    """
    touch ${meta.id}_preprocessing.bam
    """
}