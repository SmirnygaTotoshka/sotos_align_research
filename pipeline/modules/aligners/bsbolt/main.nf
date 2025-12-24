process ALIGN {
    tag "${meta.id}"
    label "bsbolt"

    container "eod-tools.med-gen.ru/bsbolt:latest"

    input:
        tuple val(meta), path(reads)
        path reference_dir
    output:
        tuple val(meta), path("${meta.id}.bam")
        
    script:
    if (meta.reverse == ""){
        def reads_params = "-F1 ${reads[0]}"
    } else {
        def reads_params = "-F1 ${reads[0]} -F2 ${reads[1]}"
    }
    """
    bsbolt Align ${reads_params} -DB ${reference_dir} -t ${task.cpus} -R "@RG ID:${meta.id} SM:${meta.id}" -UN   -O ${meta.id}.bam
    """
    stub:
    """
    touch ${meta.id}.bam
    """
}