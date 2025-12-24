process METH_CALL {
    tag "${meta.id}"
    label "bsbolt"
    container "eod-tools.med-gen.ru/bsbolt:latest"

    input:
        tuple val(meta), path(alignment, arity : 1) 
        path reference_dir
    output:
        tuple val(meta), path("*.txt"), emit: statistics
        tuple val(meta), path("*.bedGraph"), emit: tracks 

    script:
    """
    bsbolt CallMethylation -I ${alignment} -DB ${reference_dir} -O ${meta.id} -verbose -t ${task.cpus} -BG    
    mv ${meta.id}.CGmap ${meta.id}.txt
    """
    stub:
    """
    touch ${meta.id}.txt
    touch ${meta.id}.bedGraph
    """
}