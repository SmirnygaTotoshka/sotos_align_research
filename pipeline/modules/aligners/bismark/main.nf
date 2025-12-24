process ALIGN {
    tag "${meta.id}"
    label 'bismark_align'
    
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"
    input:
        tuple val(meta), path(fastq, arity: 1..2) 
        path reference_folder
    output:
        tuple val(meta), path("${meta.id}_bismark_bt2.bam", arity : 1)
    
    script:
        if (meta.reverse == ""){
            """
                bismark ${task.ext.bismark_torrent_str} --basename ${meta.id}_bismark_bt2 ${reference_folder} ${fastq[0]}
            """
        }
        else{
            """
                bismark ${task.ext.bismark_illumina_str} --basename ${meta.id}_bismark_bt2 ${reference_folder} -1 ${fastq[0]} -2 ${fastq[1]}
                mv ${meta.id}_bismark_bt2_pe.bam ${meta.id}_bismark_bt2.bam
            """
        }
    stub:
    """
    touch ${meta.id}_bismark_bt2.bam
    """
}