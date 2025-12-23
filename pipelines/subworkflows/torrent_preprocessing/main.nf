include {FASTQC as FASTQC_BEFORE} from "./../../modules/fastqc"
include {FASTQC as FASTQC_AFTER} from "./../../modules/fastqc"

process BAM_TO_FASTQ{
    tag "${meta.id}"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"
    containerOptions "${ workflow.containerEngine == 'singularity' ? "-B ${meta.forward}" : ""}"
    
    input:
        val meta
    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    """
        samtools fastq ${meta.forward} | gzip -c > ${meta.id}_raw.fastq.gz
    """
    stub:
    """
    touch ${meta.id}_raw.fastq.gz
    """
}

process TORRENT_TRIMMING {
    tag "$meta.id"
    label "process_medium"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path('*_R1.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')        , emit: json
    when:
        params.trim_fastq
    script:
        """
        cutadapt ${task.ext.args} \
        --json cutadapt.json \
        --cores ${task.cpus} \
        -o ${meta.id}_R1.fastq.gz ${reads[0]}
        """
    stub:
    """
        touch cutadapt.json
        touch ${meta.id}_R1.fastq.gz
    """
}

workflow TORRENT_PREPROCESSING {
    take:
        sample_data
    main:
        merged_fastq = sample_data | BAM_TO_FASTQ
        qc_before_trim = FASTQC_BEFORE(merged_fastq)
        trimming_result =  Channel.empty()
        qc_after_trim = Channel.empty()
        if (params.trim_fastq){
            trimming_result = TORRENT_TRIMMING(qc_before_trim.reads)
            qc_after_trim = FASTQC_AFTER(trimming_result.reads)
        }
        output = Channel.empty()
        if (params.trim_fastq){
            output = trimming_result.reads
        }
        else{
            output = merged_fastq
        }
    emit:
        output
        qc_before = qc_before_trim
        qc_after = qc_after_trim
        trim_report = trimming_result?.json // null-safe
}
