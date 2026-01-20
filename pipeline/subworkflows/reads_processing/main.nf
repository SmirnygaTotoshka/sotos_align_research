process COMBINE_FASTQ {
    tag "${meta.id}"

    input:
        val meta

    output:
        tuple val(meta), path("${meta.id}_R?_raw.fastq.gz", arity: 2)
    
    script:
    def forward = meta.forward.join(' ')
    def reverse = meta.reverse.join(' ')
    """
        zcat ${forward} | gzip -c > ${meta.id}_R1_raw.fastq.gz
        zcat ${reverse} | gzip -c > ${meta.id}_R2_raw.fastq.gz
    """
    stub:
    """
    touch ${meta.id}_R1_raw.fastq.gz
    touch ${meta.id}_R2_raw.fastq.gz
    """
}

process BAM_TO_FASTQ{
    tag "${meta.id}"
    
    input:
        val meta
    output:
        tuple val(meta), path("*.fastq.gz", arity: 1)

    script:
    def forward = meta.forward.join(' ')
    """
        samtools fastq ${forward} | gzip -c > ${meta.id}_raw.fastq.gz
    """
    stub:
    """
    touch ${meta.id}_raw.fastq.gz
    """
}

process FASTP {
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('*_processed.fastq.gz', arity: 1..2), emit: reads
        tuple val(meta), path('*.{json,html}'), emit: report

    script:
    def io_options = meta.reverse.size() == 0 ? "--in1 ${reads[0]} --out1 ${meta.id}_R1_processed.fastq.gz" : "--in1 ${reads[0]} --out1 ${meta.id}_R1_processed.fastq.gz --in2 ${reads[1]} --out2 ${meta.id}_R2_processed.fastq.gz"      
        """
        fastp ${io_options} \
        --adapter_fasta ${params.primers_fasta} \
        --cut_front ${params.cut_front} \
        --cut_tail ${params.cut_tail} \
        --n_base_limit ${params.n_base_limit} \
        --length_required ${params.read_min_len} \
        --thread ${task.cpus}
        """

    stub:
    if (meta.reverse.size() == 0){
        """
        touch ${meta.id}_R1_processed.fastq.gz
        touch "${meta.id}.fastp.json"
        touch "${meta.id}.fastp.html"
        touch "${meta.id}.fastp.log"
        """ 
    }
    else {
        """
        touch ${meta.id}_R1_processed.fastq.gz
        touch ${meta.id}_R2_processed.fastq.gz
        touch "${meta.id}.fastp.json"
        touch "${meta.id}.fastp.html"
        touch "${meta.id}.fastp.log"
        """
    }
}

process BBDUK_REPAIR {
    tag "${meta.id}"
    
    input:
        tuple val(meta), path(reads, arity: 1..2)

    output:
        tuple val(meta), path('*.fastq.gz')
    
    script:
    def io_options = meta.reverse.size() == 0 ? "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz" : "--in ${reads[0]} --out ${meta.id}_R1.fastq.gz --in2 ${reads[1]} --out2 ${meta.id}_R2.fastq.gz"      
    """
        bash /bbmap/repair.sh ${io_options} 
    """
    stub:
    if (meta.reverse.size() == 0){
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

process FASTQC {
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.{zip,html}')

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

workflow PROCESS_READS{
    take:
        samples_info
    main:
        split_by_tech = samples_info.branch{
            s ->
                ILLUMINA: s.reverse.size() != 0
                TORRENT: s.reverse.size() == 0
        }

        illumina_raw = COMBINE_FASTQ(split_by_tech.ILLUMINA)
        torrent_raw = BAM_TO_FASTQ(split_by_tech.TORRENT)

        raw_reads = illumina_raw.mix(torrent_raw) 
        trimmed_reads = FASTP(raw_reads) 
        processed_reads = BBDUK_REPAIR(trimmed_reads.reads)
        qc = FASTQC(processed_reads.mix(raw_reads))

    emit:
        processed_reads = processed_reads
        qc = qc
        fastp = trimmed_reads.report
}

workflow PROCESS_EXTRACTED_READS{
    take:
        extracted_reads
    main:
        trimmed_reads = FASTP(extracted_reads) 
        processed_reads = BBDUK_REPAIR(trimmed_reads.reads)
        qc = FASTQC(processed_reads.mix(extracted_reads))

    emit:
        processed_reads = processed_reads
        qc = qc
        fastp = trimmed_reads.report
}