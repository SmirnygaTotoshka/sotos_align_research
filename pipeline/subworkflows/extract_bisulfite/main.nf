include { PROCESS_EXTRACTED_READS } from './../reads_processing'

process ALIGN {
    tag "${meta.id}"
    label 'align'
    container "eod-tools.med-gen.ru/sotos-align:1.0"
    input:
        tuple val(meta), path(reads)
        path reference
        path index
    output:
        tuple val(meta), path("*.bam"), path("*.bam.bai")

    script:
    def samtools_threads = Math.max(1, task.cpus.intdiv(3))
    def reads_files = reads.collect().join(' ')
    """
        bwa mem -t ${task.cpus} \
        -M -R "@RG\\tID:${meta.id}\\tSM:${meta.id}" \
        ${reference} ${reads_files} | \
        samtools view -b -h -@ ${samtools_threads} - | \
        samtools sort - -O BAM -o ${meta.id}_mixed.bam -@ ${samtools_threads}

        samtools index -M -@ ${task.cpus} ${meta.id}_mixed.bam
    """
    stub:
    """
    touch ${meta.id}_mixed.bam
    touch ${meta.id}_mixed.bam.bai
    """
}


process EXTRACT{
    tag "${meta.id}"
    label 'samtools'
container "eod-tools.med-gen.ru/sotos-align:1.0"
    input:
        tuple val(meta), path(total_alignment), path(total_index)
        path usual_panel
    output:
        tuple val(meta), path("*.fastq.gz", arity: 2)

    script:
    """

    samtools view -b -h -@ ${task.cpus} -L ${usual_panel} -U ${meta.id}_mismap.bam -o ${meta.id}_usual_raw.bam ${total_alignment}
    samtools sort -n -@ ${task.cpus} ${meta.id}_mismap.bam > ${meta.id}_mismap_sorted.bam
    samtools fixmate -@ ${task.cpus} ${meta.id}_mismap_sorted.bam - | \
    samtools view -f 0x2 -o ${meta.id}_mismap_fix.bam
    samtools fastq -1 ${meta.id}_R1.fastq -2 ${meta.id}_R2.fastq ${meta.id}_mismap_fix.bam
    gzip -c ${meta.id}_R1.fastq > ${meta.id}_R1.fastq.gz
    gzip -c ${meta.id}_R2.fastq > ${meta.id}_R2.fastq.gz
    """
    stub:
    """
    touch ${meta.id}_R1.fastq.gz
    touch ${meta.id}_R2.fastq.gz
    """
}

workflow EXTRACT_BISULFITE{
    take:
        input_reads
    main:
        index = channel.fromPath(params.usual_reference)
            .map { fasta -> file("${fasta}.{amb,ann,pac,bwt,sa,fai}")
            }.collect()
        alignments = ALIGN(input_reads, params.usual_reference, index)
        bisulfite_reads = EXTRACT(alignments, params.usual_bed)
        processed_reads = PROCESS_EXTRACTED_READS(bisulfite_reads)
    emit:
        reads = processed_reads.processed_reads
        qc = processed_reads.qc
        fastp = processed_reads.fastp
}