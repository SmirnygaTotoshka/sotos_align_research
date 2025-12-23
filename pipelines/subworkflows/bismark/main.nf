include {AMPLICON_CLIP} from "./../../modules/samtools/amplicon_clip"
include {COVERAGE_STATS} from "./../../modules/samtools/coverage_stats"

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


process METH_CALL{
    tag "${meta.id}"
    label "bismark"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta),path(alignment, arity : 1) 

    output:
        tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: alignments
        tuple val(meta), path("*.gz"), emit: tracks 
        tuple val(meta), path("*.txt"), emit: statistics

    script:
        def mode = meta.reverse == "" ? "--single-end" : "--paired-end"
        """
            samtools sort -@ ${task.cpus} -n ${alignment} > ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam
            bismark_methylation_extractor ${task.ext.bismark_param_str} ${mode} ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam 
            samtools sort -@ ${task.cpus} ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam > ${meta.id}_bismark_bt2_ready.bam
            samtools index -@ ${task.cpus} ${meta.id}_bismark_bt2_ready.bam
            rm ${meta.id}_bismark_bt2_cleaned_sorted_byname.bam
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

process UNIFY_METH_FORMAT{
    tag "${meta.id}"
    label "bismark"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta),path(meth_data, arity : 1) 

    output:
        tuple val(meta), path("*.csv"), emit: tracks 

    script:
    """

    """
    stub:
    """
    
    """
}

workflow BISMARK{
    take:
        processed_reads
    main:
        raw_alignment = ALIGN(processed_reads, params.reference_folder)
        clipped_alignment = AMPLICON_CLIP(raw_alignment)
        meth_call_files = METH_CALL(clipped_alignment)
        coverage_stats = COVERAGE_STATS(clipped_alignment)
        meth_data = UNIFY_METH_FORMAT(meth_call_files)
    emit:
        raw_alignment = raw_alignment
        clipped_alignment = clipped_alignment
        meth_call_files = meth_call_files
        coverage_stats = coverage_stats
        meth_data = meth_data
}