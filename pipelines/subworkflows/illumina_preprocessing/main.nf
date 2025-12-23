include {FASTQC as FASTQC_BEFORE} from "./../../modules/fastqc"
include {FASTQC as FASTQC_AFTER} from "./../../modules/fastqc"


process COMBINE_FASTQ {
    tag "${meta.id}"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        val meta

    output:
        tuple val(meta), path("${meta.id}_R?_raw.fastq.gz")
    
    script:

        """
            forward_files=\$(ls ${meta.forward})
            reverse_files=\$(ls ${meta.reverse})
            zcat \${forward_files} | gzip -c > ${meta.id}_R1_raw.fastq.gz
            zcat \${reverse_files} | gzip -c > ${meta.id}_R2_raw.fastq.gz
        """

    stub:
    """
    touch ${meta.id}_R1_raw.fastq.gz
    touch ${meta.id}_R2_raw.fastq.gz
    """

}

process ILLUMINA_TRIMMING {
    tag "$meta.id"
    label 'process_medium'

    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
    tuple val(meta), path(reads)
    path output_dir
    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    when:
    params.trim_fastq
    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def adapter_list = task.ext.adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
        def fail_fastq = task.ext.save_trimmed_fail ? "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        def out_fq1 = task.ext.discard_trimmed_pass ?: "--out1 ${prefix}_1.fastp.fastq.gz" 
        def out_fq2 = task.ext.discard_trimmed_pass ?: "--out2 ${prefix}_2.fastp.fastq.gz"
        def merge_fastq = task.ext.save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
       
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            $out_fq1 \\
            $out_fq2 \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args

        """
    stub:
        def prefix              = task.ext.prefix ?: "${meta.id}"
        def touch_reads         = (discard_trimmed_pass) ? "" : (is_single_output) ? "echo '' | gzip > ${prefix}.fastp.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastp.fastq.gz ; echo '' | gzip > ${prefix}_2.fastp.fastq.gz"
        def touch_merged        = (!is_single_output && save_merged) ? "echo '' | gzip >  ${prefix}.merged.fastq.gz" : ""
        def touch_fail_fastq    = (!save_trimmed_fail) ? "" : meta.single_end ? "echo '' | gzip > ${prefix}.fail.fastq.gz" : "echo '' | gzip > ${prefix}.paired.fail.fastq.gz ; echo '' | gzip > ${prefix}_1.fail.fastq.gz ; echo '' | gzip > ${prefix}_2.fail.fastq.gz"
        """
        $touch_reads
        $touch_fail_fastq
        $touch_merged
        touch "${prefix}.fastp.json"
        touch "${prefix}.fastp.html"
        touch "${prefix}.fastp.log"
        """
}

process REPAIR {
    tag "${meta.id}"
    label 'process_medium'

    container "eod-tools.med-gen.ru/nf-lis-worker:latest"
    
    input:
        tuple val(meta), path(reads, arity: 2)

    output:
        tuple val(meta), path('*_repaired.fastq.gz') , emit: reads
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
        """
        [ ! -f  ${prefix}_R1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_R1.fastq.gz
        [ ! -f  ${prefix}_R2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_R2.fastq.gz

            bash /bbmap/repair.sh \
            --in=${prefix}_R1.fastq.gz \
            --in2=${prefix}_R2.fastq.gz \
            --out=${prefix}_R1_repaired.fastq.gz \
            --out2=${prefix}_R2_repaired.fastq.gz \
            $args

        """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_R1_repaired.fastq.gz"
    touch "${prefix}_R2_repaired.fastq.gz"

    """
}

workflow ILLUMINA_PREPROCESSING {
    take:
        sample_data
    main:
        merged_fastq = sample_data | COMBINE_FASTQ
        qc_before_trim = FASTQC_BEFORE(merged_fastq, params.output + "/QC/pre_fastqc")
        ready = Channel.empty()
        qc_after_trim = Channel.empty()
        trimming_result = Channel.empty()
        if (params.trim_fastq){
            trimming_result = ILLUMINA_TRIMMING(qc_before_trim.reads, params.output + "/QC/fastp")
            ready = REPAIR(trimming_result.reads)
            qc_after_trim = FASTQC_AFTER(ready.reads, params.output + "/QC/post_fastqc")
        }
        else {
            ready = REPAIR(merged_fastq)
        }
    emit:
        processed_reads = ready.reads
        qc_before = qc_before_trim
        qc_after = qc_after_trim
        trim_report = trimming_result?.html // null-safe
}

output{
    qc_before {
        path "${params.output}/QC/preprocessing"
    }
    qc_after {
        path "${params.output}/QC/postprocessing"
    }
    trim_report {
        path "${params.output}/QC/trimming"
    }
}
