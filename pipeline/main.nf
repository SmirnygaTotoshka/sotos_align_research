/**
@author - Антон Смирнов
*/
nextflow.enable.dsl = 2
// Корректные названия функций из nf-schema v.2.6.1 для валидации samplesheets и параметров
include { paramsSummaryLog } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'

include { ILLUMINA_PREPROCESSING } from 'subworkflows/illumina_preprocessing'
include { TORRENT_PREPROCESSING } from 'subworkflows/torrent_preprocessing'
include { BAM_PROCESS } from "subworkflows/bam_process"
include { COVERAGE_STATS } from "subworkflows/coverage_stats"

include { ALIGN as BISMARK } from 'modules/aligners/bismark'
include { ALIGN as BWA_METH } from 'modules/aligners/bwa-meth'
include { ALIGN as BSBOLT } from 'modules/aligners/bsbolt'

include { METHYLATION_EXTRACTOR as BISMARK_CALL } from 'modules/aligners/bismark'
include { ALIGN as BWA_METH } from 'modules/aligners/bwa-meth'
include { ALIGN as BSBOLT } from 'modules/aligners/bsbolt'



process FORM_METH_MAP{
    tag "${meta.id}"
    label "python"
    container "eod-tools.med-gen.ru/nf-lis-worker:latest"

    input:
        tuple val(meta), path(meth_data, arity : 1) 
        path cpg_map

    output:
        tuple val(meta), path("*.csv"), emit: map 

    script:
    """
    python ${workflow.projectDir}/scripts/form_panel_methylation_panel.py --input ${meth_data} --output ${meta.id}_map.csv
    """
    stub:
    """
    touch ${meta.id}_map.csv
    """
}

workflow {
    main: 
        log.info paramsSummaryLog(workflow)
        validateParameters()
        
        validate_samplesheet = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        raw_reads = validate_samplesheet.branch{
            ILLUMINA: validate_samplesheet.reverse != ""
            TORRENT: validate_samplesheet.reverse == ""
        }

        illumina_processed = ILLUMINA_PREPROCESSING(raw_reads.ILLUMINA)
        torrent_preprocessing = TORRENT_PREPROCESSING(raw_reads.TORRENT)
        processed_reads = illumina_processed.mix(torrent_preprocessing)

        raw_alignment = Channel.empty()

        if (params.aligner == "bismark") {
            raw_alignment = BISMARK(processed_reads, params.reference_folder)
        } else if (params.aligner == "bwa-meth") {
            raw_alignment = BWA_METH(processed_reads, params.reference_folder)
        } else if (params.aligner == "bsbolt") {
            raw_alignment = BSBOLT(processed_reads, params.reference_folder)
        }

        processed_alignment = BAM_PROCESS(raw_alignment, params.panel_bed)

        meth_call_files = Channel.empty()
        if (params.aligner == "bismark") {
            raw_alignment = BISMARK_CALL(processed_alignment)
        } else if (params.aligner == "bwa-meth") {
            raw_alignment = BWA_METH_CALL(processed_alignment)
        } else if (params.aligner == "bsbolt") {
            raw_alignment = BSBOLT_CALL(processed_alignment)
        }

        meth_map = UNIFY_METH_FORMAT(meth_call_files)
        
    publish:
        qc_before = processed_reads.qc_before
        qc_after = processed_reads?.qc_after
        trim_report = processed_reads?.trimming_result // null-safe
        raw_alignment = raw_alignment
        processed_alignment = processed_alignment
        meth_call_files = meth_call_files
        coverage_stats = coverage_stats
        meth_map = meth_map
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
    coverage_stats {
        path "${params.output}/QC/coverage"
    }
    raw_alignment{
        path "${params.output}/intermediate/raw_alignment"
    }
    meth_call_files {
        path "${params.output}/intermediate/meth_calling"
    }
    processed_alignment {
        path "${params.output}/alignments"
    }
    meth_map{
        path "${params.output}/methylation"
    }
}
