/**
@author - Антон Смирнов
*/
nextflow.enable.dsl = 2
// Корректные названия функций из nf-schema v.2.6.1 для валидации samplesheets и параметров
include { paramsSummaryLog } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'

include { PROCESS_READS } from 'subworkflows/reads_processing'
include { PREPARE_REFERENCE } from "subworkflows/prepare_reference"
include { EXTRACT_BISULITE as EXTRACT_BISULITE_FROM_MIX } from "subworkflows/extract_bisulfite"
include { BSBOLT } from "subworkflows/bsbolt"


workflow {
    main: 
        log.info paramsSummaryLog(workflow)
        validateParameters()
        
        validate_samplesheet = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        
        reads = PROCESS_READS(validate_samplesheet)
        bisulfite_reference = PREPARE_REFERENCE(params.bisulfite_reference)

        reads_by_library_type = reads.branch{
            PURE: meta.type == "pure"
            MIXED: meta.type == "mixed"
        }

        extracted_bisulfite = EXTRACT_BISULITE_FROM_MIX(reads_by_library_type.MIXED)
        bisulfite_reads = reads_by_library_type.mix(extracted_bisulfite.reads)
        meth_results = BSBOLT(bisulfite_reads, bisulfite_reference)
        
    publish:
        qc = reads.qc
        fastp = reads.fastp
        extracted_qc = extracted_bisulfite.qc
        extracted_fastp = extracted_fastp.fastp
        reference = bisulfite_reference.bisulfite_reference
        cpg_context = bisulfite_reference.cpg
        conversion_context = bisulfite_reference.conversion
}
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
