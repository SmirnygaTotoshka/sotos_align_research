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
include { BISMARK } from 'subworkflows/bismark'
include { BWA_METH } from 'subworkflows/bwa-meth'
include { BSBOLT } from 'subworkflows/bsbolt'


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

        aligned_reads = Channel.empty()

        if (params.aligner == "bismark") {
            aligned_reads = BISMARK(processed_reads)
        } else if (params.aligner == "bwa-meth") {
            aligned_reads = BWA_METH(processed_reads)
        } else if (params.aligner == "bsbolt") {
            aligned_reads = BSBOLT(processed_reads)
        }
        
    publish:
        qc_before = processed_reads.qc_before
        qc_after = processed_reads?.qc_after
        trim_report = processed_reads?.trimming_result // null-safe
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
