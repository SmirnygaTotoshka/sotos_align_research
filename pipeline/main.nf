/**
@author - Антон Смирнов
*/
nextflow.enable.dsl = 2
// Корректные названия функций из nf-schema v.2.6.1 для валидации samplesheets и параметров
include { paramsSummaryLog } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { samplesheetToList } from 'plugin/nf-schema'
include { paramsHelp } from 'plugin/nf-schema'

include { PROCESS_READS } from './subworkflows/reads_processing'
include { PREPARE_REFERENCE } from "./subworkflows/prepare_reference"
include { EXTRACT_BISULFITE as EXTRACT_BISULITE_FROM_MIX } from "./subworkflows/extract_bisulfite"
include { BSBOLT } from "./subworkflows/bsbolt"

outputDir = params.output


workflow {
    main: 
        log.info paramsSummaryLog(workflow)
        validateParameters()

        
        validate_samplesheet = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")).flatten()
        reads = PROCESS_READS(validate_samplesheet)
        path_to_ref = params.bisulfite_reference == null ? channel.empty() : channel.fromPath(params.bisulfite_reference)
        bisulfite_reference = PREPARE_REFERENCE(path_to_ref)

        reads_by_library_type = reads.processed_reads.branch{
            meta, _r -> 
                PURE: meta.type == "pure"
                MIXED: meta.type == "mixed"
        }

        extracted_bisulfite = EXTRACT_BISULITE_FROM_MIX(reads_by_library_type.MIXED)
        bisulfite_reads = reads_by_library_type.PURE.mix(extracted_bisulfite.reads)

        qc = reads.qc.mix(extracted_bisulfite.qc)
        fastp = reads.fastp.mix(extracted_bisulfite.fastp)
        meth_results = BSBOLT(bisulfite_reads, 
                              bisulfite_reference.bisulfite_reference, 
                              bisulfite_reference.cpg_context, 
                              bisulfite_reference.conversion_context)

        workflow.onComplete = {
            def cmd = ["chmod", "-R", "755", outputDir]
            def cmd2 = ["chown", "-R", "stotoshka:domain^users", outputDir]
            
            try {
                cmd.execute().waitFor()
                cmd2.execute().waitFor()
                log.info "Output directory permissions updated"
            } catch (Exception e) {
                log.error "Failed to update permissions: ${e.message}"
            }

            try{
                telegram_notify() 
                
            } catch (Exception e){
                log.error "Failed to notify: ${e.message}"
            }
        }
        
    publish:
        qc = qc
        fastp = fastp
        bisulfite_reference = bisulfite_reference.bisulfite_reference
        cpg_context = bisulfite_reference.cpg_context
        conversion_context = bisulfite_reference.conversion_context
        alignment = meth_results.alignment
        cgmap = meth_results.cgmap
        track = meth_results.track
        by_locus = meth_results.by_locus
}


output{
    qc{
        path "QC/fastqc"
        mode 'copy'
    }
    fastp{
        path "QC/fastp"
        mode 'copy'
    }
    bisulfite_reference{
        path  "reference" 
        enabled params.save_reference
        mode 'copy'
    }
    cpg_context{
        path "context"
        mode 'copy'
    }
    conversion_context{
        path "context"
        mode 'copy'
    }
    alignment{
        path "alignmnets"
        mode 'copy'
    }
    cgmap{
        path "methylation/cgmaps"
        mode 'copy'
    }
    track{
        path "methylation/tracks"
        mode 'copy'
    }
    by_locus{
        path "QC/statistics"
        mode 'copy'
    }
}

def telegram_notify(){
    def BOT_TOKEN = secrets.TELEGRAM_API_KEY
    def CHAT_ID = secrets.TELEGRAM_CHAT_ID

    def status_str = workflow.success ? 'OK' : 'failed'
    def error = !workflow.success ? "%0AError message is ${workflow.errorMessage}" : ""
    def run_path = "${outputDir}"

    def message = "Run name: ${workflow.runName} %0APath: ${run_path}%0APipeline completed at: ${workflow.complete}%0ADuration: ${workflow.duration}%0AExecution status: ${status_str}.${error}".replace(" ","+")
    def sout = new StringBuilder()
    def serr = new StringBuilder()
    def proc = """
        curl -s -X POST -H "Content-Type: application/txt" https://api.telegram.org/bot${BOT_TOKEN}/sendMessage \
            -d chat_id=-${CHAT_ID} \
            -d text="${message}"
    """.execute()
    proc.consumeProcessOutput(sout, serr)
    proc.waitForOrKill(600000)
}


