/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSeekanddestroy.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.multiqc_config
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.skip_scouting && !params.scout_database) { exit 1. 'No scouting database was chosen!'}
if (params.skip_host_removal && !params.host_database) { exit 1. 'No host database was chosen!'}

if (params.scout_database && !params.skip_scouting) { ch_scout_database = Channel.fromPath(params.scout_database) }
if (params.host_database && !params.skip_host_removal) { ch_host_database = Channel.fromPath(params.host_database) }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULE: Installed directly from nf-core/modules
include { FASTQC as FASTQC_PREVIOUS               } from '../modules/nf-core/modules/fastqc/main'
include { CUTADAPT                                } from '../modules/nf-core/modules/cutadapt/main'
include { FASTP                                   } from '../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQC_POSTERIOR              } from '../modules/nf-core/modules/fastqc/main'
include { UNTAR as UNTAR_SCOUTING_DB              } from '../modules/nf-core/modules/untar/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_SCOUTING     } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { PREPARE_KRAKEN_REPORT                   } from '../modules/local/prepare_kraken_report/main'
include { KRONA_KRONADB                           } from '../modules/nf-core/modules/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY                  } from '../modules/nf-core/modules/krona/ktimporttaxonomy/main'
include { UNTAR as UNTAR_HOST_DB                  } from '../modules/nf-core/modules/untar/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_HOST_REMOVAL } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { MULTIQC                                 } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS             } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SEEK_AND_DESTROY {

    ch_versions = Channel.empty()

    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // MODULE: Run FastQC: check initial quality
    FASTQC_PREVIOUS (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_PREVIOUS.out.versions)
    
    // QUALITY CONTROL SUBWORKFLOW??
    // MODULE: Run CutAdapt: remove adapters
    CUTADAPT (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // MODULE: Run FastP: trim sequences
    FASTP (
        CUTADAPT.out.reads,
        params.save_trimmed_reads,
        params.save_merged_reads
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // MODULE: Run FastQC: check quality after quality control

    FASTQC_POSTERIOR (
        FASTP.out.reads
    )

    // MODULE: Run Kraken2: exploration with the first database
    // Should this be a subworkflow?

    if (params.scout_database.endsWith("tar.gz") or params.scout_database.endsWith(".tgz")){
        ch_krakendb_scout = [[], file(params.scout_database)]
        UNTAR_SCOUTING_DB (ch_krakendb_scout)
        ch_scout_database = UNTAR_SCOUTING_DB.out.untar.map{ it[1] }

        ch_scout_database_test = UNTAR_SCOUTING_DB.out.untar.map{ it[1] }
        ch_scout_database_test.view()
    }
    
    KRAKEN2_SCOUTING (
        FASTP.out.reads,
        ch_scout_database,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_SCOUTING.out.versions)
    
    PREPARE_KRAKEN_REPORT (
        KRAKEN2_SCOUTING.out.report
    )
    ch_versions = ch_versions.mix(PREPARE_KRAKEN_REPORT.out.versions)

    // MODULE: Run Krona: visualization of the kraken2 scouting results
    KRONA_KRONADB ()
    ch_versions = ch_versions.mix(KRONA_KRONADB.out.versions)

    KRONA_KTIMPORTTAXONOMY (
        PREPARE_KRAKEN_REPORT.out.krona_report,
        KRONA_KRONADB.out.db
    )

    // if estipulated not to extract data, do not do this part
    // MODULE: Run Kraken2: to remove host reads
    // HOST REMOVAL SUBWORKFLOW??
    
    if (params.host_database.endsWith("tar.gz") or params.host_database.endsWith(".tgz")) {
        ch_krakendb_host = [[], file(params.host_database)]
        UNTAR_HOST_DB (ch_krakendb_host)
        ch_host_database = UNTAR_HOST_DB.out.untar.map{ it[1] }

        ch_host_database_test = UNTAR_HOST_DB.out.untar.map{ it[1] }
        ch_host_database_test.view()
    }
    
    
    KRAKEN2_HOST_REMOVAL (
        FASTP.out.reads,
        ch_host_database,
        true,
        false
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: MultiQC
    workflow_summary    = WorkflowSeekanddestroy.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PREVIOUS.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
