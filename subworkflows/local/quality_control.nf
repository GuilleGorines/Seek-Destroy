//
// Perform the quality control
//

include { CUTADAPT                                } from '../modules/nf-core/modules/cutadapt/main'
include { FASTP                                   } from '../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQC_POSTERIOR              } from '../modules/nf-core/modules/fastqc/main'

workflow QUALITY_CONTROL {
    take:
    reads                // channel: [ val(meta), path(reads) ]
    save_trimmed_reads   // value: boolean
    save_merged_reads    // value: boolean

    main:
    ch_versions =  Channel.empty

    CUTADAPT (
        reads
    )

    FASTP (
        CUTADAPT.out.reads
        save_trimmed_reads
        save_merged_reads
    )

    FASTQC_POSTERIOR (
        FASTP.out.reads
    )

    emit: 

}   