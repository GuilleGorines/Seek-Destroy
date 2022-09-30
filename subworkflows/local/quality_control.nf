//
// Perform the quality control
//

include { CUTADAPT                                } from '../modules/nf-core/modules/cutadapt/main'
include { FASTP                                   } from '../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQC_POSTERIOR              } from '../modules/nf-core/modules/fastqc/main'