//
// Perform the exploratoy analysis using kraken
//

include { UNTAR as UNTAR_SCOUTING_DB              } from '../modules/nf-core/modules/untar/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_SCOUTING     } from '../modules/nf-core/modules/kraken2/kraken2/main'
include { PREPARE_KRAKEN_REPORT                   } from '../modules/local/prepare_kraken_report/main'
include { KRONA_KRONADB                           } from '../modules/nf-core/modules/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY                  } from '../modules/nf-core/modules/krona/ktimporttaxonomy/main'

workflow SEEK {
    take:
    reads     // channel: [ val(meta), path(reads) ]
    scout_db  // path: /path/to/scout/database


    main:
    ch_versions =  Channel.empty

    if (scout_db.endsWith("tar.gz") or scout_db.endsWith(".tgz")) {

        ch_krakendb_scout = [[], file(params.scout_db)]
        
        UNTAR_HOST_DB (ch_krakendb_scout)

        ch_scout_db      = UNTAR_SCOUTING_DB.out.untar.map{ it[1] }
        ch_versions      = ch_versions.mix(UNTAR_HOST_DB.out.versions)

    } else {
        ch_scout_db = file(scout_database)
    }

    KRAKEN2_SCOUTING (
        reads,
        ch_scout_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_SCOUTING.out.versions)

    KRONA_KRONADB ()

    PREPARE_KRAKEN_REPORT (
        KRAKEN2_SCOUTING.out.report
    )

    KRONA_KTIMPORTTAXONOMY (
        PREPARE_KRAKEN_REPORT.out.krona_report,
        KRONA_KRONADB.out.db
    )
}