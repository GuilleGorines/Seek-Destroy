//
// Remove host using kraken
//

include { UNTAR as UNTAR_HOST_DB } from '../../modules/'
include { KRAKEN2 as KRAKEN2_HOST_REMOVAL } from '../../modules/'


workflow DESTROY {
    take:
    reads   // channel: [ val(meta), path(reads) ]
    host_db // path : /path/to/host/database

    main:
    ch_versions =  Channel.empty()

    if (host_database.endsWith("tar.gz") or params.host_database.endsWith(".tgz")) {

        ch_krakendb_host = [[], file(params.host_database)]
        
        UNTAR_HOST_DB (ch_krakendb_host)

        ch_host_database = UNTAR_SCOUTING_DB.out.untar.map{ it[1] }
        ch_versions      = ch_versions.mix(UNTAR_HOST_DB.out.versions)

    } else {
        ch_host_database = Channel.fromPath(host_db)
    }

    KRAKEN2_HOST_REMOVAL (
        reads
        ch_host_database
        true,
        false
    )



    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}