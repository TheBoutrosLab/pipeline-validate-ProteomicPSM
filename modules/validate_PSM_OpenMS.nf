
include { index_peptides } from "${moduleDir}/index_peptides"
include { score_FDR_OpenMS } from "${moduleDir}/score_FDR_OpenMS"
include { merge_PSM } from "${moduleDir}/merge_PSM"
include { aggregate_consensus_id } from "${moduleDir}/aggregate_consensus_id"
include { score_FDR_OpenMS as score_FDR_OpenMS_consensus } from "${moduleDir}/score_FDR_OpenMS"
include { filter_PSM } from "${moduleDir}/filter_PSM"

workflow validate_PSM_OpenMS {
    take:
        idXML
    main:
        index_peptides(idXML, file(params.decoy_database))
        // Calculate FDR using OpenMS
        score_FDR_OpenMS(index_peptides.out[0])
        idXMLs = score_FDR_OpenMS.out[0].collect()
        // Merge idXML from different search engines
        merge_PSM(idXMLs)
        aggregate_consensus_id(merge_PSM.out[0])
        score_FDR_OpenMS_consensus(aggregate_consensus_id.out[0])
        filter_PSM(score_FDR_OpenMS_consensus.out[0])
}
