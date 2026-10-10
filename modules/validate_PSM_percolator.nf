
include { index_peptides } from "${moduleDir}/index_peptides"
include { score_FDR_OpenMS } from "${moduleDir}/score_FDR_OpenMS"
include { extract_PSM_features } from "${moduleDir}/extract_PSM_features"
include { score_PEP_percolator } from "${moduleDir}/score_PEP_percolator"
include { aggregate_consensus_id } from "${moduleDir}/aggregate_consensus_id"
include { switch_score } from "${moduleDir}/switch_score"
include { filter_PSM } from "${moduleDir}/filter_PSM"

workflow validate_PSM_percolator {
    take:
        idXML
    main:
        index_peptides(idXML, file(params.decoy_database))
        // Calculate Posterior Error Probability using percolator
        extract_PSM_features(index_peptides.out[0])
        score_PEP_percolator(extract_PSM_features.out[0])

        idXMLs = score_PEP_percolator.out[0].collect()
        aggregate_consensus_id(idXMLs)
        score_FDR_OpenMS(aggregate_consensus_id.out[0])
        filter_PSM(score_FDR_OpenMS.out[0])
}
