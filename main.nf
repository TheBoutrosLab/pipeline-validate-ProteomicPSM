nextflow.enable.dsl = 2

include { printPrelogue } from "${projectDir}/modules/common"

include { index_peptides; index_peptides as index_peptide_aggregated } from "${projectDir}/modules/index_peptides"
include { extract_PSM_features } from "${projectDir}/modules/extract_PSM_features"
include { score_PEP_percolator } from "${projectDir}/modules/score_PEP_percolator"
include { score_FDR_OpenMS } from "${projectDir}/modules/score_FDR_OpenMS"
include { switch_score_2_search_engine } from "${projectDir}/modules/switch_score_2_search_engine"
include { switch_score_2_pep } from "${projectDir}/modules/switch_score_2_pep"
include { fix_idXMLSpectraDataPath } from "${projectDir}/modules/fix_idXMLSpectraDataPath.nf"
include { merge_PSM; merge_PSM as merge_PSM_fractions } from "${projectDir}/modules/merge_PSM"
include { aggregate_consensus_id } from "${projectDir}/modules/aggregate_consensus_id"
include { score_FDR_OpenMS as score_FDR_OpenMS_consensus } from "${projectDir}/modules/score_FDR_OpenMS"
include { filter_PSM } from "${projectDir}/modules/filter_PSM"
include { filter_PSM as filter_PSM_consensus } from "${projectDir}/modules/filter_PSM"
include { split_idXML } from "${projectDir}/modules/split_idXML"
include { export_mzTab; export_mzTab as export_mzTab_consensus } from "${projectDir}/modules/export_mzTab"

printPrelogue()

workflow {
    ich = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.idXML), row.search_engine] }

    // PeptideIndexer can be skipped if the input idXML is from OpenMS search engine adapters.
    // Because PeptideIndexer is called internally already.
    if (params.skip_peptide_indexing) {
        idXML_indexed = ich
    } else {
        index_peptides(ich, file(params.decoy_database))
        idXML_indexed = index_peptides.out[0]
    }

    if (params.run_percolator == true) {
        // OpenMS PSMFeatureExtractor has no Sage-specific feature handler.
        indexed_by_engine = idXML_indexed.branch {
            sage: it[1].equalsIgnoreCase('Sage')
            supported: true
        }
        extract_PSM_features(indexed_by_engine.supported)
        score_PEP_percolator(extract_PSM_features.out[0].mix(indexed_by_engine.sage))

        percolator_by_engine = score_PEP_percolator.out[0].branch {
            sage: it[1].equalsIgnoreCase('Sage')
            supported: true
        }
        switch_score_2_search_engine(percolator_by_engine.supported)
        // FalseDiscoveryRate can score Sage PSMs directly from their Percolator PEPs.
        idXML_to_score_FDR = switch_score_2_search_engine.out[0].mix(percolator_by_engine.sage)
    } else {
        idXML_to_score_FDR = idXML_indexed
    }

    score_FDR_OpenMS(idXML_to_score_FDR)

    // Filter FDR for individual search engine results
    filter_PSM(
        score_FDR_OpenMS.out[0],
        file(params.filter_psm_whitelist_proteins),
        file(params.filter_psm_whitelist_peptides),
        file(params.filter_psm_blacklist_proteins),
        file(params.filter_psm_blacklist_peptides),
        file(params.filter_psm_in_silico_digestion)
    )

    export_mzTab(filter_PSM.out[0])

    // Swtich score of the unfiltered idXML back to PEP so it can be used to run iProphet &
    // ProteinProphet in quantitate-ProteomicProtein pipeline.
    if (params.run_percolator == true) {
        switch_score_2_pep(score_FDR_OpenMS.out[0])
        ch_idXML_individual = switch_score_2_pep.out[0]
    } else {
        ch_idXML_individual = score_FDR_OpenMS.out[0]
    }

    if (params.merge_consensus) {
        if (params.fractionated) {
            // Split idXML to fraction specific
            split_idXML(ch_idXML_individual)
            ch_idXML_split = split_idXML.out[0]
                .map { it ->
                    res = []
                    search_engine = it[1]
                    it[0].each { idXML_file ->
                        res.add( [
                            'search_engine': search_engine,
                            'idXML': idXML_file
                        ] )
                    }
                    return res
                }
                .flatten()
                .map { [it['idXML'], it['search_engine']] }

            fix_idXMLSpectraDataPath(ch_idXML_split, true, file("${projectDir}/src/SetIdXMLPrimaryMSRunPath.py"))

            ch_idXML_spectra_fixed = fix_idXMLSpectraDataPath.out[0]
                .map { idXML, se ->
                    // Get the base file name with sample name and fraction ID.
                    // E.g., CCLE0100_fr1_MSFragger.idXML -> CCLE0100_fr1
                    file_origin = idXML.baseName.split("_${se}")[0]
                    return [file_origin, idXML]
                }
                .groupTuple(by:0)
                .map { it[1] }
            merge_PSM(ch_idXML_spectra_fixed, false, true)
            aggregate_consensus_id(merge_PSM.out[0], true)
            merge_PSM_fractions(aggregate_consensus_id.out[0].collect(), false, false)
            ch_merged_consensus_idXML = merge_PSM_fractions.out[0]
        } else {
            fix_idXMLSpectraDataPath(ch_idXML_individual, false, file("${projectDir}/src/SetIdXMLPrimaryMSRunPath.py"))
            // Merge search engine results and call consensus
            merge_PSM(fix_idXMLSpectraDataPath.out[0].map{it[0]}.collect(), false, false)
            aggregate_consensus_id(merge_PSM.out[0], false)
            ch_merged_consensus_idXML = aggregate_consensus_id.out[0]
        }

        if (params.skip_peptide_indexing) {
            ch_indexed_consensus_idXML = ch_merged_consensus_idXML.map{ [it, 'consensus'] }
        } else {
            index_peptide_aggregated(ch_merged_consensus_idXML.map{ [it, 'consensus'] }, file(params.decoy_database))
            ch_indexed_consensus_idXML = index_peptide_aggregated.out[0]
        }

        // Rescore and filter FDR with the consensus ID
        score_FDR_OpenMS_consensus(ch_indexed_consensus_idXML)
        filter_PSM_consensus(
            score_FDR_OpenMS_consensus.out[0],
            file(params.filter_psm_whitelist_proteins),
            file(params.filter_psm_whitelist_peptides),
            file(params.filter_psm_blacklist_proteins),
            file(params.filter_psm_blacklist_peptides),
            file(params.filter_psm_in_silico_digestion)
        )
        export_mzTab_consensus(filter_PSM_consensus.out[0])
    }
}
