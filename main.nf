nextflow.enable.dsl = 2

include { printPrelogue } from "${projectDir}/modules/common"

include { index_peptides; index_peptides as index_peptide_aggregated } from "${projectDir}/modules/index_peptides"
include { extract_PSM_features } from "${projectDir}/modules/extract_PSM_features"
include { score_PEP_percolator } from "${projectDir}/modules/score_PEP_percolator"
include { export_MSBooster_PIN } from "${projectDir}/modules/export_MSBooster_PIN"
include { run_MSBooster } from "${projectDir}/modules/run_MSBooster"
include { run_percolator_PIN } from "${projectDir}/modules/run_percolator_PIN"
include { import_percolator_POUT } from "${projectDir}/modules/import_percolator_POUT"
include { score_FDR_OpenMS } from "${projectDir}/modules/score_FDR_OpenMS"
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
    if (params.run_percolator && !(params.percolator_workflow in ['openms', 'msbooster'])) {
        error "percolator_workflow must be 'openms' or 'msbooster'."
    }
    if (params.run_percolator && params.percolator_workflow == 'msbooster' && !params.msbooster_params) {
        error "msbooster_params is required for percolator_workflow = 'msbooster'."
    }
    if (params.run_percolator && params.percolator_workflow == 'msbooster' && !params.input_mzml_csv) {
        error "input_mzml_csv is required for percolator_workflow = 'msbooster'."
    }
    idxml_rows = Channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .ifEmpty { error "Input CSV '${params.input_csv}' has no data rows." }
        .collect()
        .map { rows ->
            if (!rows[0].containsKey('sample_name')) {
                error "Input CSV '${params.input_csv}' is missing required sample_name column."
            }
            def seen = [] as Set
            rows.collect { row ->
                def sample_name = row.sample_name?.trim()
                def search_engine = row.search_engine?.trim()
                def idxml = row.idXML?.trim()
                if (!(sample_name ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                    error "Invalid or missing sample_name in input CSV: '${sample_name}'."
                }
                if (!search_engine || !idxml) {
                    error "Each input CSV row needs search_engine and idXML."
                }
                if (!seen.add([sample_name, search_engine.toLowerCase()])) {
                    error "Duplicate input CSV row for sample '${sample_name}' and search engine '${search_engine}'."
                }
                [sample_name, file(idxml), search_engine]
            }
        }
    ich = idxml_rows.flatMap { it }

    // PeptideIndexer can be skipped if the input idXML is from OpenMS search engine adapters.
    // Because PeptideIndexer is called internally already.
    if (params.skip_peptide_indexing) {
        idXML_indexed = ich
    } else {
        index_peptides(ich, file(params.decoy_database))
        idXML_indexed = index_peptides.out[0]
    }

    if (params.run_percolator == true) {
        if (params.percolator_workflow == 'msbooster') {
            // MSBooster supplies the Percolator features; PSMFeatureExtractor
            // cannot process idXML files containing multiple ProteinID runs.
            feature_idXML = idXML_indexed
            mzML_by_sample = Channel.fromPath(params.input_mzml_csv)
                .splitCsv(header:true)
                .ifEmpty { error "mzML CSV '${params.input_mzml_csv}' has no data rows." }
                .collect()
                .map { [rows: it] }
                .combine(idxml_rows.map { [rows: it] })
                .flatMap { mzml_payload, search_payload ->
                    def mzml_rows = mzml_payload.rows
                    def search_rows = search_payload.rows
                    def required_columns = ['sample_name', 'fraction', 'mzML']
                    def missing_columns = required_columns.findAll { !mzml_rows[0].containsKey(it) }
                    if (missing_columns) {
                        error "mzML CSV '${params.input_mzml_csv}' is missing columns: ${missing_columns.join(', ')}."
                    }
                    def raw_by_sample = [:].withDefault { [] }
                    def seen_fractions = [] as Set
                    def seen_mzml_paths = [] as Set
                    mzml_rows.each { row ->
                        def sample_name = row.sample_name?.trim()
                        def fraction = row.fraction?.trim()
                        def mzml = row.mzML?.trim()
                        if (!(sample_name ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                            error "Invalid or missing sample_name in mzML CSV: '${sample_name}'."
                        }
                        if (!fraction || !mzml) {
                            error "Each mzML CSV row needs fraction and mzML for sample '${sample_name}'."
                        }
                        if (!seen_fractions.add([sample_name, fraction])) {
                            error "Duplicate mzML CSV row for sample '${sample_name}' and fraction '${fraction}'."
                        }
                        if (!seen_mzml_paths.add([sample_name, mzml])) {
                            error "Duplicate mzML path for sample '${sample_name}': '${mzml}'."
                        }
                        raw_by_sample[sample_name] << [fraction, file(mzml)]
                    }
                    def search_samples = search_rows.collect { it[0] }.toSet()
                    def missing_samples = search_samples - raw_by_sample.keySet()
                    def extra_samples = raw_by_sample.keySet() - search_samples
                    if (missing_samples || extra_samples) {
                        error "mzML CSV sample names do not match input CSV (missing: ${missing_samples}; extra: ${extra_samples})."
                    }
                    raw_by_sample.collect { sample_name, fractions ->
                        def ordered = fractions.sort { left, right ->
                            def a = left[0].isInteger() ? left[0].toInteger() : null
                            def b = right[0].isInteger() ? right[0].toInteger() : null
                            (a != null && b != null) ? (a <=> b) : (left[0] <=> right[0])
                        }
                        [sample_name, ordered.collect { it[1] }]
                    }
                }
            paired_input = feature_idXML
                .map { sample_name, idxml, search_engine -> [sample_name, idxml, search_engine] }
                .combine(mzML_by_sample, by: 0)
            export_MSBooster_PIN(paired_input, file("${projectDir}/src/IdXMLPinBridge.py"))
            run_MSBooster(export_MSBooster_PIN.out[0], file(params.msbooster_params))
            run_percolator_PIN(run_MSBooster.out[0])
            import_percolator_POUT(run_percolator_PIN.out[0], file("${projectDir}/src/IdXMLPinBridge.py"))
            idXML_to_score_FDR = import_percolator_POUT.out[0]
        } else {
            // OpenMS PSMFeatureExtractor has no Sage-specific feature handler.
            indexed_by_engine = idXML_indexed.branch {
                sage: it[2].equalsIgnoreCase('Sage')
                supported: true
            }
            extract_PSM_features(indexed_by_engine.supported)
            feature_idXML = extract_PSM_features.out[0].mix(indexed_by_engine.sage)
            score_PEP_percolator(feature_idXML)
            idXML_to_score_FDR = score_PEP_percolator.out[0]
        }
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
                .flatMap { sample_name, idxml_files, search_engine ->
                    idxml_files.collect { idxml_file ->
                        [sample_name, idxml_file, search_engine]
                    }
                }

            fix_idXMLSpectraDataPath(ch_idXML_split, true, file("${projectDir}/src/SetIdXMLPrimaryMSRunPath.py"))

            ch_idXML_spectra_fixed = fix_idXMLSpectraDataPath.out[0]
                .map { sample_name, idXML, se ->
                    // Get the base file name with sample name and fraction ID.
                    // E.g., CCLE0100_fr1_MSFragger.idXML -> CCLE0100_fr1
                    file_origin = idXML.baseName.replaceFirst(/_${se}_FixSpectraDataPath$/, '')
                    return [sample_name, file_origin, idXML]
                }
                .groupTuple(by:[0,1])
                .map { sample_name, file_origin, idxml_files -> [sample_name, idxml_files] }
            merge_PSM(ch_idXML_spectra_fixed, false, true)
            aggregate_consensus_id(merge_PSM.out[0], true)
            merge_PSM_fractions(aggregate_consensus_id.out[0].groupTuple(), false, false)
            ch_merged_consensus_idXML = merge_PSM_fractions.out[0]
        } else {
            fix_idXMLSpectraDataPath(ch_idXML_individual, false, file("${projectDir}/src/SetIdXMLPrimaryMSRunPath.py"))
            // Merge search engine results and call consensus
            merge_PSM(fix_idXMLSpectraDataPath.out[0]
                .map { sample_name, idxml, search_engine -> [sample_name, idxml] }
                .groupTuple(), false, false)
            aggregate_consensus_id(merge_PSM.out[0], false)
            ch_merged_consensus_idXML = aggregate_consensus_id.out[0]
        }

        if (params.skip_peptide_indexing) {
            ch_indexed_consensus_idXML = ch_merged_consensus_idXML.map { sample_name, idxml -> [sample_name, idxml, 'consensus'] }
        } else {
            index_peptide_aggregated(ch_merged_consensus_idXML.map { sample_name, idxml -> [sample_name, idxml, 'consensus'] }, file(params.decoy_database))
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
