nextflow.enable.dsl = 2

include { printPrelogue } from "${projectDir}/modules/common"

include { index_peptides; index_peptides as index_peptide_aggregated } from "${projectDir}/modules/index_peptides"
include { extract_PSM_features } from "${projectDir}/modules/extract_PSM_features"
include { score_PEP_percolator } from "${projectDir}/modules/score_PEP_percolator"
include { export_MSBooster_PIN } from "${projectDir}/modules/export_MSBooster_PIN"
include { extract_PSMFeatures_MSBooster } from "${projectDir}/modules/extract_PSMFeatures_MSBooster"
include { score_PSM_PercolatorPIN } from "${projectDir}/modules/score_PSM_PercolatorPIN"
include { import_percolator_POUT } from "${projectDir}/modules/import_percolator_POUT"
include { score_FDR_OpenMS } from "${projectDir}/modules/score_FDR_OpenMS"
include { switch_score_2_pep } from "${projectDir}/modules/switch_score_2_pep"
include { fix_idXMLSpectraDataPath } from "${projectDir}/modules/fix_idXMLSpectraDataPath.nf"
include { merge_PSM; merge_PSM as merge_PSM_fractions } from "${projectDir}/modules/merge_PSM"
include { aggregate_consensus_id } from "${projectDir}/modules/aggregate_consensus_id"
include { score_FDR_OpenMS as score_FDR_OpenMS_consensus } from "${projectDir}/modules/score_FDR_OpenMS"
include { filter_PSM } from "${projectDir}/modules/filter_PSM"
include { filter_PSM as filter_PSM_consensus } from "${projectDir}/modules/filter_PSM"
include { filter_PSM as filter_PSM_protloom } from "${projectDir}/modules/filter_PSM"
include { split_idXML } from "${projectDir}/modules/split_idXML"
include { export_mzTab; export_mzTab as export_mzTab_consensus } from "${projectDir}/modules/export_mzTab"
include { infer_ProteinGroups_ProtLoomLocal } from "${projectDir}/modules/infer_ProteinGroups_ProtLoomLocal"
include { build_CohortGraph_ProtLoom } from "${projectDir}/modules/build_CohortGraph_ProtLoom"
include { infer_ProteinGroups_ProtLoomJoint } from "${projectDir}/modules/infer_ProteinGroups_ProtLoomJoint"
include { score_FDR_ProtLoom } from "${projectDir}/modules/score_FDR_ProtLoom"
include { annotate_ProtLoom_idXML } from "${projectDir}/modules/annotate_ProtLoom_idXML"
include { switch_ProtLoom_qvalue } from "${projectDir}/modules/switch_ProtLoom_qvalue"
include { export_EasyPQP_reports } from "${projectDir}/modules/export_EasyPQP_reports"
include { calculate_SHA512 } from "${projectDir}/modules/calculate_SHA512"
include { prepare_ProteinProphet_idXML } from "${projectDir}/modules/prepare_ProteinProphet_idXML"
include { extract_ProteinProphet_EnzymeRule } from "${projectDir}/modules/extract_ProteinProphet_EnzymeRule"
include { convert_ProteinProphet_pepXML } from "${projectDir}/modules/convert_ProteinProphet_pepXML"
include { repair_ProteinProphet_pepXML } from "${projectDir}/modules/repair_ProteinProphet_pepXML"
include { integrate_PSM_iProphet } from "${projectDir}/modules/integrate_PSM_iProphet"
include { infer_ProteinGroups_ProteinProphet } from "${projectDir}/modules/infer_ProteinGroups_ProteinProphet"
include { convert_ProteinProphet_idXML } from "${projectDir}/modules/convert_ProteinProphet_idXML"
include { export_Fido_graph } from "${projectDir}/modules/export_Fido_graph"
include { infer_ProteinGroups_Fido } from "${projectDir}/modules/infer_ProteinGroups_Fido"
include { parse_Fido_groups } from "${projectDir}/modules/parse_Fido_groups"
include { infer_ProteinGroups_Epifany } from "${projectDir}/modules/infer_ProteinGroups_Epifany"
include { export_ProteinGroup_EasyPQP } from "${projectDir}/modules/export_ProteinGroup_EasyPQP"

printPrelogue()

workflow {
    // Check the selected rescoring route before creating any input channels.
    if (params.run_percolator && !(params.percolator_workflow in ['openms', 'msbooster'])) {
        error "percolator_workflow must be 'openms' or 'msbooster'."
    }
    if (params.run_percolator && params.percolator_workflow == 'msbooster' && !params.msbooster_params) {
        error "msbooster_params is required for percolator_workflow = 'msbooster'."
    }
    if (params.run_percolator && params.percolator_workflow == 'msbooster' && !params.input_mzml_csv) {
        error "input_mzml_csv is required for percolator_workflow = 'msbooster'."
    }
    def protloom_project = params.protloom_project

    // Keep all search-engine rows together while checking sample names and
    // duplicate sample/engine pairs; downstream tasks use [sample, idXML, engine].
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

    // mzML lists are shared by all search engines for each sample. Both
    // MSBooster and the EasyPQP exports need these fraction-aware paths.
    if ((params.run_percolator && params.percolator_workflow == 'msbooster') ||
        params.export_easypqp_reports) {
        // One mzML list is shared by every search engine for a sample;
        // fraction labels determine its stable order.
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
    }

    // OpenMS search adapters may already have run PeptideIndexer.
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
            paired_input = feature_idXML
                .map { sample_name, idxml, search_engine -> [sample_name, idxml, search_engine] }
                .combine(mzML_by_sample, by: 0)
            // Convert idXML to PIN, add MSBooster features, rescore in
            // Percolator, then transfer the scores back to idXML.
            export_MSBooster_PIN(paired_input)
            extract_PSMFeatures_MSBooster(export_MSBooster_PIN.out[0], file(params.msbooster_params))
            score_PSM_PercolatorPIN(extract_PSMFeatures_MSBooster.out[0])
            import_percolator_POUT(score_PSM_PercolatorPIN.out[0])
            ch_rescored_idXML = import_percolator_POUT.out[0]
        } else {
            // The original OpenMS feature extraction and Percolator route.
            // OpenMS PSMFeatureExtractor has no Sage-specific feature handler.
            indexed_by_engine = idXML_indexed.branch {
                sage: it[2].equalsIgnoreCase('Sage')
                supported: true
            }
            extract_PSM_features(indexed_by_engine.supported)
            feature_idXML = extract_PSM_features.out[0].mix(indexed_by_engine.sage)
            score_PEP_percolator(feature_idXML)
            ch_rescored_idXML = score_PEP_percolator.out[0]
        }
    } else {
        ch_rescored_idXML = idXML_indexed
    }

    if (params.run_protloom) {
        // Keep Percolator PEP as the active score for ProtLoom inference.
        ch_idXML_individual = ch_rescored_idXML
    } else {
        // The OpenMS route computes and filters PSM q-values per engine.
        score_FDR_OpenMS(ch_rescored_idXML)
        filter_PSM(
            score_FDR_OpenMS.out[0],
            file(params.filter_psm_whitelist_proteins),
            file(params.filter_psm_whitelist_peptides),
            file(params.filter_psm_blacklist_proteins),
            file(params.filter_psm_blacklist_peptides),
            file(params.filter_psm_in_silico_digestion)
        )
        ch_final_mztab_idxml = filter_PSM.out[0]
        ch_checksum_idxml = filter_PSM.out[0]

        // Restore PEP in the unfiltered idXML for downstream inference;
        // the filtered OpenMS output retains q-values.
        if (params.run_percolator) {
            switch_score_2_pep(score_FDR_OpenMS.out[0])
            ch_idXML_individual = switch_score_2_pep.out[0]
        } else {
            ch_idXML_individual = score_FDR_OpenMS.out[0]
        }

        if (params.export_easypqp_reports) {
            // Group the unfiltered Percolator PEP results cohort-wide. Only
            // ProteinProphet needs the pepXML / iProphet conversion path.
            // EasyPQP PSM selection remains based on OpenMS q-values;
            // Epifany's output supplies associations when greedy resolution is enabled.
            if (params.protein_inference_algorithm == 'ProteinProphet') {
                prepare_ProteinProphet_idXML(ch_rescored_idXML, file(params.decoy_database))
                // OpenMS writes incorrect pepXML termini for non-Trypsin enzymes;
                // repair each converted file before cohort-wide iProphet.
                extract_ProteinProphet_EnzymeRule(prepare_ProteinProphet_idXML.out[0])
                convert_ProteinProphet_pepXML(extract_ProteinProphet_EnzymeRule.out[0])
                repair_ProteinProphet_pepXML(convert_ProteinProphet_pepXML.out[0])
                prophet_pepxmls = repair_ProteinProphet_pepXML.out[0]
                    .map { sample_name, pepxml, search_engine -> pepxml }
                    .collect()
                integrate_PSM_iProphet(prophet_pepxmls, file(params.decoy_database))
                infer_ProteinGroups_ProteinProphet(integrate_PSM_iProphet.out.pepxml, file(params.decoy_database))
                convert_ProteinProphet_idXML(infer_ProteinGroups_ProteinProphet.out.protxml)
                ch_protein_groups = convert_ProteinProphet_idXML.out.groups
            } else {
                cohort_pep_idxmls = ch_rescored_idXML
                    .map { sample_name, idxml, search_engine -> idxml }
                    .collect()
                if (params.protein_inference_algorithm == 'Fido') {
                    export_Fido_graph(cohort_pep_idxmls)
                    infer_ProteinGroups_Fido(export_Fido_graph.out.graph, export_Fido_graph.out.proteins)
                    parse_Fido_groups(infer_ProteinGroups_Fido.out.results, export_Fido_graph.out.mapping)
                    ch_protein_groups = parse_Fido_groups.out.groups
                } else {
                    infer_ProteinGroups_Epifany(cohort_pep_idxmls)
                    ch_protein_groups = infer_ProteinGroups_Epifany.out.groups
                }
            }

            group_easypqp_input = score_FDR_OpenMS.out[0]
                .join(mzML_by_sample, by: 0)
                .collect(flat: false)
                .map { rows ->
                    def ordered = rows.sort { a, b ->
                        def sample_order = a[0] <=> b[0]
                        sample_order != 0 ? sample_order : (a[2] <=> b[2])
                    }
                    def raw_files = ordered.collectMany { it[3] }.unique { it.toString() }
                    def staged = ordered.collect { it[1] } + raw_files
                    def names = staged.collect { it.name }
                    if (names.size() != names.toSet().size()) {
                        error 'Protein grouping EasyPQP inputs must have unique basenames for staging.'
                    }
                    def manifest = ordered.collect { row ->
                        [sample_name: row[0], idxml: row[1].name,
                         mzmls: row[3].collect { it.name }]
                    }
                    [groovy.json.JsonOutput.toJson(manifest),
                     ordered.collect { it[1] }, raw_files]
                }
            export_ProteinGroup_EasyPQP(group_easypqp_input, ch_protein_groups)
        }
    }

    if (params.merge_consensus) {
        // Consensus is scoped to one sample. Fractionated runs first combine
        // engines within each fraction, then combine the fraction results.
        if (params.fractionated) {
            // Split idXML by fraction before fixing each primary MS run path.
            split_idXML(ch_idXML_individual)
            ch_idXML_split = split_idXML.out[0]
                .flatMap { sample_name, idxml_files, search_engine ->
                    idxml_files.collect { idxml_file ->
                        [sample_name, idxml_file, search_engine]
                    }
                }

            fix_idXMLSpectraDataPath(ch_idXML_split, true)

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
            fix_idXMLSpectraDataPath(ch_idXML_individual, false)
            // Merge search engine results and call consensus
            merge_PSM(fix_idXMLSpectraDataPath.out[0]
                .map { sample_name, idxml, search_engine -> [sample_name, idxml] }
                .groupTuple(), false, false)
            aggregate_consensus_id(merge_PSM.out[0], false)
            ch_merged_consensus_idXML = aggregate_consensus_id.out[0]
        }

        if (!params.run_protloom) {
            if (params.skip_peptide_indexing) {
                ch_indexed_consensus_idXML = ch_merged_consensus_idXML.map { sample_name, idxml -> [sample_name, idxml, 'consensus'] }
            } else {
                index_peptide_aggregated(ch_merged_consensus_idXML.map { sample_name, idxml -> [sample_name, idxml, 'consensus'] }, file(params.decoy_database))
                ch_indexed_consensus_idXML = index_peptide_aggregated.out[0]
            }

            // Score and filter the merged consensus separately from engine results.
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
            ch_checksum_idxml = ch_checksum_idxml.mix(filter_PSM_consensus.out[0])
        }
    }

    if (params.run_protloom) {
        // ProtLoom uses the PEP-scored, unfiltered identifications. Its FDR
        // tables annotate PSMs and groups here; they do not impose a cutoff.
        // Use one identification set per sample. Multiple search engines require
        // the consensus result so that a spectrum is not counted several times.
        if (params.merge_consensus) {
            protloom_source = ch_merged_consensus_idXML
        } else {
            protloom_source = ch_idXML_individual.groupTuple(by: 0)
                .map { sample_name, idxmls, engines ->
                    if (idxmls.size() != 1) {
                        error "ProtLoom needs merge_consensus = true when sample '${sample_name}' has multiple search engines."
                    }
                    [sample_name, idxmls[0]]
                }
        }

        protloom_idxml = protloom_source.map { sample_name, idxml ->
            def annotation_id = protloom_project.sample_annotations[sample_name]
            if (!annotation_id) {
                error "ProtLoom project YAML has no annotation_id for input sample '${sample_name}'."
            }
            [sample_name, idxml, annotation_id]
        }
        protloom_samples = protloom_idxml.collect(flat: false)
            .map { rows -> rows.sort { a, b -> a[0] <=> b[0] } }
        // Database paths and sample-to-annotation mappings were resolved by
        // config/methods.config, so tasks only need staged files here.
        protloom_database = file(protloom_project.fasta)
        protloom_entries = file(protloom_project.entries)
        protloom_annotation_ids = protloom_project.annotations.collect { it.annotation_id }
        protloom_gtfs = protloom_project.annotations.collect { file(it.gtf_path) }
        protloom_markers = protloom_project.markers

        local_input = protloom_samples.filter { it.size() == 1 }.flatMap { it }
        infer_ProteinGroups_ProtLoomLocal(local_input, protloom_database, protloom_entries,
            protloom_annotation_ids, protloom_gtfs, protloom_markers)

        joint_input = protloom_samples.filter { it.size() > 1 }
            .map { rows -> [rows.collect { it[0] }, rows.collect { it[1] }] }
        // Multi-sample inference shares one cohort graph but produces a
        // separate assignment and FDR result for each sample.
        build_CohortGraph_ProtLoom(joint_input, protloom_database, protloom_entries,
            protloom_markers)
        joint_infer_input = build_CohortGraph_ProtLoom.out.graph
            .combine(protloom_idxml)
            .map { graph, sample_name, idxml, annotation_id ->
                [sample_name, idxml, graph, annotation_id]
            }
        infer_ProteinGroups_ProtLoomJoint(joint_infer_input, protloom_database, protloom_entries,
            protloom_annotation_ids, protloom_gtfs, protloom_markers)

        inference_results = infer_ProteinGroups_ProtLoomLocal.out.inference
            .mix(infer_ProteinGroups_ProtLoomJoint.out.inference)
        score_FDR_ProtLoom(inference_results, protloom_markers)

        // Join by sample_name: inference contributes idXML/assignments, while
        // computeFDR contributes PSM and protein-group q-value tables.
        protloom_annotation_input = inference_results
            .join(score_FDR_ProtLoom.out[0], by: 0)
            .map { sample_name, idxml, group_table, assignments, psm_fdr, group_fdr ->
                [sample_name, idxml, assignments, psm_fdr, group_fdr]
            }
        annotate_ProtLoom_idXML(protloom_annotation_input)
        // The published annotation may retain unmatched hits. Only matched
        // hits have a ProtLoom q-value and can enter the final q-value filter.
        protloom_mapped_idxml = annotate_ProtLoom_idXML.out.annotated
            .map { sample_name, annotated, mapped -> [sample_name, mapped] }
        if (params.protloom_idxml_primary_score in ['q_value', 'q-value']) {
            protloom_qvalue_idxml = protloom_mapped_idxml
        } else {
            switch_ProtLoom_qvalue(protloom_mapped_idxml)
            protloom_qvalue_idxml = switch_ProtLoom_qvalue.out.switched
        }
        filter_PSM_protloom(
            protloom_qvalue_idxml.map { sample_name, idxml -> [sample_name, idxml, 'ProtLoom'] },
            file(params.filter_psm_whitelist_proteins),
            file(params.filter_psm_whitelist_peptides),
            file(params.filter_psm_blacklist_proteins),
            file(params.filter_psm_blacklist_peptides),
            file(params.filter_psm_in_silico_digestion)
        )
        ch_final_mztab_idxml = filter_PSM_protloom.out[0]
        // Preserve a checksum for the full ProtLoom annotation as well as
        // the final IDFilter result.
        protloom_annotated_idxml = annotate_ProtLoom_idXML.out.annotated
            .map { sample_name, annotated, mapped -> [sample_name, annotated] }
        ch_checksum_idxml = ch_final_mztab_idxml.mix(protloom_annotated_idxml)

        if (params.export_easypqp_reports) {
            // EasyPQP reports are cohort-wide. Collect all sample inputs into
            // one task after joining their FDR tables and matching mzML lists.
            easypqp_input = inference_results
                .join(score_FDR_ProtLoom.out[0], by: 0)
                .join(mzML_by_sample, by: 0)
                .collect(flat: false)
                .map { rows ->
                    def sorted_rows = rows.sort { a, b -> a[0] <=> b[0] }
                    def staged_paths = sorted_rows.collectMany { row ->
                        [row[1], row[3], row[4], row[5]] + row[6]
                    }
                    def names = staged_paths.collect { it.name }
                    if (names.size() != names.toSet().size()) {
                        error 'ProtLoom EasyPQP inputs must have unique basenames for staging.'
                    }
                    def manifest = sorted_rows.collect { row ->
                        [sample_name: row[0], idxml: row[1].name,
                         assignments: row[3].name, psm_fdr: row[4].name,
                         group_fdr: row[5].name, mzmls: row[6].collect { it.name }]
                    }
                    [groovy.json.JsonOutput.toJson(manifest),
                     sorted_rows.collect { it[1] }, sorted_rows.collect { it[3] },
                     sorted_rows.collect { it[4] }, sorted_rows.collect { it[5] },
                     sorted_rows.collectMany { it[6] }]
                }
            export_EasyPQP_reports(easypqp_input)
        }
    }

    // Both routes export mzTab from their final idXML.
    export_mzTab(ch_final_mztab_idxml)
    // Checksum final idXMLs and the unfiltered ProtLoom annotation.
    calculate_SHA512(ch_checksum_idxml)
}
