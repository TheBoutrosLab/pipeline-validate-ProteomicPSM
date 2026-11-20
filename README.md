# pipeline-validate-ProteomicPSM

[![Prepare release](https://img.shields.io/badge/Action-Create%20New%20Release-blue)](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/actions/workflows/prepare-release.yaml)

- [pipeline-validate-ProteomicPSM](#pipeline-validate-proteomicpsm)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Pipeline Steps](#pipeline-steps)
  - [Inputs](#inputs)
    - [Input CSV](#input-csv)
    - [mzML CSV for MSBooster](#mzml-csv-for-msbooster)
    - [Config](#config)
  - [Outputs](#outputs)
  - [References](#references)
  - [License](#license)

## Overview

This Nextflow DSL2 pipeline validates peptide-spectrum matches (PSMs) in idXML search results. One run can process multiple label-free samples or isobaric plexes, each with results from one or more search engines. It can index peptides against a target-decoy FASTA, rescore PSMs with [Percolator](https://doi.org/10.1038/nmeth1113), estimate PSM-level false discovery rates with OpenMS, filter identifications, and export mzTab. It can also combine search engines with ConsensusID, including results from multiple fractions.

---

## How To Run

1. Create a `.config` using the [template](config/template.config). Set `input_csv`, `output_dir`, and `decoy_database`. Match `decoy_string` to the target-decoy FASTA and review the indexing, Percolator, filtering, and consensus settings.

2. Create an input CSV using the [template](input/input.csv). List the sample or plex name and idXML result for each search engine. Use absolute idXML paths. For MSBooster, also create an [mzML CSV](input/mzml.csv) with one row per fraction and set `input_mzml_csv`.

3. Initialize the Git submodules and run the pipeline with Nextflow, for example:

   ```bash
   export NXF_VER=24.10.6
   export NXF_SYNTAX_PARSER=v1
   nextflow run TheBoutrosLab/pipeline-validate-ProteomicPSM -r <version> -c config/my-run.config
   ```

   The supplied configuration template defaults to Docker. Use `-profile apptainer` or `-profile singularity` when appropriate for the execution environment. The `local` profile uses installed executables for development. The [`.env-template`](.env-template) records the Nextflow version used by the shared SBP configuration modules.

---

## Pipeline Steps

![Proteomic PSM validation workflow](img/flow-diagram.drawio.svg?raw=true)

The left lane validates each search-engine idXML. Peptide indexing and Percolator rescoring are optional. Sage retains the features already in its idXML; Comet, MSFragger, and X!Tandem use PSMFeatureExtractor first. The selected Percolator workflow then rescores the PSMs. OpenMS or ProtLoom estimates PSM FDR before IDFilter and MzTabExporter produce the filtered engine results.

The right lane runs when `merge_consensus = true`. It groups inputs by sample or plex, aligns primary MS run paths, combines search engines within each group with ConsensusID, and performs consensus FDR estimation, filtering, and mzTab export. With `fractionated = true`, the pipeline splits results by fraction before consensus and merges fraction-level consensus results afterward. ConsensusID uses `PEPMatrix` when Percolator provides posterior error probabilities and `ranks` when Percolator is disabled; `consensusid_algorithm` can override this choice.

---

## Inputs

### Input CSV

The CSV has one idXML file per sample or plex and search engine:

| Field | Required | Description |
| ----- | -------- | ----------- |
| `sample_name` | yes | Sample or plex name. Use letters, digits, dot, underscore, or hyphen; the first character must be a letter or digit. |
| `search_engine` | yes | Search engine name, such as `Comet`, `Sage`, `MSFragger`, or `XTandem`. |
| `idXML` | yes | Absolute path to that engine's idXML search result. |

```csv
sample_name,search_engine,idXML
SBP0001,Comet,/data/SBP0001_Comet.idXML
SBP0001,Sage,/data/SBP0001_Sage.idXML
SBP0002,MSFragger,/data/SBP0002_MSFragger.idXML
```

The `sample_name` column is required even when the CSV contains only one sample. Duplicate sample/search-engine pairs are rejected.

### mzML CSV for MSBooster

For MSBooster, provide a second CSV using `input_mzml_csv`. It lists each fraction once per sample; every search engine for that sample receives the same set of mzML files. Fraction labels must be unique within a sample. Numeric fractions are ordered numerically; other labels are ordered alphabetically.

```csv
sample_name,fraction,mzML
SBP0001,1,/data/SBP0001_f01.mzML
SBP0001,2,/data/SBP0001_f02.mzML
SBP0002,1,/data/SBP0002_f01.mzML
```

The sample names in the two CSV files must match. Spectrum references in idXML must contain `scan=N` or `frame=N` and identify MS2 spectra in the corresponding mzML. diaTracer mzML files use `frame=N` native IDs. When a scan number occurs in several fractions, `file_origin` or `id_merge_index` in the idXML must disambiguate it. Missing or ambiguous matches fail during PIN export.

With `percolator_workflow = 'msbooster'`, the edited PIN files sent to Percolator are always published under `<final_output_dir>/<sample_name>/MSBooster/output/engine_<search_engine>/` as `run_XX-fraction_<fraction>_edited.pin`. For example, sample `SBP0001`, engine `MSFragger`, and fraction `2` produce `SBP0001/MSBooster/output/engine_MSFragger/run_02-fraction_2_edited.pin`. `run_XX` is the one-based position of the fraction after the mzML CSV rows for that sample are sorted (numeric labels numerically, other labels alphabetically). Search-engine and fraction names are URL-encoded in these paths when they contain characters unsuitable for filenames; for example, a space becomes `+`. Fractions without exported PSMs have no PIN. These files are final outputs regardless of `save_intermediate_files`.

The idXML files must preserve target/decoy annotations for FDR estimation. If peptide indexing is enabled, the FASTA must contain matching target and decoy protein sequences. Percolator also needs enough decoy PSMs to train; a small test dataset may require `run_percolator = false`.

### Config

| Parameter | Type | Description |
| --------- | ---- | ----------- |
| `input_csv` | string(path) | Path to the input CSV. |
| `input_mzml_csv` | string(path) | Path to the `sample_name,fraction,mzML` CSV, required for MSBooster or EasyPQP export. |
| `output_dir` | string(path) | Root output directory. |
| `decoy_database` | string(path) | Target-decoy protein FASTA used for peptide indexing. Required by the configuration schema. |
| `decoy_string` | string | Decoy protein prefix used by PeptideIndexer. Defaults to `DECOY_`; must match `database.decoy_string` in the ProtLoom project YAML when ProtLoom is enabled. |
| `enzyme` | string | Digestion enzyme for PeptideIndexer, IDFilter, and Percolator. Defaults to `Trypsin`. |
| `skip_peptide_indexing` | bool | Skip PeptideIndexer when input idXML is already indexed. Defaults to `false`. |
| `run_percolator` | bool | Rescore PSMs with PercolatorAdapter. Defaults to `true`. |
| `percolator_workflow` | string | `openms` (default) runs PercolatorAdapter; `msbooster` exports PIN, adds predicted features with MSBooster, runs Percolator directly, and imports PEPs into idXML. |
| `msbooster_params` | string(path) | MSBooster prediction settings. Defaults to [input/msbooster.params](input/msbooster.params), which uses DIA-NN. The pipeline supplies input PIN, mzML, and DIA-NN executable paths. |
| `fragpipe_version` | string | FragPipe image version. Defaults to `24.0`. |
| `docker_image_fragpipe` | string | Image containing MSBooster, DIA-NN, and Percolator. Defaults to `fcyucn/fragpipe:24.0`, as used by the protein quantification pipeline. |
| `run_protloom` | bool | Run ProtLoom inference and FDR after MSBooster and Percolator. Defaults to `false`. |
| `protloom_split_fdr` | bool | Use ProtLoom PT/PC-split PSM FDR instead of OpenMS FalseDiscoveryRate, without protein inference. Defaults to `false`; incompatible with `run_protloom`. |
| `dataset_id` | string | Directory name for cohort-level inference and EasyPQP outputs. Defaults to `cohort`. |
| `export_easypqp_reports` | bool | Export cohort-level `<dataset_id>-psm.tsv`, `<dataset_id>-peptides.tsv`, and exclusion counts using ProtLoom or the selected protein inference backend. Requires Percolator and `input_mzml_csv`. Defaults to `false`. |
| `protein_inference_algorithm` | string | Protein grouping backend when ProtLoom is off: `ProteinProphet` (default), `Fido`, or `Epifany` (`Epiphany` is accepted as an alias). |
| `proteinprophet_min_probability` | number | Minimum ProteinProphet group probability for the non-ProtLoom EasyPQP route; defaults to `0.0`. |
| `fido_min_probability` | number | Minimum Fido group posterior probability for EasyPQP export; defaults to `0.01`. |
| `epifany_min_probability` | number | Minimum Epifany group posterior probability for EasyPQP export; defaults to `0.01`. |
| `epifany_greedy_group_resolution` | string | Epifany shared-peptide resolution: `none` (default), `remove_associations_only`, or `remove_proteins_wo_evidence`. |
| `docker_image_fido` | string | Standalone Fido image; defaults to `ghcr.io/theboutroslab/fido:1.0`. |
| `easypqp_psm_qvalue` | number | Maximum PSM q-value for EasyPQP export in every inference route. Defaults to `0.01`. |
| `protloom_easypqp_group_qvalue` | number | Maximum ProtLoom protein-group q-value for EasyPQP export. Defaults to `0.01`. |
| `protloom_idxml_primary_score` | string | Primary score in ProtLoom-annotated idXML: `original` (default), `pep`, or `q_value` (`q-value` is also accepted). |
| `protloom_idxml_unmatched_psms` | string | Handling for PSMs without mapped ProtLoom results: `remove` (default) or `keep`. |
| `docker_image_protloom` | string | ProtLoom container image. Defaults to `ghcr.io/theboutroslab/protloom:0.3.1`. |
| `protloom_project_yaml` | string(path) | ProtLoom project YAML with built database paths and decoy markers; required for inference or split FDR. Inference also requires `annotations` and `samples`. Relative paths use the YAML file's directory; `$VARIABLE` paths use the environment. |
| `protloom_pt_order`, `protloom_pc_order` | string | Comma-separated protein-tier and peptide-class priorities. Required when ProtLoom is enabled. |
| `protloom_enzyme` | string | ProtLoom cleavage rule name; defaults to `trypsin`. |
| `percolator_post_processing_tdc` | bool | Enable Percolator target-decoy competition post-processing. Defaults to `true`. |
| `merge_consensus` | bool | Combine search-engine results with ConsensusID. Defaults to `true`. |
| `fractionated` | bool | Split fractionated idXML before consensus merging. Defaults to `false`. |
| `filter_psm_score_peptide` | number | Maximum PSM q-value passed to `IDFilter -score:psm` when q-values are used. Defaults to `0.01`. The parameter name is retained for compatibility. |
| `filter_psm_remove_decoys` | bool | Remove decoy identifications from filtered outputs. Defaults to `true`. |
| `filter_psm_whitelist_proteins`, `filter_psm_whitelist_peptides` | string(path) | Optional FASTA or idXML files restricting retained proteins or peptides. |
| `filter_psm_blacklist_proteins`, `filter_psm_blacklist_peptides` | string(path) | Optional FASTA or idXML files excluding proteins or peptides. |
| `filter_psm_in_silico_digestion` | string(path) | Optional FASTA for IDFilter digestion-based filtering. |
| `save_intermediate_files` | bool | Publish intermediate idXML files. Defaults to `false`. |
| `ucla_cds` | bool | Compatibility flag required by the shared execution setup. Keep `false` outside the legacy UCLA-CDS environment. |
| `work_dir` | string(path) | Optional location for Nextflow intermediate files. |
| `base_resource_update` | namespace | Optional per-process CPU and memory multipliers applied to default allocations. |

Additional OpenMS options and their pipeline parameter names are listed in [docs/params.rst](docs/params.rst).

The MSBooster route uses the FragPipe image for MSBooster, DIA-NN prediction, and direct Percolator execution. It writes PEP, q-value, and SVM scores back to the original idXML PSMs. The original `openms` workflow remains available and is the default. PIN export converts MSFragger and Comet `expect` scores to `log10_evalue`, and Sage `ln(hyperscore)` scores to raw `hyperscore`. These scores select PSMs for MSBooster's RT calibration; `useRT` controls whether RT features are computed.

When `run_protloom = true`, set `protloom_project_yaml` to the project YAML used to build the search database. The configuration setup loads and validates that YAML before the workflow starts. Its `database.output_fasta` and `database.output_entries_tsv` must point to the existing `buildDB` outputs. Each `samples` row must have a `sample_id` matching `sample_name` in the input CSV and an `annotation_id` listed under `annotations`. The YAML's `id_files` values are ignored: the pipeline inserts the idXML files produced during this run into per-task project manifests. All referenced FASTA, entries TSV, and GTF files are staged into the ProtLoom tasks. Database decoy and contaminant markers come from the YAML; set the pipeline's `decoy_string` to the same decoy marker. See [input/protloom-project.example.yaml](input/protloom-project.example.yaml) for the required fields.

The pipeline passes unfiltered Percolator idXML with lower-is-better PEP scores directly to ProtLoom 0.3.1. One sample runs `localInfer --disable-ccs`. Multiple samples run `jointGraph` once, followed by `jointInfer --disable-ccs` for each sample. `computeFDR` runs after either inference path. Samples with multiple search engines use the consensus result; set `merge_consensus = true` for those inputs. The default container is `ghcr.io/theboutroslab/protloom:0.3.1`; set `docker_image_protloom` to override it.

For split FDR without inference, set `run_protloom = false` and `protloom_split_fdr = true`. Supply `protloom_project_yaml` with the built database FASTA, entries TSV, and decoy markers, plus `protloom_pt_order` and `protloom_pc_order`. ProtLoom scores each engine and optional consensus idXML using its current primary score (Percolator PEP when enabled), then the usual `IDFilter` cutoff applies to its q-values. This mode scores the first hit per identification, preserves the original idXML provenance, and publishes PSM FDR TSVs, q-value and filtered idXMLs with checksums, and mzTab under each sample's `ProtLoom-<version>/output/`.

---

## Outputs

With both ProtLoom options off, each sample or plex has final OpenMS results under `<output_dir>/validate-ProteomicPSM-<version>/<sample_name>/OpenMS-<version>/output/`:

| Output | Description |
| ------ | ----------- |
| `<sample_name>_<engine>_filter.idXML` | Search-engine identifications after PSM-level FDR estimation and IDFilter. |
| `<sample_name>_<engine>_filter.idXML.sha512` | SHA-512 checksum of the filtered search-engine idXML. |
| `<sample_name>_<engine>_filter.mzTab` | mzTab export of the filtered engine result. |
| `<sample_name>_<engine>_idswitch_pep.idXML` | Unfiltered identifications with Percolator PEP restored as the active score; written when Percolator is enabled. |
| `<sample_name>_consensus_filter.idXML` | Filtered consensus identifications; written when `merge_consensus = true`. |
| `<sample_name>_consensus_filter.idXML.sha512` | SHA-512 checksum of the filtered consensus idXML. |
| `<sample_name>_consensus_filter.mzTab` | mzTab export of the filtered consensus result. |

With `run_protloom = true`, OpenMS FDR calculation is skipped because ProtLoom supplies q-values. Each sample has `<sample_name>/ProtLoom-<version>/output/` with `<sample_name>_groups.tsv`, `<sample_name>_peptides.tsv`, `<sample_name>_assignments.tsv`, `<sample_name>_run.json`, `<sample_name>_psm_fdr.tsv`, `<sample_name>_protein_group_fdr.tsv`, `<sample_name>_protloom.idXML` (annotated, before final filtering) and its `.sha512` checksum, `<sample_name>_ProtLoom_filter.idXML` and its `.sha512` checksum, and `<sample_name>_ProtLoom_filter.mzTab`. The final idXML is filtered with the same `IDFilter` settings as the OpenMS route, including `filter_psm_score_peptide` (default `0.01`). mzTab uses that final idXML; SHA-512 checksums cover both the annotated and final idXML files. Multi-sample runs additionally publish `<dataset_id>/ProtLoom-<version>/output/<dataset_id>-cohort_graph.json` and `<dataset_id>-cohort_manifest.json`. Set `dataset_id` to the desired cohort directory name (default `cohort`).

Each `.sha512` file contains the SHA-512 digest and idXML filename in standard checksum format. It is generated from the bytes of the corresponding idXML file.

Set `export_easypqp_reports = true` to publish cohort-wide `<dataset_id>-psm.tsv`, `<dataset_id>-peptides.tsv`, and `<dataset_id>-easypqp_export_summary.tsv` under `<dataset_id>/ProtLoom-<version>/output/` for EasyPQP library generation; export is disabled by default. The reports include target PSMs with q-values at or below `easypqp_psm_qvalue` and protein groups with q-values at or below `protloom_easypqp_group_qvalue` (both default to `0.01`). Since EasyPQP maps proteins by unmodified peptide sequence, the exporter excludes sequences assigned to different ProtLoom groups across the cohort and reports the excluded sequence and PSM counts in the summary. `Gene` uses the leading protein accession as a stable placeholder because the ProtLoom assignment table has no gene field. Supply the matching idXML and mzML files when running EasyPQP; the TSVs alone are filters and protein annotations, not spectra.

When `run_protloom = false`, select `protein_inference_algorithm = 'ProteinProphet'`, `'Fido'`, or `'Epifany'`. ProteinProphet (the default) runs iProphet on the Percolator pepXML files and publishes pepXML, protXML, and idXML. Fido converts the unfiltered Percolator PEP idXMLs into Fido’s native peptide–protein graph and target/decoy list, runs standalone Fido, and publishes its text output plus a parsed group JSON file. Epifany consumes the PEP idXMLs directly and publishes its cohort idXML using the current OpenMS image. The three EasyPQP TSVs and native grouping outputs are published under `<dataset_id>/<algorithm>/output/`. Fido and Epifany export inferred protein groups above their respective posterior probability cutoffs (`0.01` by default); ProteinProphet uses `proteinprophet_min_probability`. The PSM table uses the selected FDR route's per-engine q-value idXML at or below `easypqp_psm_qvalue`; the peptide table assigns each sequence to one accepted group from the selected algorithm. Decoys, PSMs without a unique accepted group, conflicting spectrum matches, and peptide sequences assigned to multiple groups are excluded and counted in the summary. `Gene` and `Protein ID` use the group's leading accession. `Mapped Proteins` lists other target accessions recorded in the retained PSMs' peptide evidence, separated by semicolons; it does not infer mappings from every group member. This route requires `run_percolator = true`, `input_mzml_csv`, and q-values (`fdr_algorithm_no_qvalues = false`); it works with either Percolator workflow.

Before iProphet, the ProteinProphet route reads each idXML's digestion-enzyme rule and repairs OpenMS-converted pepXML. It recalculates `num_tol_term` separately for the primary and alternative protein matches and makes the PEP-derived probability available at every NTT value. One-residue C-terminal and N-terminal cleavage rules are supported, as is `no cleavage` (only protein boundaries count as termini); other rules fail explicitly. Matches with unknown flanking residues retain their existing NTT and are counted in the repair log. The separate `num_missed_cleavages` field remains as written by OpenMS.

For Epifany, set `epifany_greedy_group_resolution = 'remove_associations_only'` to resolve shared peptide associations without removing protein hits, or `remove_proteins_wo_evidence` to also remove proteins left without evidence. The pipeline runs Epifany with `-algorithm:keep_best_PSM_only false` so nonbest PSMs remain in its output for EasyPQP matching; Epifany still uses only the best PSM per peptide for protein inference. In every Epifany mode, EasyPQP assigns passing PSMs to accepted groups using peptide evidence in Epifany's output; q-values and spectrum mappings still come from the per-engine q-value idXML files. Evidence matching resolves each idXML's mzML origin before comparing spectrum reference, modified sequence, charge, RT, and m/z; Epifany's merged `id_merge_index` is not compared directly with a per-engine index. The export summary reports PSMs matching zero, one, or multiple accepted groups, PSMs with one group after resolution (before spectrum deduplication and sequence ambiguity filtering), and PSMs without a matching Epifany identification. With the default `none`, Epifany leaves shared associations unresolved, so those PSMs remain excluded when they match multiple accepted groups. Epifany can omit or alter a peptide identification during inference; an unmatched q-value PSM has no resolved association and is excluded from EasyPQP export. Greedy resolution chooses an association from inferred protein probabilities; it does not establish that a shared peptide is biologically unique.

The annotation step does not apply a q-value cutoff: every input PSM with a matching ProtLoom PSM and protein-group result is annotated, regardless of its q-value. PSMs without a mapped result are removed from the annotated file by default; set `protloom_idxml_unmatched_psms = 'keep'` to retain them there with their original scores and protein evidence. The final `IDFilter` input always contains only mapped PSMs, since unmatched hits have no ProtLoom q-value to filter. When keeping unmatched PSMs in the annotated file, the writer places matched and unmatched hits in separate identification runs so their score types remain accurate. Matched peptide hits carry `protein_tier`, `peptide_class`, `q_value`, `probability`, `protloom_pep` (1 − probability), and `protloom_group_id` metadata. Their protein evidence points to the ProtLoom leading protein. Protein hits carry group q-value, probability, and PEP metadata, and `ProteinIdentification` records each group's full protein membership. By default, matched peptide hits keep their original Percolator score as the primary score in the annotated file and matched protein hits use the ProtLoom group q-value. Set `protloom_idxml_primary_score = 'pep'` or `'q_value'` to use that ProtLoom score as the primary score for matched peptide and protein hits. When available, the previous hit score is saved under its original score-type name in hit metadata. Before final filtering, `IDScoreSwitcher` changes `pep` or `original` primary PSM scores to ProtLoom q-value; `q_value` goes directly to `IDFilter`.

Intermediate idXML files are published under each sample’s `intermediate/` directory when `save_intermediate_files = true`. Pipeline reports are stored in a run-level timestamped `log-validate-ProteomicPSM-*` directory; process logs are grouped there by sample. [nftest.yml](nftest.yml) defines OpenMS and MSBooster multi-sample integration cases plus an MSFragger DIA case. Fill in the [DIA idXML CSV](test/test_integration_msbooster_msfragger_dia/input.csv) and [DIA mzML CSV](test/test_integration_msbooster_msfragger_dia/mzml.csv) before running the DIA case. From the repository root, run `python3 test/run-nftest.py test-integration-OpenMS`; set `NFT_CLUSTER=pines` for Apptainer, or `NFT_CLUSTER=boutroslab_aws` inside an AWS Slurm job. `NFT_RUN_ROOT` selects the output base directory.

---

## References

1. Pfeuffer et al. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. *Nature Methods*, 21, 365–367 (2024). [doi:10.1038/s41592-024-02197-7](https://doi.org/10.1038/s41592-024-02197-7)
2. Percolator: Käll et al. Semi-supervised learning for peptide identification from shotgun proteomics datasets. *Nature Methods*, 4, 923–925 (2007). [doi:10.1038/nmeth1113](https://doi.org/10.1038/nmeth1113)

---

## License

The code in pipeline-validate-ProteomicPSM is licensed under the
BSD 3-Clause License (`BSD-3-Clause`). See [LICENSE](LICENSE) for the terms.

Shared code and third-party tools retain their respective licenses.

Copyright (C) 2026 Sanford Burnham Prebys Medical Discovery Institute ("Boutros Lab")
