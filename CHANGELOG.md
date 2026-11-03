# Changelog

## [Unreleased]

### Fixed

- Match preserved pepXML native spectrum IDs directly to mzML MS2 spectra during MSBooster PIN export, including diaTracer `frame=N` IDs.

## [2.4.0] - 2026-11-02

### Added

- Use ProtLoom 0.3.1 PT/PC-split PSM FDR as an alternative to OpenMS FalseDiscoveryRate without protein inference, preserving source idXML metadata.

## [2.3.2] - 2026-10-27

### Fixed

- Export valid MSBooster calibration scores: log10 expectation values for MSFragger and Comet, and raw hyperscores for Sage (#31).

## [2.3.1] - 2026-10-27

### Added

- Publish MSBooster edited PIN files used by Percolator with sample, search-engine, and fraction identities for downstream QC (#28).
- Publish OpenMS PSMFeatureExtractor idXML files as standard final outputs.

## [2.3.0] - 2026-10-25

### Fixed

- Repair enzyme-specific pepXML cleavage termini and PEP-derived probabilities before iProphet and ProteinProphet (#25).

## [2.2.2] - 2026-10-23

### Added

- Support Epifany's greedy shared-peptide resolution and use its resolved peptide evidence for EasyPQP group assignment; report zero, one, and multiple-group PSM counts.

## [2.2.1] - 2026-10-23

### Changed

- Share the EasyPQP PSM q-value cutoff across all inference routes and give ProtLoom protein-group selection its own cutoff (#20).

## [2.2.0] - 2026-10-23

### Added

- Add selectable standalone Fido and OpenMS Epifany cohort protein grouping for non-ProtLoom EasyPQP export, retaining ProteinProphet as the default.

### Changed

- Use `ghcr.io/theboutroslab/fido:1.0` for standalone Fido protein inference.
- Rename the MSBooster, Percolator PIN, iProphet, PeptideProphet, ProtLoom, and protein-grouping processes to describe their actions; update the workflow schematic to show the ProteinProphet, Fido, Epifany, and ProtLoom routes.

## [2.1.0] - 2026-10-22

### Added

- Export EasyPQP PSM and peptide tables without ProtLoom by running cohort iProphet and ProteinProphet on Percolator results; retain OpenMS q-value selection, report excluded ambiguous mappings, and publish the Prophet intermediate files.

## [2.0.0] - 2026-10-19

### Changed

- Prefix cohort-level ProtLoom graph, manifest, and EasyPQP reports with `dataset_id`.
- Filter ProtLoom-annotated PSMs with the existing OpenMS `IDFilter` criteria after switching PEP or original scores to ProtLoom q-value; export mzTab and SHA-512 from the filtered result.
- Allow ProtLoom-annotated idXML to use ProtLoom PEP or q-value as the primary peptide and protein score while preserving the previous score in metadata.
- Write SHA-512 checksums beside each final idXML, including OpenMS consensus results, and beside the unfiltered ProtLoom-annotated idXML.
- Annotate Percolator idXML with all mapped ProtLoom PSM and protein-group FDR results without a q-value cutoff; allow unmatched PSMs to be kept or removed, and publish per-sample ProtLoom idXML files.
- Add an opt-in EasyPQP export switch and a dataset ID for cohort-level ProtLoom output paths.
- Export cohort-wide minimal EasyPQP `<dataset_id>-psm.tsv` and `<dataset_id>-peptides.tsv` from ProtLoom FDR and assignment outputs, omitting peptides assigned to multiple groups and reporting exclusion counts.
- Accept a ProtLoom project YAML for built database paths and sample-specific GTF annotations; generate task-local project manifests with current Percolator idXML files.
- Optional ProtLoom protein inference after MSBooster and Percolator: local inference for one sample, cohort graph and joint inference for multiple samples, followed by ProtLoom PSM and protein-group FDR; CCS is disabled.
- Optional MSBooster feature extraction before direct Percolator rescoring using the FragPipe image and bundled DIA-NN, with a separate sample/fraction/mzML CSV and idXML score import; the OpenMS PercolatorAdapter route remains the default ([#9](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/9)).
- Support multiple samples or plexes in one run by requiring `sample_name` in every input CSV row and keeping validation, consensus, and outputs separate by sample ([#8](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/8)).

### Fixed

- Preserve replaced idXML scores under their original score-type names without extra score-type or orientation metadata.
- Use pyOpenMS 3.5-compatible peptide-list and protein-group APIs when writing ProtLoom-annotated idXML.
- Skip per-engine and consensus OpenMS FDR calculation when ProtLoom is enabled; use ProtLoom q-values for final `IDFilter` and mzTab export.
- Load and validate the ProtLoom project YAML during configuration setup, then pass its resolved settings to the workflow.
- Use ProtLoom 0.2.3 and pass lower-is-better Percolator PEP idXML directly to inference; ProtLoom database preparation now uses `buildDB`.
- Bypass OpenMS PSMFeatureExtractor in the MSBooster route, allowing DIA idXML inputs whose Percolator features are supplied by MSBooster.
- Use the rank-based ConsensusID algorithm when Percolator is disabled, because the default PEP-based algorithm cannot score q-value inputs. Retain PEPMatrix when Percolator is enabled.

## [1.0.1] - 2026-10-11

### Fixed

- Restore Percolator PEP scores after FDR using the `Posterior Error Probability_score` metadata key in `IDScoreSwitcher` ([#4](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/4)).
- Use `PeptideIdentificationList` when loading idXML in the primary MS run path helper with pyOpenMS 3.5 ([#5](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/5)).

## [1.0.0] - 2026-10-10

### Added

- Nextflow DSL2 workflow that reads search-engine idXML paths from a sample-level CSV and optionally indexes peptide hits against a target-decoy FASTA.
- Optional Percolator rescoring followed by OpenMS PSM-level false discovery rate estimation, configurable PSM filtering, and mzTab export. Runs without Percolator are supported.
- Sage-specific routing that sends its feature-annotated idXML directly to PercolatorAdapter because PSMFeatureExtractor does not recognize Sage. Other supported search engines use PSMFeatureExtractor before rescoring.
- Optional consensus identification across search engines, including fraction-aware splitting, primary MS run path correction, merging, and consensus FDR filtering.
- Boutros Lab configuration, parameter schema, container profiles, resource allocation, versioned output directories, process logs, and workflow reports.
- NFTest integration cases for the standard workflow, skipped peptide indexing, disabled Percolator, and unfractionated consensus; test CSV paths support environment-variable expansion.

[1.0.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/releases/tag/v1.0.0
[1.0.1]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v1.0.0...v1.0.1
[2.0.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v1.0.1...v2.0.0
[2.1.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.0.0...v2.1.0
[2.2.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.1.0...v2.2.0
[2.2.1]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.2.0...v2.2.1
[2.2.2]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.2.1...v2.2.2
[2.3.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.2.2...v2.3.0
[2.3.1]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.3.0...v2.3.1
[2.3.2]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.3.1...v2.3.2
[2.4.0]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.3.2...v2.4.0
[unreleased]: https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/compare/v2.4.0...HEAD
