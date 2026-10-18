# Changelog

## [Unreleased]

### Added
- Optional MSBooster feature extraction before direct Percolator rescoring using the FragPipe image and bundled DIA-NN, with a separate sample/fraction/mzML CSV and idXML score import; the OpenMS PercolatorAdapter route remains the default ([#9](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/9)).
- Support multiple samples or plexes in one run by requiring `sample_name` in every input CSV row and keeping validation, consensus, and outputs separate by sample ([#8](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/issues/8)).

### Fixed
- Bypass OpenMS PSMFeatureExtractor in the MSBooster route, allowing DIA idXML inputs whose Percolator features are supplied by MSBooster.
- Use the rank-based ConsensusID algorithm when Percolator is disabled, because the default PEP-based algorithm cannot score q-value inputs. Retain PEPMatrix when Percolator is enabled.

### Tests
- Keep the OpenMS integration case, add a matching MSBooster Percolator case and an MSFragger DIA MSBooster case with empty input templates; remove the skip-indexing, no-Percolator, and unfractionated cases.

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
