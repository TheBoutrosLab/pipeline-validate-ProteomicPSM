# Changelog

## [1.0.0] - Pending release

### Added
- Nextflow DSL2 workflow that reads search-engine idXML paths from a sample-level CSV and optionally indexes peptide hits against a target-decoy FASTA.
- Optional Percolator rescoring followed by OpenMS PSM-level false discovery rate estimation, configurable PSM filtering, and mzTab export. Runs without Percolator are supported.
- Sage-specific routing that sends its feature-annotated idXML directly to PercolatorAdapter because PSMFeatureExtractor does not recognize Sage. Other supported search engines use PSMFeatureExtractor before rescoring.
- Optional consensus identification across search engines, including fraction-aware splitting, primary MS run path correction, merging, and consensus FDR filtering.
- Boutros Lab configuration, parameter schema, container profiles, resource allocation, versioned output directories, process logs, and workflow reports.
- NFTest integration cases for the standard workflow, skipped peptide indexing, disabled Percolator, and unfractionated consensus; test CSV paths support environment-variable expansion.
