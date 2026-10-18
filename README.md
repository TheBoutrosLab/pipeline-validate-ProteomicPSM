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

The left lane validates each search-engine idXML. Peptide indexing and Percolator rescoring are optional. Sage retains the features already in its idXML; Comet, MSFragger, and X!Tandem use PSMFeatureExtractor first. The selected Percolator workflow then rescores the PSMs. FalseDiscoveryRate, IDFilter, and MzTabExporter produce the filtered engine results.

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

The sample names in the two CSV files must match. Spectrum references in idXML must contain `scan=N` and identify MS2 scans in the corresponding mzML. When a scan number occurs in several fractions, `file_origin` or `id_merge_index` in the idXML must disambiguate it. Missing or ambiguous matches fail during PIN export.

The idXML files must preserve target/decoy annotations for FDR estimation. If peptide indexing is enabled, the FASTA must contain matching target and decoy protein sequences. Percolator also needs enough decoy PSMs to train; a small test dataset may require `run_percolator = false`.

### Config

| Parameter | Type | Description |
| --------- | ---- | ----------- |
| `input_csv` | string(path) | Path to the input CSV. |
| `input_mzml_csv` | string(path) | Path to the `sample_name,fraction,mzML` CSV, required for the MSBooster route. |
| `output_dir` | string(path) | Root output directory. |
| `decoy_database` | string(path) | Target-decoy protein FASTA used for peptide indexing. Required by the configuration schema. |
| `decoy_string` | string | Decoy protein prefix used by PeptideIndexer. Defaults to `DECOY_`. |
| `enzyme` | string | Digestion enzyme for PeptideIndexer, IDFilter, and Percolator. Defaults to `Trypsin`. |
| `skip_peptide_indexing` | bool | Skip PeptideIndexer when input idXML is already indexed. Defaults to `false`. |
| `run_percolator` | bool | Rescore PSMs with PercolatorAdapter. Defaults to `true`. |
| `percolator_workflow` | string | `openms` (default) runs PercolatorAdapter; `msbooster` exports PIN, adds predicted features with MSBooster, runs Percolator directly, and imports PEPs into idXML. |
| `msbooster_params` | string(path) | MSBooster prediction settings. Defaults to [input/msbooster.params](input/msbooster.params), which uses DIA-NN. The pipeline supplies input PIN, mzML, and DIA-NN executable paths. |
| `fragpipe_version` | string | FragPipe image version. Defaults to `24.0`. |
| `docker_image_fragpipe` | string | Image containing MSBooster, DIA-NN, and Percolator. Defaults to `fcyucn/fragpipe:24.0`, as used by the protein quantification pipeline. |
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

The MSBooster route uses the FragPipe image for MSBooster, DIA-NN prediction, and direct Percolator execution. It writes PEP, q-value, and SVM scores back to the original idXML PSMs. The original `openms` workflow remains available and is the default.

---

## Outputs

Each sample or plex has final results under `<output_dir>/validate-ProteomicPSM-<version>/<sample_name>/OpenMS-<version>/output/`:

| Output | Description |
| ------ | ----------- |
| `<sample_name>_<engine>_filter.idXML` | Search-engine identifications after PSM-level FDR estimation and IDFilter. |
| `<sample_name>_<engine>_filter.mzTab` | mzTab export of the filtered engine result. |
| `<sample_name>_<engine>_idswitch_pep.idXML` | Unfiltered identifications with Percolator PEP restored as the active score; written when Percolator is enabled. |
| `<sample_name>_consensus_filter.idXML` | Filtered consensus identifications; written when `merge_consensus = true`. |
| `<sample_name>_consensus_filter.mzTab` | mzTab export of the filtered consensus result. |

Intermediate idXML files are published under each sample’s `intermediate/` directory when `save_intermediate_files = true`. Pipeline reports are stored in a run-level timestamped `log-validate-ProteomicPSM-*` directory; process logs are grouped there by sample. [nftest.yml](nftest.yml) defines OpenMS and MSBooster multi-sample integration cases plus an MSFragger DIA case. Fill in the [DIA idXML CSV](test/test_integration_msbooster_msfragger_dia/input.csv) and [DIA mzML CSV](test/test_integration_msbooster_msfragger_dia/mzml.csv) before running the DIA case.

---

## References

1. Pfeuffer et al. OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. *Nature Methods*, 21, 365–367 (2024). [doi:10.1038/s41592-024-02197-7](https://doi.org/10.1038/s41592-024-02197-7)
2. Percolator: Käll et al. Semi-supervised learning for peptide identification from shotgun proteomics datasets. *Nature Methods*, 4, 923–925 (2007). [doi:10.1038/nmeth1113](https://doi.org/10.1038/nmeth1113)

---

## License

pipeline-validate-ProteomicPSM is licensed under the GNU General Public License version 2. See [LICENSE](LICENSE) for the terms of the GNU GPL license.

pipeline-validate-ProteomicPSM

Validate peptide-spectrum matches with OpenMS and Percolator.

Copyright (C) 2026 Sanford Burnham Prebys Medical Discovery Institute ("Boutros Lab")

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
