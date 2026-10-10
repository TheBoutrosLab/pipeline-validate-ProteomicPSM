# pipeline-validate-ProteomicPSM

[![Prepare release](https://img.shields.io/badge/Action-Create%20New%20Release-blue)](https://github.com/TheBoutrosLab/pipeline-validate-ProteomicPSM/actions/workflows/prepare-release.yaml)

- [pipeline-validate-ProteomicPSM](#pipeline-validate-proteomicpsm)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Pipeline Steps](#pipeline-steps)
  - [Inputs](#inputs)
    - [Input CSV](#input-csv)
    - [Config](#config)
  - [Outputs](#outputs)
  - [References](#references)
  - [License](#license)

## Overview

This Nextflow DSL2 pipeline validates peptide-spectrum matches (PSMs) in idXML search results. One run processes one label-free sample or one isobaric plex, with results from one or more search engines. It can index peptides against a target-decoy FASTA, rescore PSMs with [Percolator](https://doi.org/10.1038/nmeth1113), estimate PSM-level false discovery rates with OpenMS, filter identifications, and export mzTab. It can also combine search engines with ConsensusID, including results from multiple fractions.

---

## How To Run

1. Create a `.config` using the [template](config/template.config). Set `input_csv`, `output_dir`, `sample_name`, and `decoy_database`. Match `decoy_string` to the target-decoy FASTA and review the indexing, Percolator, filtering, and consensus settings.

2. Create an input CSV using the [template](input/input.csv). List the idXML result for each search engine used on the same sample or plex. Use absolute idXML paths.

3. Initialize the Git submodules and run the pipeline with Nextflow, for example:

   ```bash
   export NXF_VER=24.10.6
   export NXF_SYNTAX_PARSER=v1
   nextflow run TheBoutrosLab/pipeline-validate-ProteomicPSM -r <version> -c config/my-run.config
   ```

   The supplied configuration template defaults to Docker. Use `-profile apptainer` or `-profile singularity` when appropriate for the execution environment. The `local` profile uses installed executables for development. The [`.env-template`](.env-template) records the Nextflow version used by the shared SBP configuration modules.

> :warning: Run one sample or plex at a time. With consensus enabled, all idXML rows in the CSV are combined into one sample-level result.

---

## Pipeline Steps

![Proteomic PSM validation workflow](img/flow-diagram.drawio.svg?raw=true)

The left lane validates each search-engine idXML. Peptide indexing and Percolator rescoring are optional. Sage goes directly to PercolatorAdapter with the features already in its idXML; Comet, MSFragger, and X!Tandem use PSMFeatureExtractor first. FalseDiscoveryRate, IDFilter, and MzTabExporter produce the filtered engine results.

The right lane runs when `merge_consensus = true`. It aligns primary MS run paths, combines search engines with ConsensusID, and performs consensus FDR estimation, filtering, and mzTab export. With `fractionated = true`, the pipeline splits results by fraction before consensus and merges fraction-level consensus results afterward.

---

## Inputs

### Input CSV

The CSV has one idXML file per search engine, with all rows belonging to the same sample or plex:

| Field | Required | Description |
| ----- | -------- | ----------- |
| `search_engine` | yes | Search engine name, such as `Comet`, `Sage`, `MSFragger`, or `XTandem`. |
| `idXML` | yes | Absolute path to that engine's idXML search result. |

```csv
search_engine,idXML
Comet,/data/SBP0001_Comet.idXML
Sage,/data/SBP0001_Sage.idXML
MSFragger,/data/SBP0001_MSFragger.idXML
```

The idXML files must preserve target/decoy annotations for FDR estimation. If peptide indexing is enabled, the FASTA must contain matching target and decoy protein sequences. Percolator also needs enough decoy PSMs to train; a small test dataset may require `run_percolator = false`.

### Config

| Parameter | Type | Description |
| --------- | ---- | ----------- |
| `input_csv` | string(path) | Path to the input CSV. |
| `output_dir` | string(path) | Root output directory. |
| `sample_name` | string | Name of the one sample or plex processed in this run. |
| `decoy_database` | string(path) | Target-decoy protein FASTA used for peptide indexing. Required by the configuration schema. |
| `decoy_string` | string | Decoy protein prefix used by PeptideIndexer. Defaults to `DECOY_`. |
| `enzyme` | string | Digestion enzyme for PeptideIndexer, IDFilter, and Percolator. Defaults to `Trypsin`. |
| `skip_peptide_indexing` | bool | Skip PeptideIndexer when input idXML is already indexed. Defaults to `false`. |
| `run_percolator` | bool | Rescore PSMs with PercolatorAdapter. Defaults to `true`. |
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

---

## Outputs

Final results are saved under `<output_dir>/validate-ProteomicPSM-<version>/<sample_name>/OpenMS-<version>/output/`:

| Output | Description |
| ------ | ----------- |
| `<sample_name>_<engine>_filter.idXML` | Search-engine identifications after PSM-level FDR estimation and IDFilter. |
| `<sample_name>_<engine>_filter.mzTab` | mzTab export of the filtered engine result. |
| `<sample_name>_<engine>_idswitch_pep.idXML` | Unfiltered identifications with Percolator PEP restored as the active score; written when Percolator is enabled. |
| `<sample_name>_consensus_filter.idXML` | Filtered consensus identifications; written when `merge_consensus = true`. |
| `<sample_name>_consensus_filter.mzTab` | mzTab export of the filtered consensus result. |

Intermediate idXML files are published under `intermediate/` when `save_intermediate_files = true`. Pipeline reports and process logs are stored in the corresponding timestamped `log-validate-ProteomicPSM-*` directory. Four NFTest integration cases are defined in [nftest.yml](nftest.yml), covering the standard workflow, skipped indexing, disabled Percolator, and unfractionated consensus.

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
