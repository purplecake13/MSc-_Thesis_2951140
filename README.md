# MSc-_Thesis_2951140
Exploring Gut Microbiome Signatures in Parkinson’s and Alzheimer’s Diseases for Potential Machine Learning-Based Classification

This repository contains the full processing and analysis pipeline used to investigate gut microbiome signatures in Parkinson's disease (PD) and Alzheimer's disease (AD) via cross-platform meta-analysis and machine learning classification. It is provided for examiner review and reproducibility reference. **Raw sequencing data is not included**.

## Study design

- **PD arm**: shotgun metagenomics, 5 studies —
  Bedarf 2017, Boktor 2023 (TBC and RUMC cohorts, treated as separate `study_id`s due to different sequencers/encodings), Clasen 2024, Mao 2021, Wallen 2022.
- **AD arm**: 16S rRNA amplicon sequencing, 7 studies retained for ML after exclusions (see below).
- **Common ground**: genus-level taxonomic collapse, used for all cross-platform comparison.
- **Mixed-platform design is a deliberate, documented limitation**, not an oversight. No usable AD shotgun dataset existed at sufficient sample size and accessibility (see Methods in the thesis for the full justification). Because of this, PD and AD are **never compared directly** — only platform-matched disease-vs-healthy-control comparisons are reported (PD-vs-HC, AD-vs-HC), each analysed and modelled separately.

## Repository structure

```
00_shared_preprocessing/     Download, merge, QC (FastQC/MultiQC/fastp), antibiotic filtering, sample counting
01_ad_16s_pipeline/          ad_* — QIIME2 import, DADA2, region-matched SILVA 138.2 taxonomy assignment, genus collapse, PICRUSt2
02_pd_shotgun_pipeline/      pd_* — host removal, MetaPhlAn4, HUMAnN3, genus-level table construction
03_cross_platform_analysis/  cx_* — metadata joining, platform merge, batch correction, diversity, ANCOM-BC2, MaAsLin2, Cytoscape network export

ml/                          Machine learning (run locally, not on HPC)
  final/                       Two-model architecture that produced the thesis's reported numbers
  earlier_iterations/          Superseded three-class-classifier approach, kept for provenance

metadata/                      master_metadata.csv, all_metadata.xlsx, per-study sample lists

environment/                   environment.yml (conda env export)
```

### Naming convention

- Numeric prefixes (`01_`, `02a_`, ...) — shared preprocessing steps, run
  before the AD/PD arms diverge.
- `ad_NN_` — AD 16S pipeline, in execution order.
- `pd_NN_` — PD shotgun pipeline, in execution order.
- `cx_NN_` — cross-platform steps, after both arms reach genus level.
- Paired `.sh` / `.py` (or `.sh` / `.R`) files: the `.sh` is the SLURM submission wrapper for the corresponding `.py`/`.R` script.
- Files suffixed `_verify` or `verify_*` are post-hoc sanity checks (job completion, read counts, file integrity) run after each major step, not part of the main analytical path.

## Pipeline overview (execution order)

1. **Shared preprocessing** — download, merge split files, FastQC/MultiQC, fastp/DADA2-ready trimming, antibiotic-use filtering.
2. **AD arm** — QIIME2 import → DADA2 denoising → taxonomy assignment with region-matched SILVA 138.2 classifiers (custom-trained V1-V2 classifier for Yamashiro 2024) → genus collapse → PICRUSt2 functional prediction.
3. **PD arm** — host read removal (Bowtie2) → MetaPhlAn4 taxonomic profiling → genus-level table construction → HUMAnN3 functional profiling (run directly on all datasets, including the full 724-sample Wallen 2022 cohort, per supervisor instruction).
4. **Cross-platform** — merge PD + AD genus tables → batch correction → diversity analysis → differential abundance → pathway analysis → network export.
5. **Machine learning** (local, `notebooks/ml/final/`) — LASSO feature selection → two separate Random Forest models (PD-vs-HC, AD-vs-HC) with `StratifiedGroupKFold`/Leave-One-Study-Out cross-validation → SHAP (TreeExplainer) interpretation.

## Data availability

Raw sequencing data is not redistributed in this repository. All datasets are public; accession numbers and sources are listed in `metadata/master_metadata.csv` and in the thesis Methods section.
Intermediate large files (`.qza`/`.qzv`/`.biom`, MetaPhlAn/HUMAnN3 outputs) are also excluded — regenerate them by running the pipeline in order against the original accessions.

## Environment

- HPC: BlueBEAR (University of Birmingham), SLURM scheduler.
- QIIME2 2025.4, DADA2, PICRUSt2, SILVA 138.2 classifiers.
- MetaPhlAn4 (`mpa_vOct22_CHOCOPhlAnSGB_202403`), HUMAnN3 (`humann/4.0.0a1-foss-2023a`, EC-filtered UniRef90).
- R: ConQuR, MMUPHin, ANCOM-BC2, MaAsLin2, ComplexHeatmap.
- Python (local): scikit-learn, SHAP, pandas — see `environment/environment.yml`.

## Limitations

See the thesis Discussion/Limitations section for the full treatment. In brief: mixed-platform design (no shotgun-shotgun AD comparison possible), cross-sectional data, potential medication confounding, and the inherent assumptions of batch-correction methods across sequencing platforms.

## Contact

Questions about this code can be directed to aishwarya13p@gmail.com or raised as a GitHub issue on this repository.
