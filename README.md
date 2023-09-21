# Radarstudy Scripts
This repository contains all analysis scripts used in the processing of the Radar, Chernobyl and Inova cohorts for the Radarstudy. The repository consists of two parts, the workflows in their subdirectory were used for larger processing steps, largely with already established tools. The analysis folder contains scripts that were written for this study to perform the _de novo_ detection and statistical analysis for the study.

Please refer to the [docs](docs/index.md) folder for a detailed description of each workflow and tool.

## Setup

The [workflows](workflows/) directory includes four workflows that are referenced in the main paper and supplemental material of our study. In no particular order:
* [snakemake-dragen](workflows/snakemake-dragen/): The workflow for running Illumina DRAGEN on AWS
* [snakemake-parabricks](workflows/snakemake-parabricks/): The workflow for using NVIDIA Parabricks GPU accelerated mapping and variant calling
* [snakemake-phasing](workflows/snakemake-phasing/): Phasing of _de novo_ mutations using `unfazed` and `WhatsHap`
* [snakemake-quality-control](workflows/snakemake-quality-control/): Quality control for NGS experiments, monitoring coverage, contamination etc.
Alongside these workflows, the [snakemake-wrapper](workflows/snakemake-wrapper/) directory features the snakemake wrappers used to execute each tool in a HPC setting.

The [analysis](analysis/) directory contains the sripts for the detection of (clustered) _de novo_ mutations at [cDNM Detection](analysis/find_dnm/) and the statistical analysis at [cDNM Analysis](analysis/plot).