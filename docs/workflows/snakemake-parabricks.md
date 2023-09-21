# snakemake-parabricks
The snakemake-parabricks is the second workflow used to process `.fastq` or `.bam` files to VCF. It uses the NVIDIA Parabricks toolkit to perform GPU accelerated mapping and variant calling in a manner closely resembling GATK best practices pipelines.

Tested Parabricks Versions: `3.6.1`, `3.7.1` (discontinued by NVIDIA, License Required), `4.0.0`, `4.1.0`

## Hardware Requirements
The execution of most parabricks tools for alignment, variant calling and quality control requires access to a NVIDIA enterprise-grade GPU with enough video memory. We tested and executed this workflow NVIDIA A100 and V100S video cards, with the former providing a significant speedup for the workflows. From our experience, roughly 24GB of GPU memory are sufficient to run a parabricks alignment job, and excess memory provides little benefit in terms of computation speed.

## Usage
This workflow follows the conventions layed out in [[Common]]. To adapt this pipeline for your use, or to reproduce our processing results, you have to adapt the `cluster.yml` and `config.yml` files to your environment. `cluster.yml` needs to be adjusted, s.t. the correct resources are requested from your scheduler to gain access to GPU's on each processing node. A schema for `config.yml` is given in the `schemas/` folder, but generally all paths to reference files like genome builds and known variants should be adjusted to fit your local setup.

All input files and sample relationships are passed in using the `files.tsv` and `samples.tsv` tables. Each of those also has schemas in the `schemas/` folder. In short, `samples.tsv` resembles a standard pedigree file, used by the workflow primarily to retrieve sample ids. `files.tsv` details input paths for all `.fastq` or `.bam` input files. Examples for both files are present inside the root folder of the workflow.

If there are `.bam` input files for the workflow, we first use the Parabricks toolkit to convert these files back to their original `.fastq` Format, since NVIDIA Parabricks currently (v4.1.1) can not re-align BAM files directly.