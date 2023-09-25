# Common 
The workflows in this repository are all using `snakemake` to orchestrate the execution and, sometimes, the installation of different tools. Each workflow has its own description, detailing the tools used and the task achieved by the respective workflow. This document describes the common setup that all workflows adhere to and some general aspects around the reproducibility that is expected with this code.

## Usage
All workflows can be executed through the snakemake CLI, please see their documentation for an indepth look at possible options that influence runtime behaviour. The workflows written for this study do not depend on other software, besides snakemake and pysam and should be fairly version agnostic, unless otherwise noted (exception is, of course, DRAGEN on AWS). Therefore, the following conda environment should suffice to execute each workflow:
```yaml
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python
  - snakemake
  - pysam
```
Depending on your environment and assuming access to a HPC cluster and snakemake profile, the snakemake invocation can then look like the following:
```bash
snakemake \
	--use-conda --conda-frontend conda \
	--profile <YOUR_PROFILE> --cluster-config cluster.yaml \
	--keep-going -j 100
```
## Workflow Setup
All workflows share a common setup, with folders for snakemake rules, config schemas, wrappers and a report file each. At the root of each workflow exists a `Snakefile`, which only includes `rule all:`. In all cases, `rule all` is preceded by the inclusion of `rules/meta.smk`, which is a snakemake and python hybrid file taking care of reading all input sample sheets and configurations. After `rule all`, we include all other snakemake rules, to ensure that all rules are picked up in the correct order by snakemake (i.e. `rule all` comes first).

`meta.smk` is a common idiom in all workflows in this repository. It is used to perform setup tasks for the workflow, such as setting some snakemake variables (e.g. `report`, `configfile`, `wildcard_constraints`) and parse and validate any input files such as sample sheets or file tables. In our case, `meta.smk` also includes all input functions for snakemake rules used in the workflow, so it should be the first entrypoint for anyone wishing to understand how particular rules are connected in a given snakemake workflow.
## Workflow Configuration
All workflows are configured by a `config.yaml` and a `cluster.yaml`. The `cluster.yaml` is very site specific and should only be treated as an example of how such a file could look like. The format of `config.yaml` is described in a schema, named `schemas/config.schema.yaml`. The configuration read during execution of the workflow is validated against this schema, as suggested by snakemake best practices. The `config.yaml` file is also usually augmented by one or more samplesheets, typically `samples.tsv` or `files.tsv`, which are also described and validated against their respective schemas in the `schemas/` folder.

In general, to run these workflows it should be sufficient to adapt the `config.yaml`, `cluster.yaml` and any required files or samples tables to the input dataset.
## Wrappers
To aid in transportability and maintainability of these workflows, they make extensive use of snakemake wrappers. Snakemake wrappers are a mechanism by which the implementation of a rule can be disassociated from its use in a workflow. This is useful, because it allows reuse of common components, i.e. the rule body within different workflows. For example, many `bcftools` invocations are part of multiple workflows and its call can be reused through wrappers.

Each wrapper for a tool includes a `meta.yaml`, `environment.yaml` and a `wrapper.py` file. The `meta.yaml` file contains a high level description of the task input and output files and some information about the involved processing steps. `environment.yaml` contains the `conda` environment specification including all necessary software to execute the wrapper successfully from a clean install of any unix node, unless the wrapper is using `EasyBuild` based software. The `wrapper.py` file includes the python code that is necessary in order to parse the snakemake rule inputs, outputs and custom parameters, check the input for potential errors and execute the tool the wrapper is written for. Simple examples of wrappers can be found in the `bcftools/` folder, which includes one wrapper for each supported `bcftools` subcommand (e.g. `view`, `isec` ).

Our own `snakemake-wrappers` repository is present at `<PATH>`. The wrappers are meant to be used with `conda`, or in some cases with `EasyBuild` environments. Wrappers using an `EasyBuild` module usually can be distinguished from standard `conda`-based Wrappers by the fact that their `environment.yaml` is empty and they feature a param called `module_name`.

Some wrappers are written with a specific hardware setup in mind. While they _should_ work well under most circumstances, there might be issues with portability due to missing modules or environment variables. One example for an environment variable that is used by many wrappers  is `$SCRATCH_DIR`, which is used by our HPC provider to give each job a temporary directory on fast, node-local SSDs to write to. Therefore, many wrappers expect there to be a fast and potentially large local temporary storage directory (if `$SCRATCH_DIR` is unset it defaults to the system temp) for execution, which can lead to performance degradations if this assumption is broken.

## Scripts
The `scripts/` folder is used in some workflows to store custom programs written in Python or other programming languages that are used in that workflow. These scripts usually provide some missing functionality or perform more complex tasks than the common wrappers. Nevertheless, each `script` follows the layout of the snakemake wrapper folders, all having at least a `environment.yaml`, `meta.yaml` and `wrapper.py` file. In addition, there is usually a python or other command line tool inside the folders.

## Processing Environments
As discussed in the [Wrappers](../../docs/workflows/index.md#wrappers) section of this document, the workflows are usually written with the local HPC system in mind. In particular `cluster.yaml` and the snakemake profile used for deploying the snakemake jobs to the cluster are highly specific to our local environment. To be portable to other systems, some changes are likely needed in these files. Additionally, the wrappers that make excessive use of the `$SCRATCH_DIR` environment variable possibly need to be overhauled to avoid excessive runtime or memory usage during execution. Workflows without any special hardware requirements (e.g. GPU, FPGA) should work on any commodity hardware, the largest node these workflows had access to was 128 Cores, 512 GB RAM but smaller ones do work for the large majority of snakemake jobs.

The `snakemake-dragen` workflow is an exception of this rule, as it was written with an AWS setup in mind. These files were used to run the DRAGEN pipeline from Illumina on AWS `f1.4xlarge` nodes, as recommended by Illumina. The setup on AWS is similar to the local setup. On AWS, node management was taken care of by the cluster scaling features of the Slurm Workload Manager, which has provisions for starting and deleting nodes for new jobs on demand. The Terraform infrastructure as code snippets that were used for this setup can be made available on request, but have not been included in this repository for brevity and due to security concerns. Due to missing a shared home directory to put files into, this workflow also has to be run with the snakemake options `--default-remote-provider S3` and `--no-shared-fs`.