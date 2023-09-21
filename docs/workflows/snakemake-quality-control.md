# snakemake-quality-control
The quality control workflow can be used to generate a report containing the most important quality metrics for any sequencing run. These include the whole genome coverage, contamination estimates, transition transversion ratios and pedigree checks. Afterwards, the output of all quality control tools is collected by the MultiQC tools and rendered into a nice report, which is also collected with the `snakemake --report <OUTPUT>.html` feature.

## Usage
This workflow is designed to be dynamic, and feature a variable set of input files. Nevertheless, as all workflows do, it requires a configuration in `cluster.yaml` and `config.yaml`. The cluster configuration file has to be adapted to fit your computational setup, if you want to run the Parabricks versions of the `MultipleMetrics` and `WGSMetrics` tools for their speedup, make sure the processing nodes for these jobs have appropriate GPUs available.

The main configuration file of the workflow supplies several input and reference file paths for tools like VerifyBamID, which should be adapted to fit your local system. Input for the pipeline is found by searching through the set of directories linked to by `input-dirs`, where each file is picked up if its filename matches any of the glob patterns given in the `search-patterns` option. To afford this flexibility in providing input to the workflow, it makes some assumptions about the structure of input file names, which are mirrored in `rules/meta.smk`. In particular, the sample name should be part of the filename and, unless its a fastq file, is assumed to be the filename without the file type suffix component.  The workflow also assumes that BAM files are named `.bam`, and `vcf` files contain the string `vcf` somewhere in the filename. After finding all possible input files, the workflow will use the aforementioned heuristics to supply the input files to each quality control tool.

The output of the workflow primarily consists of the `MultiQC` output `.html` and table files. We recommend to generate the snakemake report as well, in order to document the names and versions of all tools that were included during workflow execution.
## Tools
The pipeline includes many tools, which are executed on different files types. Currently, the list of tools executed by file type is as follows:

| File Type | Tool                   | Purpose                                                                                                               |
| --------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------- |
| \*.bam    | mosdepth               | Compute mean and median coverage values over the whole genome and individual chromosomes                              |
|           | samtools stats         | A combination of samtools stats, flagstat, and idxstats to collect coverage and other meta information from bam files |
|           | picard WGSMetrics      | Collect WGS Metrics from GATK or NVIDIA Parabricks                                                                    |
|           | picard MultipleMetrics | GATK or NVIDIA MultipleMetrics invocation producing many alignment- and base quality metrics for each bam file        |
|           | VerifyBamID            | Compute contamination estimates based on 1KG variants                                                                 |
| \*.fastq  | FastQC                 | .fastq File quality and duplication metrics                                                                           |
| \*.vcf    | bcftools stats         | Compute transition transversion ratios and base exchange histograms                                                   |
|           | peddy                  | Sex and Pedigree and Ethnicity checking for the cohort                                                                |
|           | vcftools relatedness2  | Compute a relationship matrix for the whole cohort given (rendered nicely by MultiQC)                                 | 
The output of all tools is collected and nicely bound together by `MultiQC`.