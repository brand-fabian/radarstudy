# cDNM Detection
The main hail script to transform any set of cohort VCF files into a single cohort dataframe and perform _de novo_ filtering on the resulting dataset. The script is intended to be run from a Unix commandline that has access to the necessary environment variables to connect to a Spark Cluster that has hail utilities installed.

```
Script:        scripts/find_dnm/find_dnms.py
Spark Version: v3.1.1
Hail Version:  v0.2.89
```
## Prerequisites
To run this script, you will need the hail binaries, and Apache Spark and its associated dependencies. The script has been tested with Hail v0.2.89 and Spark v3.1.1. It is not recommended to run this script with hail's builtin local spark executor, since running times might grow excessively on larger cohorts. This script has been tested and run on an HPC system that utilizes the `EasyBuild` build system with the modules `Miniconda3/4.12.0` and `Spark/3.1.1-foss-2020b` loaded. Please refer to your local installation guides to see how to install this software on your system if you want to run these tools.

The python environment additionally must have the `pandas` and `numpy` libraries installed to be able to run the script successfully.
## Usage
To run this in a valid environment (see [Prerequisites](./cdnm-detection.md#prerequisites)), the following options can be used. The options `-f` (a pedigree file), `-R` (.fasta reference) and any number of named VCF files are required.

```bash
usage: FIND-MSDNS [-h] [-p OUTPUT] -f FAM -R REFERENCE [--max-ac MAX_AC]
                  [--min-p-de-novo MIN_P_DE_NOVO] [--min-aaf MIN_AAF]
                  [--min-transversion-aaf MIN_TRANSVERSION_AAF]
                  [--min-parent-dp MIN_PARENT_DP] [--min-dp MIN_DP]
                  [--max-parent-alleles MAX_PARENT_ALLELES]
                  [--window-size WINDOW_SIZE]
                  [--checkpoint {refined_dnm,data,msdns,r_de_novo,ddn}]
                  [--save-checkpoints] [--no-save-checkpoints] [-t TMP_DIR]
                  [--verbose {debug,info,warning,error}]
                  VCF [VCF ...]

Find msdns in a given vcf file.

positional arguments:
  VCF                   VCF input file (prefer bgziped files) and optionally cohort
                        name separated by colons. Example: radar:/path/to/vcf.bgz

options:
  -h, --help            show this help message and exit
  -p OUTPUT, --output OUTPUT
                        Output file prefix directory (default: output)
  -f FAM, --fam FAM     Pedigree file (plink .fam format) (default: None)
  -R REFERENCE, --reference REFERENCE
                        Reference sequence (default: None)
  --max-ac MAX_AC       Max. allele count in the cohort for a variant to be
                        considered. (default: 10)
  --min-p-de-novo MIN_P_DE_NOVO
                        Min. probability of de novo event (default: 0.8)
  --min-aaf MIN_AAF     Min. alternate allele frequency to call het (default: 0.3)
  --min-transversion-aaf MIN_TRANSVERSION_AAF
                        Min. alternate allele frequencies for transversions (default:
                        0.3)
  --min-parent-dp MIN_PARENT_DP
                        Min. sequencing depth in parents (default: 10)
  --min-dp MIN_DP       Min. sequencing depth in index (default: 15)
  --max-parent-alleles MAX_PARENT_ALLELES
                        Max. reads in parents supporting alt allele (default: 1)
  --window-size WINDOW_SIZE
                        MSDN window size (default: 20)
  --checkpoint {refined_dnm,data,msdns,r_de_novo,ddn}
                        Enable checkpoints by name. (default: None)
  --save-checkpoints
  --no-save-checkpoints
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary directory (default:
                        /home/brand/scratch/spark/5433404/tmp/tmp.EIqNjFnFDh)
  --verbose {debug,info,warning,error}, -v {debug,info,warning,error}
                        Set verbosity (default: info)
```

## Methods
The script regularly outputs checkpoints, which can not be resumed from but denote a specific state in the processing of all input data. In total, there are 6 checkpoints:
1. `data.mt`: The first checkpoint includes all data after the VCF's have been read and merged but no QC or processing has been performed.
2. `data2.mt`: The second checkpoint includes all data from the VCF's after correction of the `PL` values and `split_multi_hts` as well as annotation of QC values (`hl.variant_qc, hl.sample_qc` )
3. `ddn.mt`: The third checkpoint still contains all VCF data, but now data has been annotated with _de novo_ scores from `hl.de_novo`. Up until this checkpoint, the format was not changed from the original `hl.import_vcf` matrix table schema.
4. `r_de_novo.mt`: The first checkpoint containing only _de novo_ mutation data. A static cutoff of `0.7` is used for the _de novo_ probability to generate this table. The format of the table has now been changed to a matrix table with each row denoting one _de novo_ variant with columns containing paternal and maternal variant information at the given site.
5. `refined_dnm.mt`: This checkpoint contains refined _de novo_ information. At this stage, all filters described in the paper and supplemental material were applied. This checkpoint is also written as a `.tsv` file for easier consumption by downstream tools.
6. `msdns.ht`: The last checkpoint containing only cDNM calls. Calls for a cDNM are made if multiple _de novo_ mutations are included within 20 base pairs. All cDNMs are given an ID at this stage. Each row in the table now also contains columns prefixed with `previous.`, which contain the information from the preceding DNM in each cluster.

 The input for this script is any number of named VCF files in the format `<NAME>:<PATH>`, where the name will be added to the hail matrix table as row field `cohort`. In this way, it is possible to merge multiple cohort VCF's that have been joint-called (e.g. using GLnexus) to the same dataset for use in this analysis. At the beginning of the script, each file is read using the `hl.import_vcf` function and merged using `union_cols` with `row_join_type="outer"`. 

To fix some quirks with the GLnexus VCF output format, missing `PL` values in some VCF rows are set to `1000`. In a similar fashion, `split_multi_hts` is used to normalize all rows in the matrix table to exactly one alternate allele.

The _de novo_ filtering part of the script can be influenced by variables passed in from the command line. In general, the default values for each option are used. Other than that, please refer to the main text and supplemental material of the paper to find the options used for the processing of the three study cohorts. Since the first _de novo_ identification step in this script uses a hard cutoff of `p_de_novo >= 0.7`, this is also the lower bound for the `--min-p-de-novo` option of the script.