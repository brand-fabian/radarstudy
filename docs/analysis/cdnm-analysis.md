# cDNM Analysis
Description of the main statistical analysis pipeline for the radarstudy.

For the detection of variants, see: [[Radarstudie/cDNM Detection]]

The analysis is run with four main python scripts:
* [[cDNM Analysis#DNM Analysis]] dnm.py
* [[cDNM Analysis#cDNM Analysis]] msdn.py
* [[cDNM Analysis#Positive Predictive Value]] ppv.py
* [[cDNM Analysis#Phasing]] phasing.py
## Analysis
The following steps are run through the main analysis pipeline scripts. All of those are setup, s.t. the data can be generated on one run, and reused on the next run through the rudimentary `NamedCheckpoint` class. Currently, each analysis _should ideally_ be implemented through the `AnalysisFactory` abstract class and its `AnalysisFactory.get` function. However, many are still present as simple function in all the different scripts. The statistics are computed, retrieved and saved by using the `StatisticsFactory.register` decorator. Each function should be annotated with a `@NamedCheckpoint.checkpoint()` and, if it outputs something where a statistical comparison should be made (i.e. t-test, Mann-Whitney-U-Test), a second annotation `@StatisticsFactory.register` should be made.

>[!info] Nomenclature
>During the writing of the manuscript, the term for clustered _de novo_ mutations (cDNMs) was changed from the earlier multisite _de novo_ mutation (MSDN). Since scripts were written prior to this change, most of them still use the old MSDN moniker for clusters.

## Common Components

### Prerequisites
The analysis scripts is `hail` and `pandas` to represent data, `plotly` for creating figures and `statsmodels` for the statistical analysis. All of these are python libraries, an example of a compatible environment is given in `environment.yaml`. Note that, while possible, we still recommend not to run hail with the builtin local spark executor, but instead use a custom Spark Cluster with hail installed for greater stability and performance.
All analyses have been performed on an HPC cluster, which had built software using the `EasyBuild` framework. To run the scripts, the `EasyBuild` modules `Miniconda3/4.12.0` and `Spark/3.1.1-foss-2020b` have been loaded, and additionally we used a minimal conda environment with these python dependencies:
```yaml
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.10.2
  - pip
  - pip:
    - numpy==1.22.4
    - pandas==1.3.5
    - patsy==0.5.2
    - scipy==1.7.3
    - statsmodels==0.13.2
```

### Config
All scripts should make use of the constants in `meta.config` and `meta.plotting` to automatically configure most parts of the plotly library. The most influential variables for the resulting plots and statistics are `LEVELS_BASE` and `CATEGORY_ORDER_BASE`, which are used to determine the order in which the cohorts appear in the plots and which cohort is used for the Intercept in generalized linear models (i.e. gets the value `0` in categorical data).

### Checkpoints
The scripts in this folder feature a powerful checkpointing mechanism, that is used to write out almost all output files. In particular, refined data tables (from hail and pandas) are cached, and can be reused in subsequent invocations of these scripts. This reduces runtime in most common scenarios if only one or two plots were changed. However, data is only read for these tabular checkpoints (e.g. hail or pandas tables), the raw data for each plot is newly computed  each time the script is executed, to guarantee that the plot reflects the data and changes to the code.
Checkpoints are also used to save plotly plots and statsmodels tables to the filesystem. These are the primary results of these scripts and generated during the runtime of the script. If a `plotly.Figure` instance is passed to a checkpoint, the checkpoint will store the raw data as `.csv` and `.tex` file and the figure itself in `.png`, `.svg` and `.html` format. When a `statsmodels.Results` instance is checkpointed, it is saved as `.txt`, `.csv` and `.tex` file and contains the estimated model parameters and some descriptive parameters in a tabular format.

### Statistics Factory
The statistics factory is used to perform simple tests, like the t-test or Mann-Whitney U-tests on each dataset. Any function that is annotated with a call to `StatisticsFactory.register` will be added to a table which is printed at the end of each run of any of the scripts in this folder. Most commonly, the statistical tests are performed on the exact data tables used for plotting, by stratifying all samples by the `cohort` column. All these options are configurable in the call to the register function, which passes any unused arguments on to the function used to execute each statistical test.

>[!note]
>The Statistics Factory is only included in this code repository for legacy reasons. All statistical tests reported on in the paper and supplemental material feature negative binomial regression models, which we deemed a better fit for our data after assessing different options.

### Data Loading
All statistical analysis scripts feature common code for loading the different DNM and cDNM callsets. Additionally, the scripts for analyzing phasing and positive predictive value data, may include some other files that have been generated by other tools.
Each loader performs similar steps, namely reading in the `hail` (matrix) table or `pandas` DataFrame, and performing some additional quality control checks. In the case of DNM and cDNM calls, these checks include discarding some samples from the analysis that have failed the quality control and applying the allele frequency filters to correct for sequencing specific artifacts. The loader scripts also apply the Age Matching and downsampling of all three cohorts to match the Radar cohort in size and distribution of paternal and maternal age at conception.
### Arguments
The analysis script have very similar command line options. Aside from required arguments that are needed to perform the processing task of each given script, most option are concerned with metadata loading and how matching is performed. Common Options include those related to loading the metadata, which can be done either by specifying the original tables (`--[radar,inova,cru,pilot]-meta`) or by specifying the processed metadata from an earlier run (`--metadata`). Additionally, some attributes of the DNM and cDNM loaders can be changed by setting the `--isolated` (only count DNMs if they are not part of a cluster), `--[no]-apply-graphtyper-filter` (filter DNM set to a set of validated ones) or `--[no]-control-matching` (behaviour of the age matching, if enabled case and controls are paired for all tests) flags. The common `--language` option allows changing the output language of all plots produced by the scripts.

## DNM Analysis
The DNM analysis script will filter DNM's to remove potentiala sequencing artifacts and exclude some samples based on QC parameters (see [[cDNM Analysis#Data Loading]]). Afterwards, this script computes the (isolated) DNM rates for all samples in all cohorts and compares them using GLM's and the `StatisticsFactory` methods. The program also computes the paternal age effect, by fitting the model `iDNMs ~ cohort + father_age + (cohort * father_age)` to the input data.

Note that isolated DNMs are only considered if the flag `--isolated` is passed, otherwise all variants that are part of clusters are also included here and therefore counted twice in the statistics.
```bash
Script:        scripts/plot/dnm.py
Spark Version: v3.1.1
Hail Version:  v0.2.89
```
### Usage
```bash
usage: PLOT-DNM [-h] [--msdn-ht MSDN_HT] [--isolated] -R REFERENCE
                [--radar-meta RADAR_META] [--inova-meta INOVA_META]
                [--cru-meta CRU_META] [--pilot-meta PILOT_META] [--metadata METADATA]
                [-f FACTOR] [--language {de,en}] [--graphtyper GRAPHTYPER]
                [--apply-graphtyper-filter] [--no-apply-graphtyper-filter]
                [--apply-control-matching] [--no-apply-control-matching]
                [--load-cache] [--no-load-cache] [-t TMP_DIR]
                [--verbose {debug,info,warning,error}]
                MT

Plot descriptive statistics for dnms.

positional arguments:
  MT                    Input matrix table (refined_dnm) containing de novo mutations
                        for both cohorts.

options:
  -h, --help            show this help message and exit
  --msdn-ht MSDN_HT     Input hail table (msdn_ht) containing multisite de novo
                        mutations for both cohorts. Required for --isolated (default:
                        None)
  --isolated            If true, compute isolated dnm rates (i.e. exlude msdn
                        variants). (default: False)
  -R REFERENCE, --reference REFERENCE
                        Reference sequence (default: None)
  --radar-meta RADAR_META
                        Path to radarstudy metadata table (default: None)
  --inova-meta INOVA_META
                        Path to inova metadata table (default: None)
  --cru-meta CRU_META   Path to trio-cru metadata table (default: None)
  --pilot-meta PILOT_META
                        Path to the pilot study metadata table (default: None)
  --metadata METADATA   Path to a pandas dataframe .pickle file containing study
                        metadata (default: None)
  -f FACTOR, --factor FACTOR
                        Control cohort downsampling factor. If set to 0 (the
                        default), no matching is performed. (default: 0)
  --language {de,en}    Language of axis descriptions in generated plots (default:
                        de)
  --graphtyper GRAPHTYPER
                        Path to pickle files with graphtyper information for the
                        validation of DNMs. (default: None)
  --apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --no-apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --no-apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --load-cache
  --no-load-cache
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary directory (default:
                        /home/brand/scratch/spark/5433404/tmp/tmp.EIqNjFnFDh)
  --verbose {debug,info,warning,error}, -v {debug,info,warning,error}
                        Set verbosity (default: info)
```
## cDNM Analysis
The cDNM analysis script is used for all comparisons of cDNM rates between the cohorts and to model the effects of ionizing radiation on the average number of clusters. Like the DNM Analysis script, some preprocessing of the input data is performed, to exclude QC fails and likely sequencing artifacts from the data. Then, this filtered data is used in negative binomial regression models and the `StatisticsFactory` to assess the influence of paternal exposure to ionizing radiation. 
```
Script:        scripts/plot/msdn.py
Spark Version: v3.1.1
Hail Version:  v0.2.89
```
### Usage
```bash
usage: PLOT-MSDN [-h] -R REFERENCE [--radar-meta RADAR_META]
                 [--inova-meta INOVA_META] [--cru-meta CRU_META]
                 [--pilot-meta PILOT_META] [--metadata METADATA] [-f FACTOR]
                 [--language {de,en}] [--graphtyper GRAPHTYPER]
                 [--apply-graphtyper-filter] [--no-apply-graphtyper-filter]
                 [--apply-control-matching] [--no-apply-control-matching]
                 [--load-cache] [--no-load-cache] [-t TMP_DIR]
                 [--verbose {debug,info,warning,error}]
                 HT

Plot descriptive statistics for msdns.

positional arguments:
  HT                    Input hail table (msdn_ht) containing multisite de novo
                        mutations for both cohorts.

options:
  -h, --help            show this help message and exit
  -R REFERENCE, --reference REFERENCE
                        Reference sequence (default: None)
  --radar-meta RADAR_META
                        Path to radarstudy metadata table (default: None)
  --inova-meta INOVA_META
                        Path to inova metadata table (default: None)
  --cru-meta CRU_META   Path to trio-cru metadata table (default: None)
  --pilot-meta PILOT_META
                        Path to the pilot study metadata table (default: None)
  --metadata METADATA   Path to a pandas dataframe .pickle file containing study
                        metadata (default: None)
  -f FACTOR, --factor FACTOR
                        Control cohort downsampling factor. If set to 0 (the
                        default), no matching is performed. (default: 0)
  --language {de,en}    Language of axis descriptions in generated plots (default:
                        de)
  --graphtyper GRAPHTYPER
                        Path to pickle files with graphtyper information for the
                        validation of DNMs. (default: None)
  --apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --no-apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --no-apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --load-cache
  --no-load-cache
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary directory (default:
                        /home/brand/scratch/spark/5433404/tmp/tmp.EIqNjFnFDh)
  --verbose {debug,info,warning,error}, -v {debug,info,warning,error}
                        Set verbosity (default: info)
```
## Positive Predictive Value
To identify the potential consequences of the low positive predictive value found by our validation efforts, we decided to add a subsampling-based simulation script to our work. This script takes in the positive predictive value (`0.25`  by default) and performs some of the statistical analysis from both, the DNM analysis and cDNM analysis scripts, again. Each simulation run performs the statistical tests and computations of the cDNM Analysis `cDNMs ~ cohort` on subsampled data and the script by default computes 1000 iterations of these simulations.
```bash
Script:        scripts/plot/msdn.py
Spark Version: v3.1.1
Hail Version:  v0.2.89
```
### Usage
```bash
usage: PLOT-PPV [-h] [-n NUM_SIMULATIONS] [-p PPV] -R REFERENCE
                [--radar-meta RADAR_META] [--inova-meta INOVA_META]
                [--cru-meta CRU_META] [--pilot-meta PILOT_META] [--metadata METADATA]
                [-f FACTOR] [--language {de,en}] [--graphtyper GRAPHTYPER]
                [--apply-graphtyper-filter] [--no-apply-graphtyper-filter]
                [--apply-control-matching] [--no-apply-control-matching]
                [--load-cache] [--no-load-cache] [-t TMP_DIR]
                [--verbose {debug,info,warning,error}]
                HT

Simulate the effect of reduced positive predictive values on the statistics.

positional arguments:
  HT                    Input hail table (msdn_ht) containing multisite de novo
                        mutations for both cohorts.

options:
  -h, --help            show this help message and exit
  -n NUM_SIMULATIONS, --num-simulations NUM_SIMULATIONS
                        Number of simulation runs. (default: 1000)
  -p PPV, --ppv PPV     Positive predictive value to simulate (default: 0.25)
  -R REFERENCE, --reference REFERENCE
                        Reference sequence (default: None)
  --radar-meta RADAR_META
                        Path to radarstudy metadata table (default: None)
  --inova-meta INOVA_META
                        Path to inova metadata table (default: None)
  --cru-meta CRU_META   Path to trio-cru metadata table (default: None)
  --pilot-meta PILOT_META
                        Path to the pilot study metadata table (default: None)
  --metadata METADATA   Path to a pandas dataframe .pickle file containing study
                        metadata (default: None)
  -f FACTOR, --factor FACTOR
                        Control cohort downsampling factor. If set to 0 (the
                        default), no matching is performed. (default: 0)
  --language {de,en}    Language of axis descriptions in generated plots (default:
                        de)
  --graphtyper GRAPHTYPER
                        Path to pickle files with graphtyper information for the
                        validation of DNMs. (default: None)
  --apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --no-apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --no-apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --load-cache
  --no-load-cache
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary directory (default:
                        /home/brand/scratch/spark/5433404/tmp/tmp.EIqNjFnFDh)
  --verbose {debug,info,warning,error}, -v {debug,info,warning,error}
                        Set verbosity (default: info)
```
## Phasing
Using `unfazed` and `WhatsHap`, we tried to determine the parental origin of all cDNM clusters identified by the cDNM detection scripts. Our pipeline produces a table as a result, which is analyzed using this tool. Each DNM is annotated with any of three values (`maternal`, `paternal` or `unknown`), while cDNMs can also have a state of `contradiction`, if we had encountered any clusters where one lesions was inherited from the father and one from the father. As detailed in the paper and supplement, this was not the case and this script therefore performs just the basic checks. Due to the low number of phased clusters in the control cohort, the statistics of this script have very little power and should be read with great caution.
```bash
Script:        scripts/plot/msdn.py
Spark Version: v3.1.1
Hail Version:  v0.2.89
```
### Usage
The usage of this tool depends on the output of the phasing workflow, which is also part of this repository. The output file from the workflow must be passed in using the option `-p, --phasing`.
Link: 
```bash
usage: PHASING [-h] -p PHASING -R REFERENCE [--radar-meta RADAR_META]
               [--inova-meta INOVA_META] [--cru-meta CRU_META]
               [--pilot-meta PILOT_META] [--metadata METADATA] [-f FACTOR]
               [--language {de,en}] [--graphtyper GRAPHTYPER]
               [--apply-graphtyper-filter] [--no-apply-graphtyper-filter]
               [--apply-control-matching] [--no-apply-control-matching]
               [--load-cache] [--no-load-cache] [-t TMP_DIR]
               [--verbose {debug,info,warning,error}]
               DNM_MT MSDN_HT

Plot phasing analysis.

positional arguments:
  DNM_MT                Input matrix table (refined_dnm) containing de novo mutations
                        for both cohorts.
  MSDN_HT               Input hail table (msdn_ht) containing multisite de novo
                        mutations for both cohorts.

options:
  -h, --help            show this help message and exit
  -p PHASING, --phasing PHASING
                        .pickle file containing phase information for dnms. (default:
                        None)
  -R REFERENCE, --reference REFERENCE
                        Reference sequence (default: None)
  --radar-meta RADAR_META
                        Path to radarstudy metadata table (default: None)
  --inova-meta INOVA_META
                        Path to inova metadata table (default: None)
  --cru-meta CRU_META   Path to trio-cru metadata table (default: None)
  --pilot-meta PILOT_META
                        Path to the pilot study metadata table (default: None)
  --metadata METADATA   Path to a pandas dataframe .pickle file containing study
                        metadata (default: None)
  -f FACTOR, --factor FACTOR
                        Control cohort downsampling factor. If set to 0 (the
                        default), no matching is performed. (default: 0)
  --language {de,en}    Language of axis descriptions in generated plots (default:
                        de)
  --graphtyper GRAPHTYPER
                        Path to pickle files with graphtyper information for the
                        validation of DNMs. (default: None)
  --apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --no-apply-graphtyper-filter
                        If this flag is set, variants not confirmed by graphtyper
                        will be ignored. (default: False)
  --apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --no-apply-control-matching
                        If this flag is set the control matching will split the
                        control cohort for some plots (default: False)
  --load-cache
  --no-load-cache
  -t TMP_DIR, --tmp-dir TMP_DIR
                        Temporary directory (default:
                        /home/brand/scratch/spark/5433404/tmp/tmp.EIqNjFnFDh)
  --verbose {debug,info,warning,error}, -v {debug,info,warning,error}
                        Set verbosity (default: info)
```

