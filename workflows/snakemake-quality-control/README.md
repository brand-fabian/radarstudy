# Quality Control

Snakemake pipeline to execute various quality control tools on a set of files
found in a set of input folders. Based on scripts used for the original 
processing of WGS-Data for the radarstudy.

## Tools

Currently, the following quality control tools and metrics are collected:

* fastqc
* bcftools stats
* ngs-chew
* peddy
* samtools stats
* verifybamid
* wgs-metrics
* mosdepth