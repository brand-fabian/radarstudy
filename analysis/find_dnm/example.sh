#! /bin/bash
python3 find_dnm.py \
    -p <OUTPUT PREFIX> \
    --max-ac 1 \
    -f <PEDIGREE FILE> \
    -R <REFERENCE .fasta> \
    RADAR:<PATH TO .vcf.bgz> \
    INOVA:<PAHT TO .vcf.bgz> \
    CRU:<PATH TO .vcf.bgz>