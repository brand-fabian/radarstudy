#! /bin/bash

sudo yum install -y python36 python36-pip python36-devel gcc
sudo pip3 install --upgrade pip
sudo pip install --upgrade pip
sudo pip3 install snakemake boto3 pandas cookiecutter
sudo pip install snakemake boto3 pandas

##
# Setup an AWS Slurm headnode to run the pipeline
##
mkdir -p ~/.config/snakemake
cd ~/.config/snakemake
cookiecutter https://github.com/Snakemake-Profiles/slurm.git