import pandas
from snakemake.utils import validate
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.logging import logger
import os
import boto3
import sys

##
# Snakemake Setup
#
report: '../report/dragen.rst'

##
# Configuration files
#
configfile: 'config.yaml'
validate(config, schema='../schemas/config.schema.yml')

samples = pandas.read_csv(
    config['samples'], sep='\t'
).set_index('sample_id', drop=False)
validate(samples, schema='../schemas/samples.schema.yml')

if 'files' in config and config['files']:
    input_files = pandas.read_csv(
        config['files'], sep='\t'
    ).set_index('sample_id', drop=False)
    validate(input_files, schema='../schemas/files.schema.yml')

    ## 
    # Subset samples to those with input files
    #
    input_file_sample_ids = set(input_files.sample_id)
    samples = samples[samples.sample_id.apply(lambda x: x in input_file_sample_ids)]
else:
    input_files = None

##
# Remote files setup
#
S3 = S3RemoteProvider(stay_on_remote=True, keep_local=True)
s3_boto = boto3.resource('s3')

logger.info(f"Input Bucket: {config['s3']['input'].lstrip('s3://')}")
logger.info('Output bucket: {}'.format(config['s3']['output'].lstrip('s3://')))

input_bucket = s3_boto.Bucket(config['s3']['input'].lstrip('s3://'))
output_bucket = s3_boto.Bucket(config['s3']['output'].lstrip('s3://'))
remote_files = {
    'input': list(input_bucket.objects.all()),
    'output': list(output_bucket.objects.all()),
}

logger.info('Found {} files in remote buckets'.format(sum(map(len, remote_files.values()))))

##
# Check input files present
#
missing_files = set()
file_keys = set(map(lambda x: x.key, remote_files['input']))
if 'fastq1' in samples and 'fastq2' in samples:
    required_files = set(samples.fastq1) | set(samples.fastq2)
elif 'fastq1' in samples and not 'fastq2' in samples:
    required_files = set(samples.fastq1)
else:
    if 'files' in config:
        required_files = set(input_files.path)
    else:
        raise Exception('No input files provided.')
for s in required_files:
    if not s.split('/')[3] in file_keys:
        missing_files.add(s)

if len(missing_files) > 0:
    logger.info(f"Missing {len(missing_files)} files.")
    logger.error(f"Missing files:\n\n{missing_files}")
    raise Exception("Missing input files.")


##
# DRAGEN license
#
license_text = None
with open(config['license_file'], 'r') as f:
    license_text = f.read().replace('\n', '').strip()

if license_text is None:
    raise Exception('No license key found.')

##
# Wildcard constraints
#
families = samples.groupby('family_id').groups.keys()

wildcard_constraints:
    family="|".join(families),
    sample="|".join(samples['sample_id'].dropna())

##
# Helper functions
#
def is_valid(value):
    return value and pandas.notna(value)

def transform_sex(value):
    if value == 'M':
        return 1
    elif value == 'F':
        return 2
    else:
        return 0

def find_remote_file(sample_id, bucket_files=remote_files['input'], ending=None):
    for r_file in bucket_files:
        if ending is None or r_file.key.endswith(ending):
            if r_file.key.startswith(sample_id):
                return r_file
    return None

##
# Input functions
#
def get_vc_input(wildcards):
    sample = samples.loc[wildcards.sample]
    if 'fastq1' in sample and sample.fastq1.strip().endswith('bam'):
        return {
            'bam': S3.remote(sample['fastq1'], stay_on_remote=False),
        }
    elif input_files is not None:
        return {
            'fq': list(map(lambda s: S3.remote(s.lstrip('s3://')), filter(
                lambda x: x.endswith('fq.gz') or x.endswith('fastq.gz'),
                input_files.loc[wildcards.sample].path
            )))
        }
    else:
        if ('fastq1' in sample and is_valid(sample['fastq1'])
                and 'fastq2' in sample and is_valid(sample['fastq2'])):
            return {
                'fq': [ S3.remote(sample['fastq1']), S3.remote(sample['fastq2']) ],
            }
        else:
            raise Exception("No fastqs provided. Please provide a files table or two fastq columns in the sample table.")

def gather_samples(wildcards):
    family = samples[samples['family_id'] == wildcards.family]
    variants = []
    tbi = []
    for idx, row in family.iterrows():
        variants.append(
            "output/{sample}/{sample}.hard-filtered.gvcf.gz".format(
                sample=row['sample_id']
            )
        )
        tbi.append(
            "output/{sample}/{sample}.hard-filtered.gvcf.gz.tbi".format(
                sample=row['sample_id']
            )
        )
    return {
        'variants': variants,
        'tbi': tbi,
        'ped': "output/{family}/{family}.ped".format(family=wildcards.family),
    }

def cnv_input(wildcards):
    """Find and return panel of normals and sample matched for sex."""
    if config['self_normalize']:
        return {
            'bam': 'output/{sample}/{sample}.bam'.format(sample=wildcards.sample),
            'bai': 'output/{sample}/{sample}.bam.bai'.format(sample=wildcards.sample),
        }
    else:
        raise Exception('Not supported')
        sex_matched = samples[samples['sex'] == samples.loc[wildcards.sample].sex]
        if (len([ x for x in remote_files if x.key == '{s}/{s}.target.counts.gc-corrected'.format(s=wildcards.sample)]) > 0
                and not os.path.exists('output/{s}/{s}.target.counts.gc-corrected'.format(s=wildcards.sample))):
            case = S3.remote(config['s3']['output'] + '/{s}/{s}.target.counts.gc-corrected'.format(s=wildcards.sample), stay_on_remote=False)
        else:
            case = 'output/{s}/{s}.target.counts.gc-corrected'.format(s=wildcards.sample)
        if len(sex_matched) < 2:
            return {
                'case': case,
                'normals': [],
            }
        normal_samples = [ 
            s.sample_id for _, s in sex_matched.iterrows() if not s.sample_id == wildcards.sample
        ]
        return {
            'case': case,
            'normals': list(map(
                lambda x: (
                    S3.remote(
                        config['s3']['output']
                                + '/{s}/{s}.target.counts.gc-corrected'.format(s=x),
                            stay_on_remote=False) 
                    if (len([ x for x in remote_files if x.key == '{s}/{s}.target.counts.gc-corrected'.format(s=x) ]) > 0
                        and not os.path.exists('output/{s}/{s}.target.counts.gc-corrected'.format(s=x))) 
                    else 'output/{s}/{s}.target.counts.gc-corrected'.format(s=x)
                ),
                normal_samples
            ))
        }


def get_dnm_cnv_input(wildcards):
    """For a given family, return cnv input counts files."""
    family = samples[samples['family_id'] == wildcards.family]
    counts = []

    for s, row in family.iterrows():
        counts.append(
            'output/{s}/{s}.tn.tsv'.format(s=row['sample_id'])
        )
    return {
        'counts': counts,
        'pedigree': "output/{family}/{family}.ped".format(family=wildcards.family),
    }

def get_sv_input(wildcards):
    """Given a family and case get bams and pedigree."""
    family = samples[samples['family_id'] == wildcards.family]
    bams = []
    bais = []
    for idx, row in family.iterrows():
        bams.append(
            'output/{sample}/{sample}.bam'.format(sample=row['sample_id'])
        )
        bais.append(
            'output/{sample}/{sample}.bam.bai'.format(sample=row['sample_id'])
        )
    return {
        'bams': bams,
        'bais': bais, 
        'pedigree': "output/{family}/{family}.ped".format(family=wildcards.family),
    }

