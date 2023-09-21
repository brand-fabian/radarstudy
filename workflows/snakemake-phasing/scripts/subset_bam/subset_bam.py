import pandas
import os
import argparse
import logging
import typing
import pysam
import sys
from enum import Enum
from collections import namedtuple
from itertools import chain
from time import time
from map import fix_readgroup, get_reads, merge_readgroup, GetReadsOptions
from meta import ReferenceMatcher, Forklift

###############################################################################
# Logging                                                                     #
###############################################################################
logger = logging.getLogger("subset-bam")
logger.setLevel(logging.ERROR)
channel = logging.StreamHandler()
formatter = logging.Formatter("[%(asctime)s - %(name)s - %(levelname)7s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
channel.setFormatter(formatter)
logger.addHandler(channel)

log_level = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
}

###############################################################################
# Arguments                                                                   #
###############################################################################
Bam = namedtuple("Bam", [ "sample", "path" ])
def input_bam(val: str) -> Bam:
    split_val = val.split(":", 1)
    if len(split_val) != 2:
        raise argparse.ArgumentTypeError(
            "invalid bam option \"{}\". Please use the format {{sample_name}}:{{path}}.".format(val))
    sample, path = split_val
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise Exception("file {} not found".format(path))
    return Bam(sample=sample, path=path)

parser = argparse.ArgumentParser(
    prog="SUBSET-BAM",
    description='Subset bam files and prepare them as input for WhatsHap or Unfazed pedigree based calling.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bam', metavar='BAM', type=input_bam, nargs='+',
                    help='List of input bam files for subsetting. Format: {{sample_name}}:{{path}}')
parser.add_argument('-S', '--sites', type=str, required=True,
                    help='Sites to include in the final bam file')
parser.add_argument('--fam', type=str, default=None,
                    help='Pedigree input file.')
parser.add_argument('-F', '--family-id', type=str, default=None,
                    help='Family ID from the pedigree to retrieve variants from.')
parser.add_argument('-R', '--reference', type=str, default=None, choices=[ "grch37", "grch38" ],
                    help='Reference that samples should be aligned to')
parser.add_argument('--grch37', type=str, default=None,
                    help='GRCh37 reference genome .fasta file (+ bwa indexes)')
parser.add_argument('--grch38', type=str, default=None,
                    help='GRCh38 reference genome .fasta file (+ bwa indexes)')
parser.add_argument('--sample-id', type=str, default=None,
                    help='Sample ID to retrieve variants from.')
parser.add_argument('--flanking-region', type=int, default=5000,
                    help='Size of the flanking region for each site to be included in the final bam.')
parser.add_argument('--fetch-unmapped-pairs', default=False, action='store_true',
                    help='If true, retrieve all mate pairs from the bam file if alignment is performed (very slow)')
parser.add_argument('--no-ignore-unmatched-reads', default=True, action='store_false',
                    dest="ignore_unmatched_reads",
                    help='Emit the unmatched singleton reads for all regions instead of discarding them')
parser.add_argument('-T', '--threads', type=int, default=1,
                    help="Number of threads to use for mapping")
parser.add_argument('-t', '--tmp-dir', help='Temporary directory',
                    type=str, default=os.getenv('TMPDIR', '/tmp'))
parser.add_argument('--verbose', '-v', help='Set verbosity',
                    choices=log_level.keys(), default='info')

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

###############################################################################
# Main script                                                                 #
###############################################################################
if not os.path.isfile(args.sites):
    logger.error("missing sites input file {}".format(args.sites))
    sys.exit(1)

if args.family_id is not None and not os.path.isfile(args.fam):
    logger.error("missing input pedigree file {}".format(args.fam))
    sys.exit(1)

if args.sample_id is None and args.family_id is None:
    logger.error("missing family_id and sample_id")
    sys.exit(1)

if args.reference == "grch37" and not os.path.isfile(args.grch37):
    logger.error("missing grch37 reference file {}".format(args.grch37))
    sys.exit(1)

if args.reference == "grch38" and not os.path.isfile(args.grch38):
    logger.error("missing grch38 reference file {}".format(args.grch38))
    sys.exit(1)

dnm_df = pandas.read_csv(args.sites, sep="\t")
ped = pandas.read_csv(
    args.fam,
    sep="\t",
    names=["family_id", "s", "father", "mother", "sex", "is_affected"]
).set_index("s", drop=False)

family_id = args.family_id
sample_id = args.sample_id
if args.family_id is None and args.sample_id is not None:
    family_id = ped.loc[args.sample_id].family_id

if args.family_id is not None and args.family_id not in set(ped.family_id):
    logger.error("no pedigree information for {}".format(args.family_id))
    sys.exit(1)

if args.sample_id is not None and args.sample_id not in set(dnm_df.s) and family_id is None:
    logger.error("no sample information for {}".format(args.sample_id))
    sys.exit(1)

if family_id is not None:
    trios = ped[(ped.family_id == family_id)]
    dnm_df = dnm_df[dnm_df.s.isin(set(trios.s))]
elif args.sample_id is not None:
    dnm_df = dnm_df[dnm_df.s == args.sample_id]


input_files = { bam.sample: bam.path for bam in args.bam }

get_reads_opts = GetReadsOptions(
    flanking_region=args.flanking_region,
    fetch_unmapped_reads=args.fetch_unmapped_pairs,
    ignore_unmatched_reads=args.ignore_unmatched_reads,
    threads=args.threads,
)

matcher = ReferenceMatcher()
if args.grch37 is not None:
    matcher.add("grch37", args.grch37)
if args.grch38 is not None:
    matcher.add("grch38", args.grch38)

lift = Forklift("grch37", "grch38")
dnm_df["grch38.locus"] = lift.liftover(dnm_df)

data = {
    sample: get_reads(
        dnm_df,
        bam_path,
        args.grch37 if args.reference == "grch37" else args.grch38,
        matcher,
        options=get_reads_opts,
    )
    for sample, bam_path in input_files.items()
}

fixed_rg_data = {}
n_reads = 0
for sample, d in data.items():
    hdr, reads = fix_readgroup(
        d[0], d[1], sample, n_reads
    )
    fixed_rg_data[sample] = (hdr, reads)
    n_reads += len(reads)


final_header = None
if len(fixed_rg_data) > 1:
    final_header = merge_readgroup(
        *map(lambda x: x[0], fixed_rg_data.values())
    )
else:
    final_header = fixed_rg_data[args.sample_id][0]

print(str(final_header), end="")
for read in chain(*map(lambda x: x[1], fixed_rg_data.values())):
    try:
        print(read.to_string().replace('\n', ''))
    except ValueError as err:
        logger.error("failed to marshal read <{}>: {}".format(str(read), err))