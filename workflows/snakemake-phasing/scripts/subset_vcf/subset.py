import sys
import os
import argparse
import logging
import pysam
from collections import namedtuple
import typing
from liftover import get_lifter


# Header Fields
INFO = {
    'OCONTIG': [
        ('ID', 'OCONTIG'),
        ('Type', 'String'),
        ('Number', '1'),
        ('Description', 'Original contig before applying liftover'),
    ],
    'OPOS': [
        ('ID', 'OPOS'),
        ('Type', 'String'),
        ('Number', '1'),
        ('Description', 'Original position before applying liftover'),
    ]
}


###############################################################################
# Logging                                                                     #
###############################################################################
logger = logging.getLogger("subset-vcf")
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
LIFTER_REFERENCES = {
    "grch37": "hg19",
    "hg19": "hg19",
    "grch38": "hg38",
    "hg38": "hg38",
}

parser = argparse.ArgumentParser(
    prog="SUBSET-VCF",
    description='Subset VCF files based on a given regions file and liftover coordinates.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('vcf', metavar='VCF', type=str,
                    help='VCF file to subset and transcode')
parser.add_argument('-S', '--sites', type=str, default=None,
                    help='.tsv Sites file to include in the final vcf file. Each line should include the three columns contig, start and end')
parser.add_argument('-t', '--tmp-dir', help='Temporary directory',
                    type=str, default=os.getenv('TMPDIR', '/tmp'))
parser.add_argument('--source-reference', type=str, choices=LIFTER_REFERENCES.keys(),
                    default='grch37', help='Source reference for liftover')
parser.add_argument('--target-reference', type=str, choices=LIFTER_REFERENCES.keys(),
                    default='grch38', help='Target reference for liftover')
parser.add_argument('-R', '--target-reference-fasta', type=str,
                    help='Path to the target reference fasta to retrieve contigs')
parser.add_argument('-l', '--apply-liftover', action='store_true',
                    default=False, help='Liftover variant coordinates from source to target reference.')
parser.add_argument('-o', '--output', default=None, type=str,
                    help='VCF output path (default = stdout)')
parser.add_argument('--verbose', '-v', help='Set verbosity',
                    choices=log_level.keys(), default='info')

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))


if not os.path.isfile(args.vcf):
    logger.error("VCF file {} does not exist".format(args.vcf))
    sys.exit(1)

if args.sites is not None and not os.path.isfile(args.sites):
    logger.error("Sites file {} does not exist".format(args.sites))
    sys.exit(1)


Region = namedtuple("Region", ["contig", "start", "end"])
Coordinates = namedtuple("Coordinates", ["contig", "pos"])


def target_regions(regions_file: str) -> typing.Generator[Region, None, None]:
    with open(regions_file, "r") as r_f:
        for line in r_f:
            contig, start, end = line.strip().split("\t")
            yield Region(contig=contig, start=int(start), end=int(end))


def liftover_record(lifter, record: pysam.VariantRecord) -> typing.Optional[Coordinates]:
    try:
        contig, pos, _ = lifter.query(record.contig, record.pos)[0]
        return Coordinates(contig=contig, pos=pos)
    except KeyError as e:
        logger.warning("failed liftover for {}:{}".format(record.contig, record.pos))
        return None
    except IndexError as e:
        logger.warning("no location found for {}:{}".format(record.contig, record.pos))
        return None


def get_sites(vcf_in: pysam.VariantFile, sites: typing.Optional[str] = args.sites):
    if sites is None:
        for record in vcf_in.fetch():
            yield record
    else:
        for region in target_regions(args.sites):
            for record in vcf_in.fetch(region.contig, region.start, region.end):
                yield record


class Writer():
    def __init__(
        self,
        output: typing.Optional[str],
        header: pysam.VariantHeader,
        apply_liftover = False,
        source_reference: str = 'grch37',
        target_reference: str = 'grch38',
        target_reference_fasta: typing.Optional[str] = None,
    ):
        self._output = output
        self._apply_liftover = apply_liftover
        self._lifter = get_lifter(
            LIFTER_REFERENCES[source_reference],
            LIFTER_REFERENCES[target_reference]
        )
        self._header = header
        if apply_liftover:
            for val in INFO.values():
                self._header.add_meta("INFO", items=val)
            if target_reference_fasta is not None:
                self._header.contigs.clear_header()
                fasta = pysam.FastaFile(target_reference_fasta)
                try:
                    for contig, length in zip(fasta.references, fasta.lengths):
                        self._header.contigs.add(contig, length=length)
                finally:
                    fasta.close()
        self._vcf_out = None if output is None else pysam.VariantFile(output, 'w', header=header)
    
    def write(self, record: pysam.VariantRecord):
        if self._apply_liftover:
            coords = liftover_record(self._lifter, record)
            if coords is not None:
                record.info['OCONTIG'] = record.contig
                record.info['OPOS'] = str(record.pos)
                record.contig = coords.contig
                record.pos = coords.pos
            else:
                return # ignore unfound records

        if self._vcf_out is None:
            print(str(record), end='')
        else:
            self._vcf_out.write(record)


vcf_in = pysam.VariantFile(args.vcf)
output_wr = Writer(
    args.output,
    vcf_in.header,
    apply_liftover=args.apply_liftover,
    source_reference=args.source_reference,
    target_reference=args.target_reference,
    target_reference_fasta=args.target_reference_fasta,
)
try:
    for record in get_sites(vcf_in):
        output_wr.write(record)
finally:
    vcf_in.close()
