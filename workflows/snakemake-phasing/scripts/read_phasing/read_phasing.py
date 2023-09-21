import argparse
import typing
import pandas
import pysam
import os
import sys
import logging
import pickle
from textwrap import dedent
from multiprocessing import Pool, cpu_count


LOG_INTERVAL_ROWS = 1000


###############################################################################
# Logging                                                                     #
###############################################################################
logger = logging.getLogger("read-phasing")
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
def write_pickle(df: pandas.DataFrame, x: str):
    with open(x, 'wb') as pickle_f:
        pickle.dump(df, pickle_f)

output_format = {
    'pickle': write_pickle,
    'csv': lambda df, x: df.to_csv(x, index=False),
}

class TextAndDefaultsFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    prog="READ-PHASING",
    description=dedent('''\
        Read phase information from a set of VCF files.

        Read all vcf files and turn them into a pandas table that can be loaded
        by subsequent tools (e.g. hail). The resulting table will specify the
        columns ["locus", "alleles", "sample"] for indexing and the data columns
        ["is_phased", "is_paternal"] displaying the phase information.

        Note: We assume only biallelic calls in the dataset (e.g. 0/1 GT's)'''),
    formatter_class=TextAndDefaultsFormatter)
parser.add_argument('vcf', metavar='VCF', type=str, nargs='+',
                    help='Phasing VCF to read data from.')
parser.add_argument('-S', '--sites', type=str, required=False, default=None,
                    help='Sites to include in the final table')
parser.add_argument('-O', '--output-format', type=str, default="pickle",
                    choices=output_format.keys(), help='Path to write output table to')
parser.add_argument('-o', '--output', type=str, required=True,
                    help='Path to write output table to')
parser.add_argument('-t', '--tmp-dir', help='Temporary directory',
                    type=str, default=os.getenv('TMPDIR', '/tmp'))
parser.add_argument('--threads', default=cpu_count() - 1, type=int,
                    help="Number of reader threads to use.")
parser.add_argument('--verbose', '-v', help='Set verbosity',
                    choices=log_level.keys(), default='info')

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))


###############################################################################
# Main script                                                                 #
###############################################################################
class VariantReader:
    def __init__(self, header=['locus', 'alleles', 'sample', 
                               'is_phased', 'is_paternal']):
        self._data = []
        self._header = header

    def to_dataframe(self):
        return pandas.DataFrame(self._data, columns=self._header)

    def filter_by_sites(self, sites: typing.Optional[typing.Set]):
        if sites is None:
            return self

        new_obj = self.__class__()
        new_obj._data = filter(
            lambda row: row[0] in sites,
            self._data
        )
        return new_obj

    @classmethod
    def from_vcf(cls, path: str):
        reader = cls()
        logger.info("reading vcf at {}".format(path))
        vcf = pysam.VariantFile(path)
        rows = 0
        for i, record in enumerate(vcf.fetch()):
            locus = "{}:{}".format(str(record.chrom), str(record.start+1))
            for sample, record_sample in record.samples.items():
                gt = None
                if 'GT' in record_sample and all(x is not None for x in record_sample['GT']):
                    gt = record_sample['GT']
                elif 'GT' in record_sample and any(x is None for x in record_sample['GT']):
                    logger.debug("ignoring no_call for sample {} at {}".format(
                        sample, locus
                    ))
                else:
                    raise Exception("missing gt in row {} for sample {}".format(
                        i, sample
                    ))

                if gt is not None:
                    alleles = list(record_sample.alleles)
                    if gt[0] == 1 and gt[1] == 0:
                        alleles = alleles[::-1]
                        logger.debug("paternally phased mutation - putting reference first")

                    reader._data.append([
                        locus,
                        alleles,
                        sample,
                        record_sample.phased,
                        record_sample.phased and gt is not None and gt[0] == 1,
                    ])

            if i > 0 and i % LOG_INTERVAL_ROWS == 0:
                logger.debug("read {} vcf rows".format(i))
            rows = i
        logger.debug("read {} vcf rows".format(rows))
        vcf.close()
        return reader

if not all(os.path.isfile(x) for x in args.vcf):
    logger.error("missing vcf input file {}".format(args.vcf))
    sys.exit(1)

if not os.path.isdir(os.path.dirname(args.output)):
    os.makedirs(os.path.abspath(os.path.dirname(args.output)), exist_ok=True)

sites = None
if args.sites is not None:
    if not os.path.isfile(args.sites):
        logger.error("missing sites file {}".format(args.sites))
        sys.exit(1)

    sites_df = pandas.read_csv(args.sites, sep="\t")
    sites = set(sites_df.locus.map(str))
    logger.info("filtering for {} sites".format(len(sites)))


in_paths = map(os.path.abspath, args.vcf)

def _read_vcf(path: str, sites: typing.Optional[typing.Set] = sites) -> pandas.DataFrame:
    return VariantReader.from_vcf(path).filter_by_sites(sites).to_dataframe()

logger.info("reading phase information from {} vcfs".format(len(args.vcf)))

with Pool(args.threads) as pool:
    vcf_df = pandas.concat(pool.map(_read_vcf, in_paths))

logger.info("writing output .{} to {}".format(
    args.output_format, os.path.abspath(args.output),
))
output_format[args.output_format](vcf_df, os.path.abspath(args.output))