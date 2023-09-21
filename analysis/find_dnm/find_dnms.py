#!/usr/bin/env python
# coding: utf-8

import hail as hl
import hail.expr.aggregators as agg
import sys
import os
import argparse
import logging
import pandas
import typing
from enum import Enum
from uuid import uuid4

##
# Logging
logger = logging.getLogger("find-msdns")
logger.setLevel(logging.ERROR)
channel = logging.StreamHandler()
formatter = logging.Formatter(
    "[%(asctime)s - %(name)s - %(levelname)7s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
channel.setFormatter(formatter)
logger.addHandler(channel)

log_level = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
}

##
# Checkpoint function (use write | read_matrix_table) instead
# of native hl.checkpoint to improve stability
class FILE_TYPE(Enum):
    MATRIX_TABLE = "MATRIX_TABLE"
    HAIL_TABLE = "HAIL_TABLE"
    DATAFRAME = "DATAFRAME"


def checkpoint(
    data: typing.Union[hl.MatrixTable, hl.Table, pandas.DataFrame],
    file_name: str,
    file_type: FILE_TYPE,
    overwrite: bool = False,
):
    path = os.path.abspath(file_name)
    if file_type == FILE_TYPE.MATRIX_TABLE:
        data.write(path, overwrite=overwrite)
        return hl.read_matrix_table(path)
    if file_type == FILE_TYPE.HAIL_TABLE:
        data.write(path, overwrite=overwrite)
        return hl.read_table(path)
    if file_type == FILE_TYPE.DATAFRAME:
        data.to_csv(path, sep="\t", index=False)
        return pandas.read_csv(path, sep="\t")


def vcf_input(input: str):
    if ":" in input:
        cohort, path = input.split(":", maxsplit=1)
        return (cohort, os.path.abspath(path))
    else:
        return (str(uuid4()), os.path.abspath(input))


##
# Parsing Arguments and sanitizing
#

checkpoints = set(
    [
        "data",
        "ddn",
        "r_de_novo",
        "refined_dnm",
        "msdns",
    ]
)

parser = argparse.ArgumentParser(
    prog="FIND-MSDNS",
    description="Find msdns in a given vcf file.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "vcf",
    metavar="VCF",
    type=vcf_input,
    nargs="+",
    help="VCF input file (prefer bgziped files) and optionally cohort name separated by colons. Example: radar:/path/to/vcf.bgz",
)
parser.add_argument(
    "-p", "--output", type=str, default="output", help="Output file prefix directory"
)
parser.add_argument(
    "-f", "--fam", type=str, required=True, help="Pedigree file (plink .fam format)"
)
parser.add_argument(
    "-R", "--reference", type=str, required=True, help="Reference sequence"
)
parser.add_argument(
    "--max-ac",
    type=int,
    default=10,
    help="Max. allele count in the cohort for a variant to be considered.",
)
parser.add_argument(
    "--min-p-de-novo",
    type=float,
    default=0.80,
    help="Min. probability of de novo event",
)
parser.add_argument(
    "--min-aaf",
    type=float,
    default=0.30,
    help="Min. alternate allele frequency to call het",
)
parser.add_argument(
    "--min-transversion-aaf",
    type=float,
    default=0.30,
    help="Min. alternate allele frequencies for transversions",
)
parser.add_argument(
    "--min-parent-dp", type=int, default=10, help="Min. sequencing depth in parents"
)
parser.add_argument(
    "--min-dp", type=int, default=15, help="Min. sequencing depth in index"
)
parser.add_argument(
    "--max-parent-alleles",
    type=int,
    default=1,
    help="Max. reads in parents supporting alt allele",
)
parser.add_argument("--window-size", type=int, default=20, help="MSDN window size")
parser.add_argument(
    "--checkpoint",
    type=str,
    action="append",
    default=None,
    choices=checkpoints,
    help="Enable checkpoints by name.",
)
parser.add_argument("--save-checkpoints", action="store_true")
parser.add_argument("--no-save-checkpoints", action="store_false")
parser.set_defaults(save_checkpoints=True)
parser.add_argument(
    "-t",
    "--tmp-dir",
    help="Temporary directory",
    type=str,
    default=os.getenv("TMPDIR", "/tmp"),
)
parser.add_argument(
    "--verbose", "-v", help="Set verbosity", choices=log_level.keys(), default="info"
)

args = parser.parse_args()
logger.setLevel(log_level[args.verbose])
logger.debug("Arguments: {}".format(vars(args)))

if not os.path.isfile(args.reference) or not os.path.isfile(args.reference + ".fai"):
    logger.error("Could not find reference file and index at {}".format(args.reference))
    sys.exit(1)

if any(not os.path.isfile(vcf[1]) for vcf in args.vcf):
    logger.error("Could not find vcf file at [{}]{}".format(args.vcf[0], args.vcf[1]))
    sys.exit(1)

if not os.path.isfile(args.fam):
    logger.error("Could not find pedigree file at {}".format(args.fam))
    sys.exit(1)

if args.checkpoint is not None:
    checkpoints = set(args.checkpoint)

output_prefix = os.path.abspath(args.output)

hl.init(
    master=f"spark://{os.getenv('SPARK_MASTER_IP', None)}:{os.getenv('SPARK_MASTER_PORT', None)}",
    idempotent=True,
    tmp_dir=os.getenv("SPARK_LOCAL_DIRS", os.getcwd()),
    local_tmpdir=os.getenv("SCRATCH_DIR", args.tmp_dir),
    log=os.path.join(os.getenv("SPARK_LOCAL_DIRS", os.getcwd()), "hail.log"),
    spark_conf={
        "spark.local.dir": os.getenv("SPARK_LOCAL_DIRS"),
        "spark.speculation": "true",
    },
)

rg = hl.get_reference("GRCh37")
rg.add_sequence(args.reference, "{}.fai".format(args.reference))

if len(args.vcf) > 1:
    logger.info("Merging {} vcfs using union_cols...".format(len(args.vcf)))

    partial_data_mt = map(
        lambda x: (
            x[0],
            hl.import_vcf(
                x[1],
                call_fields=["GT"],
                skip_invalid_loci=True,
                array_elements_required=False,
                find_replace=("-inf", "-100"),
            ),
        ),
        args.vcf,
    )

    data = None
    for cohort_name, mt in partial_data_mt:
        mt = mt.annotate_cols(cohort=cohort_name)
        if data is None:
            data = mt
        else:
            data = data.union_cols(mt, row_join_type="outer")
else:
    data = hl.import_vcf(
        args.vcf[0][1],
        call_fields=["GT"],
        skip_invalid_loci=True,
        array_elements_required=False,
        find_replace=("-inf", "-100"),
    )
    data.annotate_cols(cohort=args.vcf[0][0])

if args.save_checkpoints and "data" in checkpoints:
    logger.info("Writing raw data to {}.data.mt".format(output_prefix))
    data = checkpoint(
        data, "{}.data.mt".format(output_prefix), FILE_TYPE.MATRIX_TABLE, overwrite=True
    )

data = hl.sample_qc(data)
data = hl.variant_qc(data)

# Fix PL missingness values prior to split_multi_hts
# We can not set PL = hl.missing(data.PL.dtype), because hl.de_novo relies
# on PL values to compute the posterior prob. for de novo events.
data = data.annotate_entries(
    PL=hl.if_else(
        (hl.is_missing(data.PL)) & (hl.is_defined(data.GT)) & (data.GT.is_hom_ref()),
        hl.range(0, (hl.len(data.alleles) * (hl.len(data.alleles) + 1)) // 2).map(
            lambda x: hl.if_else(x == 0, 0, 1000)
        ),
        data.PL.map(lambda x: hl.if_else(hl.is_missing(x), 1000, x)),
    )
)

# Split multi hts for GLnexus based VCF files...
# Some loci occur twice in GLnexus based VCFs, which breaks split_multi_hts
multi_allelic = data.filter_rows(hl.len(data.alleles) > 2)
biallelic = data.filter_rows(hl.len(data.alleles) == 2)
biallelic = biallelic.annotate_rows(a_index=1, was_split=False)
split = hl.split_multi_hts(multi_allelic)
data = split.union_rows(biallelic)

# Annotate cohort level Allele Frequencies for the allele of each row
data = data.annotate_rows(
    iAF=data.variant_qc.AF[data.a_index],
    AC=data.variant_qc.AC[data.a_index],
    AN=data.variant_qc.AN,
)

if args.save_checkpoints and "data" in checkpoints:
    logger.info("Writing preprocessed data to {}.data2.mt".format(output_prefix))
    data = checkpoint(
        data,
        "{}.data2.mt".format(output_prefix),
        FILE_TYPE.MATRIX_TABLE,
        overwrite=True,
    )

# Read Pedigree and annotate de novo likelihoods
radar_ped = hl.Pedigree.read(args.fam)
de_novo_scores = hl.de_novo(data, radar_ped, pop_frequency_prior=data.iAF)
de_novo_mt = de_novo_scores.to_matrix_table(
    row_key=["locus", "alleles"],
    col_key=["id"],
    n_partitions=300,
)

de_novo_data = data.annotate_entries(
    p_de_novo=de_novo_mt[(data.locus, data.alleles), data.s].p_de_novo
)
if args.save_checkpoints and "ddn" in checkpoints:
    logger.info(
        "Writing data with de novo annotations to {}.ddn.mt".format(output_prefix)
    )
    de_novo_data = checkpoint(
        de_novo_data,
        "{}.ddn.mt".format(output_prefix),
        FILE_TYPE.MATRIX_TABLE,
        overwrite=True,
    )

# Filter DNM based on heuristic criteria and reduce size of the dataset
# by dropping most columns from the original matrix table
trio_mt = hl.trio_matrix(de_novo_data, radar_ped, complete_trios=True)

de_novo_data = de_novo_data.annotate_entries(
    mother=trio_mt[
        (de_novo_data.locus, de_novo_data.alleles), de_novo_data.s
    ].mother_entry,
    father=trio_mt[
        (de_novo_data.locus, de_novo_data.alleles), de_novo_data.s
    ].father_entry,
)

de_novo_data = de_novo_data.filter_entries(
    hl.is_defined(de_novo_data.GT)
    & de_novo_data.GT.is_non_ref()
    & (de_novo_data.p_de_novo > 0.7)
)

r_de_novo_mt = de_novo_data.select_cols("cohort")
r_de_novo_mt = r_de_novo_mt.select_rows("AC", "AN", "iAF")
r_de_novo_mt = r_de_novo_mt.select_entries(
    "AD",
    "DP",
    "GT",
    "p_de_novo",
    mother=hl.Struct(
        **{
            "AD": r_de_novo_mt.mother.AD,
            "AF": r_de_novo_mt.mother.AD[1] / hl.sum(r_de_novo_mt.mother.AD),
            "DP": r_de_novo_mt.mother.DP,
            "GT": r_de_novo_mt.mother.GT,
        }
    ),
    father=hl.Struct(
        **{
            "AD": r_de_novo_mt.father.AD,
            "AF": r_de_novo_mt.father.AD[1] / hl.sum(r_de_novo_mt.father.AD),
            "DP": r_de_novo_mt.father.DP,
            "GT": r_de_novo_mt.father.GT,
        }
    ),
    AF=r_de_novo_mt.AD[1] / hl.sum(r_de_novo_mt.AD),
)
if args.save_checkpoints and "r_de_novo" in checkpoints:
    logger.info("Writing raw de novo calls to {}.r_de_novo.mt".format(output_prefix))
    r_de_novo_mt = checkpoint(
        r_de_novo_mt,
        "{}.r_de_novo.mt".format(output_prefix),
        FILE_TYPE.MATRIX_TABLE,
        overwrite=True,
    )

# Further refine de novo mutations based on criteria such as cohort AC/AF, DP and others.
refined_dnm = r_de_novo_mt.filter_rows(
    (r_de_novo_mt.alleles[0].length() == 1)
    & (r_de_novo_mt.alleles[1].length() == 1)
    & (r_de_novo_mt.locus.in_autosome())
)
refined_dnm = refined_dnm.filter_entries(
    (refined_dnm.AC <= args.max_ac)
    & (refined_dnm.p_de_novo >= args.min_p_de_novo)
    & (refined_dnm.DP >= args.min_dp)
    & (refined_dnm.mother.DP >= args.min_parent_dp)
    & (refined_dnm.father.DP >= args.min_parent_dp)
    & (refined_dnm.AF >= args.min_aaf)
    & (refined_dnm.AF <= 1 - args.min_aaf)
    & (refined_dnm.father.AD[1] < args.max_parent_alleles + 1)
    & (refined_dnm.mother.AD[1] < args.max_parent_alleles + 1)
    & (
        (
            (hl.is_transversion(refined_dnm.alleles[0], refined_dnm.alleles[1]))
            & (refined_dnm.AF > args.min_transversion_aaf)
            & (refined_dnm.AF < 1 - args.min_transversion_aaf)
        )
        | ((~hl.is_transversion(refined_dnm.alleles[0], refined_dnm.alleles[1])))
    )
)

if args.save_checkpoints and "refined_dnm" in checkpoints:
    logger.info(
        "Writing refined de novo calls to {}.refined_dnm.mt".format(output_prefix)
    )
    refined_dnm = checkpoint(
        refined_dnm,
        "{}.refined_dnm.mt".format(output_prefix),
        FILE_TYPE.MATRIX_TABLE,
        overwrite=True,
    )

dnm_ht = refined_dnm.entries()
logger.info("Got {} de novo calls".format(dnm_ht.count()))

logger.info("Writing refined de novo to {}.dnm.tsv".format(output_prefix))
dnm_ht.to_pandas().to_csv("{}.dnm.tsv".format(output_prefix), index=False, sep="\t")

##
# hl.methods.window_by_locus was removed after v0.2.11, and other alternatives
# reproducibly result in errors (spark executors gettings stuck).
# Therefore, we just use pandas for the MSDN calling step...
def find_msdns(tbl, window_size=20):
    msdns = []
    msdn_id = 1
    previous_sampled = False
    for sample, group in tbl.groupby("s"):
        if previous_sampled:
            # go to next msdn if the last lesions from the previous sample
            # are a msdn
            previous_sampled = False
            msdn_id += 1
        previous = None
        prev_contig = None
        prev_position = None
        for _, row in group.iterrows():
            if previous is not None:
                if abs(row.locus.position - previous.locus.position) < window_size + 1:
                    previous_sampled = True
                    msdns.append(
                        pandas.concat(
                            [
                                row,
                                previous.add_prefix("previous."),
                                pandas.Series([msdn_id], index=["msdn_id"]),
                            ]
                        )
                    )
                elif previous_sampled:
                    previous_sampled = False
                    msdn_id += 1

            previous = row
    return pandas.DataFrame(
        msdns,
        columns=list(tbl.columns)
        + list(map(lambda s: f"previous.{s}", tbl.columns))
        + ["msdn_id"],
    )


def to_hail_table(msdn_tbl):
    msdn_ht = hl.Table.from_pandas(msdn_tbl)
    msdn_ht = msdn_ht.annotate(
        mother=hl.Struct(
            **{key: msdn_ht[f"mother.{key}"] for key in ["AD", "AF", "DP", "GT"]}
        ),
        father=hl.Struct(
            **{key: msdn_ht[f"father.{key}"] for key in ["AD", "AF", "DP", "GT"]}
        ),
        previous=hl.Struct(
            **{
                **{
                    key: msdn_ht[f"previous.{key}"]
                    for key in [
                        "locus",
                        "alleles",
                        "s",
                        "AC",
                        "AN",
                        "iAF",
                        "AD",
                        "DP",
                        "GT",
                        "p_de_novo",
                        "AF",
                        "cohort",
                    ]
                },
                **{
                    "mother": hl.Struct(
                        **{
                            key: msdn_ht[f"previous.mother.{key}"]
                            for key in ["AD", "AF", "DP", "GT"]
                        }
                    ),
                    "father": hl.Struct(
                        **{
                            key: msdn_ht[f"previous.father.{key}"]
                            for key in ["AD", "AF", "DP", "GT"]
                        }
                    ),
                },
            }
        ),
    )
    msdn_ht = msdn_ht.select(
        "locus",
        "alleles",
        "s",
        "msdn_id",
        "AC",
        "AN",
        "iAF",
        "AD",
        "DP",
        "GT",
        "p_de_novo",
        "AF",
        "mother",
        "father",
        "previous",
        "cohort",
    )
    return msdn_ht


tbl = refined_dnm.entries().to_pandas()
msdn_tbl = find_msdns(tbl, window_size=args.window_size)

logger.info("Writing msdn calls to {}.msdns.tsv".format(output_prefix))
msdn_tbl.to_csv("{}.msdns.tsv".format(output_prefix), sep="\t", index=False)

msdn_ht = to_hail_table(msdn_tbl)

if args.save_checkpoints and "msdns" in checkpoints:
    logger.info("Writing msdn calls to {}.msdns.ht".format(output_prefix))
    msdn_ht = checkpoint(
        msdn_ht,
        "{}.msdns.ht".format(output_prefix),
        FILE_TYPE.HAIL_TABLE,
        overwrite=True,
    )
