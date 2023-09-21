import argparse
import json
import logging
import os
import pickle
import subprocess
import sys
import typing
from datetime import datetime
from contextlib import contextmanager
import traceback

import hail as hl
import pandas

from meta.config import LOG_LEVEL
from meta.loaders import MsdnLoader

logger = logging.getLogger(__name__)


def configure_logger():
    """Customize the base logger instance for use with the analysis scripts."""
    logger = logging.getLogger()
    logger.setLevel(logging.ERROR)
    channel = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)7s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    channel.setFormatter(formatter)
    logger.addHandler(channel)
    return logger


def add_base_arguments(parser: argparse.ArgumentParser):
    """Add common arguments for analysis scripts to the parser."""
    parser.add_argument(
        "-R", "--reference", type=str, required=True, help="Reference sequence"
    )
    parser.add_argument(
        "--radar-meta", type=str, default=None, help="Path to radarstudy metadata table"
    )
    parser.add_argument(
        "--inova-meta", type=str, default=None, help="Path to inova metadata table"
    )
    parser.add_argument(
        "--cru-meta", type=str, default=None, help="Path to trio-cru metadata table"
    )
    parser.add_argument(
        "--pilot-meta", type=str, default=None, help="Path to the pilot study metadata table"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="Path to a pandas dataframe .pickle file containing study metadata",
    )
    parser.add_argument(
        "-f",
        "--factor",
        type=int,
        default=0,
        help="Control cohort downsampling factor. If set to 0 (the default), no matching is performed.",
    )
    parser.add_argument(
        "--language",
        type=str,
        default="de",
        choices=["de", "en"],
        help="Language of axis descriptions in generated plots",
    )
    parser.add_argument(
        "--graphtyper",
        type=str,
        default=None,
        action="append",
        help="Path to pickle files with graphtyper information for the validation of DNMs.",
    )
    parser.add_argument(
        "--apply-graphtyper-filter",
        action="store_true",
        dest="apply_graphtyper_filter",
        help="If this flag is set, variants not confirmed by graphtyper will be ignored.",
    )
    parser.add_argument(
        "--no-apply-graphtyper-filter",
        action="store_true",
        dest="apply_graphtyper_filter",
        help="If this flag is set, variants not confirmed by graphtyper will be ignored.",
    )
    parser.add_argument(
        "--apply-control-matching",
        action="store_true",
        dest="apply_control_matching",
        help="If this flag is set the control matching will split the control cohort for some plots",
    )
    parser.add_argument(
        "--no-apply-control-matching",
        action="store_false",
        dest="apply_control_matching",
        help="If this flag is set the control matching will split the control cohort for some plots",
    )
    parser.add_argument("--load-cache", action="store_true", dest="load_cache")
    parser.add_argument("--no-load-cache", action="store_false", dest="load_cache")
    parser.set_defaults(
        load_cache=True, apply_control_matching=False, apply_graphtyper_filter=False
    )
    parser.add_argument(
        "-t",
        "--tmp-dir",
        help="Temporary directory",
        type=str,
        default=os.getenv("TMPDIR", "/tmp"),
    )
    parser.add_argument(
        "--verbose",
        "-v",
        help="Set verbosity",
        choices=LOG_LEVEL.keys(),
        default="info",
    )
    return parser


def hail_init(args: argparse.Namespace):
    hl.init(
        master=f"spark://{os.getenv('SPARK_MASTER_IP', None)}:{os.getenv('SPARK_MASTER_PORT', None)}",
        idempotent=True,
        tmp_dir=os.getenv("SPARK_LOCAL_DIRS"),
        local_tmpdir=os.getenv("SCRATCH_DIR", "/tmp"),
        log=os.path.abspath(
            os.path.join(os.getenv("SPARK_LOCAL_DIRS", "."), "hail.log")
        ),
        spark_conf={
            "spark.local.dir": os.getenv("SPARK_LOCAL_DIRS"),
            "spark.speculation": "true",
        },
    )

    rg = hl.get_reference("GRCh37")
    rg.add_sequence(
        os.path.abspath(args.reference),
        "{}.fai".format(os.path.abspath(args.reference)),
    )


def load_metadata_from_file(args: argparse.Namespace):
    if (
        args.inova_meta is not None
        and args.radar_meta is not None
        and args.cru_meta is not None
    ):
        if not os.path.isfile(args.inova_meta):
            logger.error("Inova metadata not found at {}".format(args.inova_meta))
            sys.exit(1)

        if not os.path.isfile(args.radar_meta):
            logger.error("Radarstudy metadata not found at {}".format(args.radar_meta))
            sys.exit(1)

        if not os.path.isfile(args.cru_meta):
            logger.error("CRU metadata not found at {}".format(args.cru_meta))
            sys.exit(1)
    elif args.metadata is not None:
        if not os.path.isfile(args.metadata):
            logger.error("Study metadata not found at {}".format(args.metadata))
            sys.exit(1)
    else:
        logger.error(
            "Study metadata missing. Please supply either [[--radar-meta FILE] and [--inova-meta FILE] and [--cru-meta FILE]] or [--metadata FILE]."
        )
        sys.exit(1)

    if args.pilot_meta is not None and not os.path.isfile(args.pilot_meta):
        logger.error("Pilot metadata not found at {}".format(args.pilot_meta))
        sys.exit(1)

    age_data: typing.Optional[pandas.DataFrame] = None
    if (
        args.inova_meta is not None
        and args.radar_meta is not None
        and args.cru_meta is not None
    ):
        age_data = MsdnLoader.load_metadata(
            args.radar_meta, args.inova_meta, args.cru_meta, pilot=args.pilot_meta,
        )
    elif args.metadata is not None:
        with open(os.path.abspath(args.metadata), "rb") as pickle_f:
            age_data = pickle.load(pickle_f)

    if age_data is None:
        logger.error("missing metadata for analysis")
        sys.exit(1)

    return age_data


def get_call_information(
    args: argparse.Namespace,
):
    """Get call information for this invocation of the tool, including
       git revision and any call arguments.
    """
    git_revision = (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return {
        "args": vars(args),
        "time": date_str,
        "revision": git_revision,
        "script": os.path.abspath(sys.argv[0]),
    }


@contextmanager
def dump_call_information(
    output: str,
    args: argparse.Namespace,
):
    """Write information about the script invocation to a file."""
    current_info = {
        "arguments": {},
        "return": {},
        "last_update": None,
    }
    if os.path.isfile(output):
        try:
            with open(output, "r") as in_f:
                current_info = json.load(in_f)
        except Exception as e:
            logger.warning("Could not read call information file at {}: {}".format(output, e))

    call_args = get_call_information(args)
    try:
        yield call_args
        current_info = {
            "arguments": {
                **current_info["arguments"],
                call_args["time"]: call_args,
            },
            "return": {
                **current_info["return"],
                call_args["time"]: "SUCCESS",
            },
            "last_update": call_args["time"],
        }
    except Exception as e:
        logger.error("Script execution failed: {}".format(str(e)))
        logger.error(traceback.format_exc())
        current_info = {
            "arguments": {
                **current_info["arguments"],
                call_args["time"]: call_args,
            },
            "return": {
                **current_info["return"],
                call_args["time"]: "FAILED ({})".format(str(e)),
            },
            "last_update": call_args["time"],
        }
    finally:
        with open(output, "w") as out_f:
            json.dump(current_info, out_f, indent=4)
