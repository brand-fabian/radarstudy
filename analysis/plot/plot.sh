#! /bin/bash

#############
# Variables #
#############
# These variables are the main settings that are mostly fixed between
# runs of the statistical analysis pipelines. Please set them according
# to your local environment.
PYTHON="<PYTHON INTERPRETER PATH>"
REFERENCE="<GENOME REFERENCE .fasta PATH>"
# Metadata Files for the different cohorts
RADAR_META="<RADARSTUDY METADATA .csv>"
INOVA_META="<INOVA METADATA .xlsx"
CRU_META="<CRU METADATA .csv>"
PILOT_META="<PILOT METADATA .csv>"
# If phasing is to be incorporated into the analysis, set the path to
# the output table of the phasing workflow.
PHASING="<PHASING TABLE>"

# Paths to the analysis scripts in this folder, change these if you
# changed the directory layout
DNM_PY="$(realpath $PWD/dnm.py)"
MSDN_PY="$(realpath $PWD/msdn.py)"
PPV_PY="$(realpath $PWD/ppv.py)"
PHASING_PY="$(realpath $PWD/phasing.py)"

#############
# Arguments #
#############
# Use these arguments to change each script execution for the analysis
# that you want to execute through the CLI.
function usage {
  echo -e "$0\n"
  echo -e "Shortcut to execute the radar scripts in this folder on some input data.\n"
  echo -e "Options:"
  echo -e "\t-i, --input-prefix\t\tPrefix/Folder to search input data in"
  echo -e "\t--metadata\t\t\tMetadata pickle file to use instead of computing from default metadata"
  echo -e "\t-o, --output-dir\t\tOutput directory and file prefix to export results to"
  echo -e "\t-f, --factor\t\t\tDownsampling factor to apply to match n controls per case (default = 0 = no matching)"
  echo -e "\t--apply-control-matching\tSet the apply control matching flag for plotting if this is set"
  echo -e "\t--language\t\t\tSet the language for plot descriptions (choices: de, en)"
  echo -e "\t--isolated\t\t\tCompute isolated de-novo mutations"
  echo -e "\t--include-pilot\t\t\tInclude the pilot cohort in all analysis runs"
}

options=$(getopt -o hi:o:f:n:p: --long "input-prefix: help isolated metadata: output-dir: factor: apply-control-matching language: num-simulations: ppv: include-pilot" -- "$@")
[ $? -eq 0 ] || {
  usage
  echo -e "\nInvalid options provided"
  exit 1
}

INPUT_DIR="$(realpath .)"
OUTPUT_DIR="$(realpath .)"
METADATA="--radar-meta $RADAR_META --inova-meta $INOVA_META --cru-meta $CRU_META"
PILOT_METADATA="--pilot-meta $PILOT_META"
FACTOR="--factor 0"
APPLY_CONTROL_MATCHING=""
LANGUAGE="--language de"
PPV=""
NUM_SIMULATIONS="--num-simulations 1000"
ISOLATED=""
INCLUDE_PILOT=""
METADATA_SET=""

eval set -- "$options"
while true; do
  case "$1" in
    -i | --input-prefix)
      shift
      INPUT_DIR="$(realpath $1)"
      ;;
    -o | --output-dir)
      shift
      if [ ! -d $1 ]; then
        mkdir -p $1
      fi
      OUTPUT_DIR="$(realpath $1)"
      ;;
    --metadata)
      shift
      METADATA="--metadata $(realpath $1)"
      METADATA_SET="1"
      ;;
    -f | --factor)
      shift
      FACTOR="--factor $1"
      ;;
    --apply-control-matching)
      APPLY_CONTROL_MATCHING="--apply-control-matching"
      ;;
    --language)
      shift
      LANGUAGE="--language $1"
      ;;
    --num-simulations | -n)
      shift
      NUM_SIMULATIONS="--num-simulations $1"
      ;;
    --ppv | -p)
      shift
      PPV="--ppv $1"
      ;;
    --isolated)
      ISOLATED="--isolated --msdn-ht ${INPUT_DIR}.msdns.ht"
      ;;
    --include-pilot)
      INCLUDE_PILOT="1"
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
  esac
  shift
done

###############
# Main Script #
###############
if [ -z "$METADATA_SET" -a -n "$INCLUDE_PILOT" ]; then
  METADATA="$METADATA $PILOT_METADATA"
fi

# The cache directory environment is used by all scripts for their
# output, since this directory has to be known when the pyhton files
# are first interpreted.
export CACHE_DIR="${OUTPUT_DIR}"

$PYTHON $DNM_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE $ISOLATED \
  "${INPUT_DIR}.refined_dnm.mt"

$PYTHON $MSDN_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE \
  "${INPUT_DIR}.msdns.ht"

$PYTHON $PPV_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE \
  $PPV $NUM_SIMULATIONS \
  "${INPUT_DIR}.msdns.ht"

$PYTHON $PHASING_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE \
  -p "$PHASING" "${INPUT_DIR}.refined_dnm.mt" "${INPUT_DIR}.msdns.ht"
