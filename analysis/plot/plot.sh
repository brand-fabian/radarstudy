#! /bin/bash

PYTHON="/home/brand/scratch/conda/envs/xcat_spark/bin/python3"
REFERENCE="$(realpath ~/scratch/library/human_g1k_v37_decoy.fasta)"
RADAR_META="$(realpath /home/brand/Projects/msdns_inova/data/radarstudy_data.csv)"
INOVA_META="$(realpath /home/brand/Projects/msdns_inova/data/Bonn-Subset-08-18-20.xlsx)"
CRU_META="$(realpath /home/brand/Projects/radarstudy/msdn_call/data/TRIO-CRU.md.csv)"
PILOT_META="$(realpath /home/brand/Projects/radarstudy/msdn_call/data/pilot.ped)"
PHASING="/ceph01/homedirs/brand/Projects/radarstudy/phasing/scripts/read_phasing/phasing.cru.pickle"
GRAPHTYPER_TABLES=("/ceph01/homedirs/brand/Projects/radarstudy/graphtyper/output/radar-vcf/radar-vcf.pickle" "/ceph01/homedirs/brand/Projects/radarstudy/graphtyper/output/inova/inova.pickle" "/ceph01/homedirs/brand/Projects/radarstudy/graphtyper/output/trio-cru/trio-cru.pickle")

DNM_PY="/ceph01/homedirs/brand/Projects/radarstudy/msdn_call/scripts/plot/dnm.py"
MSDN_PY="/ceph01/homedirs/brand/Projects/radarstudy/msdn_call/scripts/plot/msdn.py"
PPV_PY="/ceph01/homedirs/brand/Projects/radarstudy/msdn_call/scripts/plot/ppv.py"
PHASING_PY="/ceph01/homedirs/brand/Projects/radarstudy/msdn_call/scripts/plot/phasing.py"

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
  echo -e "\t--graphtyper\t\t\tInclude graphtyper tables from validation experiments"
  echo -e "\t--include-pilot\t\t\tInclude the pilot cohort in all analysis runs"
}

options=$(getopt -o hi:o:f:n:p: --long "input-prefix: help isolated metadata: output-dir: factor: apply-control-matching language: num-simulations: ppv: graphtyper include-pilot" -- "$@")
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
GRAPHTYPER=""
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
    --graphtyper)
      GRAPHTYPER="--graphtyper $(echo "${GRAPHTYPER_TABLES[@]}" | sed 's/ / --graphtyper /g') --apply-graphtyper-filter"
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

if [ -z "$METADATA_SET" -a -n "$INCLUDE_PILOT" ]; then
  METADATA="$METADATA $PILOT_METADATA"
fi

export CACHE_DIR="${OUTPUT_DIR}"

$PYTHON $DNM_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE $ISOLATED $GRAPHTYPER \
  "${INPUT_DIR}.refined_dnm.mt"

$PYTHON $MSDN_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE $GRAPHTYPER \
  "${INPUT_DIR}.msdns.ht"

$PYTHON $PPV_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE \
  $PPV $NUM_SIMULATIONS \
  "${INPUT_DIR}.msdns.ht"

$PYTHON $PHASING_PY --verbose info -R "$REFERENCE" \
  $METADATA $APPLY_CONTROL_MATCHING $FACTOR $LANGUAGE $GRAPHTYPER \
  -p "$PHASING" "${INPUT_DIR}.refined_dnm.mt" "${INPUT_DIR}.msdns.ht"
