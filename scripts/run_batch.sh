#!/bin/bash

# File: scripts/run_batch.sh
# Usage: ./scripts/run_batch.sh path/to/dataset [--svg] [--netdeg] [--negcycles] [--negcyclescover] [--triangles] [--filter pattern] [--dry-run] [--no-build]

WITH_SVG=""
WITH_NETDEG=""
WITH_NEGCYCLES=""
WITH_NEGCYCLESCOVER=""
WITH_TRIANGLES=""
MODEL_TYPE="xy"
ROOT_DIR=""
DRY_RUN=false
NO_BUILD=false

FILTER=""
for arg in "$@"; do
  case "$arg" in
    --svg)
      WITH_SVG="--svg";;
    --netdeg)
      WITH_NETDEG="--netdeg";;
    --negcycles)
      WITH_NEGCYCLES="--negcycles";;
    --negcyclescover)
      WITH_NEGCYCLESCOVER="--negcyclescover";;
    --triangles)
      WITH_TRIANGLES="--triangles";;
    --model=*)
      MODEL_TYPE="${arg#*=}";;
    --filter)
      shift
      FILTER="$1";;
    --dry-run)
      DRY_RUN=true;;
    --no-build)
      NO_BUILD=true;;
    --*)
      echo "[Warning] Unrecognized option: $arg";;
    *)
      if [ -z "$ROOT_DIR" ]; then
        ROOT_DIR="$arg"
      fi;;
  esac
done

# Validate input
if [ $# -eq 0 ] || [[ "$ROOT_DIR" == --* ]] || [ ! -d "$ROOT_DIR" ]; then
  echo "Usage: $0 <root_dataset_dir> [--svg] [--netdeg] [--negcycles] [--negcyclescover] [--triangles] [--filter pattern] [--dry-run] [--no-build]"
  echo "--netdeg Enable net-degree cut inequalities"
  echo "--negcycles Enable negative-cycle cut inequalities"
  echo "--negcyclescover Enable negative cycle inequalities as a cover"
  echo "--triangles Enable all negative triangles cut inequalities"
  echo "--svg Export graph solution as SVG visualization"
  echo "--dry-run to list matching CSVs without executing anything"
  echo "--no-build to skip the build step and directly run frustration_index"
  echo "Example: $0 /data/dataset --netdeg --negcycles --filter pattern"
  exit 1
fi

SCRIPT_START=$(date +%s)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="logs/$TIMESTAMP"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/results_summary.txt"
RES_DIR="results/$TIMESTAMP"
mkdir -p "$RES_DIR"
SKIP_ALL=false

CSV_COUNT=$(find "$ROOT_DIR" -type f -name "*.csv" | wc -l)
if [ "$CSV_COUNT" -eq 0 ]; then
  echo "[Error] No .csv files found in $ROOT_DIR."
  exit 1
fi

PROJECT_ROOT=$(realpath "$(dirname "$0")/..")

if [ "$NO_BUILD" = false ] && [ ! -x ./frustration_index ]; then
  echo "[Info] Building executable..."
  make -C "$PROJECT_ROOT/Debug/build" -j || exit 1
fi

if [ -n "$FILTER" ]; then
  find "$ROOT_DIR" -type f -name "*${FILTER}*.csv" > "$LOG_FILE.tmp"
else
  find "$ROOT_DIR" -type f -name "*.csv" > "$LOG_FILE.tmp"
fi

if [ "$DRY_RUN" = true ]; then
  echo "[Dry Run] Matching files:"
  cat "$LOG_FILE.tmp"
  rm "$LOG_FILE.tmp"
  exit 0
fi

echo "Dataset,RelativePath,Objective,Status,Time(s),Cuts" > "$LOG_FILE"
TOTAL=$(wc -l < "$LOG_FILE.tmp")
INDEX=0

INPUT_FILES=()

while read -r FILE; do
  INDEX=$((INDEX + 1))
  echo "[Progress] File $INDEX of $TOTAL"
  BASENAME=$(basename "$FILE")
  RELPATH=$(realpath --relative-to="$ROOT_DIR" "$FILE")
  DIRNAME=$(dirname "$RELPATH")
  INPUT_FILES+=("$RELPATH")

  LOG_PATH="logs/${RELPATH%.csv}.log"
  SUMMARY_PATH="logs/$DIRNAME/summary.csv"

  if [[ -f "$LOG_PATH" || -f "$SUMMARY_PATH" ]]; then
    if [ "$SKIP_ALL" = false ]; then
      echo -n "[Skip?] $RELPATH already processed. Skip? [y/N] (Y for all): "
      read -r yn
      if [[ "$yn" =~ ^[Yy]$ ]]; then
        echo "[Skipped] $RELPATH"
        continue
      elif [[ "$yn" =~ ^Y$ ]]; then
        SKIP_ALL=true
        echo "[Skipped All] All remaining files will be skipped if already processed."
        continue
      fi
    else
      echo "[Skipped] $RELPATH"
      continue
    fi
  fi

  START=$(date +%s.%N)
  OUTPUT=$("$PROJECT_ROOT/Debug/build/frustration_index" "$FILE" --model=$MODEL_TYPE $WITH_NETDEG $WITH_NEGCYCLES $WITH_NEGCYCLESCOVER $WITH_TRIANGLES $WITH_SVG 2>&1)
  STATUS=$?
  END=$(date +%s.%N)
  DURATION=$(echo "$END - $START" | bc)

  OBJECTIVE=$(echo "$OUTPUT" | grep "Objective" | awk '{print $NF}')
  if [[ -n "$WITH_NETDEG" || -n "$WITH_NEGCYCLES" || -n "$WITH_NEGCYCLESCOVER" || -n "$WITH_TRIANGLES" ]]; then
    CUT_FLAG="yes"
  else
    CUT_FLAG="no"
  fi
  echo "$BASENAME,$RELPATH,$OBJECTIVE,$STATUS,$DURATION,$CUT_FLAG" >> "$LOG_FILE"

  mkdir -p "$LOG_DIR/$DIRNAME"
  echo "$OUTPUT" > "$LOG_DIR/${RELPATH%.csv}.log"
  SUMMARY_FILE="$LOG_DIR/$DIRNAME/summary.csv"
  if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Dataset,RelativePath,Objective,Status,Time(s)" > "$SUMMARY_FILE"
  fi
  echo "$BASENAME,$RELPATH,$OBJECTIVE,$STATUS,$DURATION" >> "$SUMMARY_FILE"

  echo "[Done] $RELPATH -> Obj: $OBJECTIVE, Time: ${DURATION}s"

done < "$LOG_FILE.tmp"

rm "$LOG_FILE.tmp"

# Build file suffix used for exported CSVs
SUFFIX=""
[ -n "$WITH_NETDEG" ] && SUFFIX+="_netdeg"
[ -n "$WITH_NEGCYCLES" ] && SUFFIX+="_negcycles"
[ -n "$WITH_NEGCYCLESCOVER" ] && SUFFIX+="_negcyclescover"
[ -n "$WITH_TRIANGLES" ] && SUFFIX+="_triangles"
[ -z "$SUFFIX" ] && SUFFIX="_nocuts"
[ "$MODEL_TYPE" != "xy" ] && SUFFIX+="_${MODEL_TYPE}"

# Move only specific output files to their corresponding subdirectories
for RELPATH in "${INPUT_FILES[@]}"; do
  DIRNAME=$(dirname "$RELPATH")
  BASENAME=$(basename "$RELPATH" .csv)
  for EXT in _summary.csv _x.csv _y.csv; do
    FILE="${BASENAME}${SUFFIX}${EXT}"
    if [ -f "$FILE" ]; then
      mkdir -p "$RES_DIR/$DIRNAME"
      mv "$FILE" "$RES_DIR/$DIRNAME/"
    fi
  done

done

echo "[All Completed] Results saved to $RES_DIR, $LOG_FILE, and logs/*/summary.csv"
SCRIPT_END=$(date +%s)
ELAPSED=$((SCRIPT_END - SCRIPT_START))
echo "[Total Time] $ELAPSED seconds elapsed."
# less "$LOG_FILE"
