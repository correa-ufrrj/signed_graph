#!/bin/bash

# Directory where your datasets are
DATASET_DIR="$1"

if [ -z "$DATASET_DIR" ]; then
  echo "Usage: $0 /path/to/datasets"
  exit 1
fi

SCRIPT_DIR="$(dirname "$0")"
RUN_BATCH="$SCRIPT_DIR/run_batch.sh"

#echo "[Batch] No cuts..."
#"$RUN_BATCH" "$DATASET_DIR" --no-build
#
#echo "[Batch] Net-degree cuts only..."
#"$RUN_BATCH" "$DATASET_DIR" --netdeg --no-build
#
#echo "[Batch] All triangles cuts only..."
#"$RUN_BATCH" "$DATASET_DIR" --triangles --no-build

echo "[Batch] Negative-cycle cuts only..."
"$RUN_BATCH" "$DATASET_DIR" --negcycles # --no-build

#echo "[Batch] Negative-cycle cover cuts only..."
#"$RUN_BATCH" "$DATASET_DIR" --negcyclescover --no-build
#
#echo "[Batch] Both net-degree + triangles cuts ..."
#"$RUN_BATCH" "$DATASET_DIR" --netdeg --triangles --no-build

#echo "[Batch] Both net-degree + negative-cycle cuts..."
#"$RUN_BATCH" "$DATASET_DIR" --netdeg --negcycles --no-build

#echo "[Batch] Both net-degree + negative-cycle cover cuts..."
#"$RUN_BATCH" "$DATASET_DIR" --netdeg --negcyclescover --no-build

#echo "[Batch] Both triangles + negative-cycle cuts..."
#"$RUN_BATCH" "$DATASET_DIR" --triangles --negcycles --no-build

#echo "[Batch] Both triangles + negative-cycle cover cuts..."
#"$RUN_BATCH" "$DATASET_DIR" --triangles --negcyclescover --no-build
#
#echo "[Batch] Both net-degree + triangles cuts + negcycles ..."
#"$RUN_BATCH" "$DATASET_DIR" --netdeg --triangles --negcycles --no-build
