#!/bin/bash
# Submit the whole pipeline as one dependency chain (phases 1 -> 4).
#
# Usage: ./submit_all.sh [--dry-run]
#   --dry-run  run every phase's checks and print all sbatch commands without submitting
#
# Phases can also be submitted one at a time (submit_1_... through submit_4_...);
# each reads the job ids the previous phase saved under LOG_DIR/job_ids/.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for a in "$@"; do [[ "$a" == "--dry-run" || "$a" == "-h" || "$a" == "--help" ]] || { echo "ERROR: only --dry-run is accepted here (got $a)" >&2; exit 1; }; done
source "${HERE}/submit_lib.sh"
parse_flags "$@"

# start from a clean slate so stale job ids from an earlier run can't be picked up
rm -f "${JOBID_DIR}"/phase*.env

for phase in submit_1_raw_to_alignment.sh submit_2_stringtie.sh submit_3_model_predictions.sh submit_4_final_outputs.sh; do
    echo "================ ${phase} ================" >&2
    "${HERE}/${phase}" "$@"
done
echo "All phases submitted." >&2
