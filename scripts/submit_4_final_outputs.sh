#!/bin/bash
# Phase 4 - intersect model calls, pull sequences/GTF, classify, join TPM + tissue tracking, write final table.
#
#   10-14 model predictions (phase 3) ─┐
#   09 TPM matrix (phase 2) ───────────┼─> 16 intersect + classify + finalize
#   15 tissue tracking (phase 2) ──────┘
#
# Usage: ./submit_4_final_outputs.sh [--dry-run] [--standalone]
#   --dry-run     run the checks and print the sbatch commands without submitting
#   --standalone  do not wait on phases 2/3 (use when all their outputs already exist)
#
# Dependencies: 16 = afterok(10,11,12,13,14) + afterok(09) + afterok(15)

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"
load_upstream phase2 phase3

# ---- preflight ----
need_file "$TE_GFF3"  "EDTA TE annotation (config: TE_GFF3)"
need_file "$REF_ANN"  "mRNA GTF (config: REF_ANN)"
need_file "$RBIO_SIF" "R tools image"
need_file "$FEELNC_SIF" "FEELnc image (classifier)"
need_dir  "$FEELNC_HOME" "FEELnc install (config: FEELNC_HOME)"
if (( STANDALONE )); then
    need_file "$FEELNC_GTF"        "step 10 output"
    need_file "$CPAT_OUT"          "step 11 output"
    need_file "$LNCBOOST_PREDS"    "step 12 output"
    need_file "$LNCFINDER_CSV"     "step 13 output"
    need_file "$DIAMOND_OUT"       "step 14 output"
    need_file "$TPM_SUMMARY"       "step 09 output (needs a tissue column in METADATA)"
    need_file "$TRACKING_SUMMARY"  "step 15 output"
fi
preflight_done

# ---- submit ----
DEP=$(dep "afterok:${FEELNC_JOB:-}" "afterok:${CPAT_JOB:-}" "afterok:${LNCBOOST_JOB:-}" \
          "afterok:${LNCFINDER_JOB:-}" "afterok:${DIAMOND_JOB:-}" \
          "afterok:${TPM_JOB:-}" "afterok:${TRACK_JOB:-}")
FINAL_JOB=$(submit 16_intersect_classify.sb final ${DEP:+"$DEP"})

save_ids phase4 FINAL_JOB
info "Phase 4 submitted. Final table will be: ${FINAL_ANNOTATION}"
monitor_hint
