#!/bin/bash
# Phase 3 - coding-potential / lncRNA prediction models. All five run in parallel on the candidate set.
#
#   07 merge (phase 2) ─┬─> 10 FEELnc filter
#                       ├─> 11 CPAT-plant
#                       ├─> 12 PlantLncBoost
#                       ├─> 13 LncFinder
#                       └─> 14 DIAMOND blastx vs UniProt
#
# Usage: ./submit_3_model_predictions.sh [--dry-run] [--standalone]
#   --dry-run     run the checks and print the sbatch commands without submitting
#   --standalone  do not wait on phase 2 (use when the candidate GTF/fasta already exist)
#
# Dependencies: every model = afterok(07 merge)

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"
load_upstream phase2

# ---- preflight (only things that exist before the jobs run) ----
need_file "$REF_ANN"             "mRNA GTF (config: REF_ANN)"
need_file "$FEELNC_SIF"          "FEELnc image"
need_dir  "$FEELNC_HOME"         "FEELnc install (config: FEELNC_HOME)"
need_file "$CPAT_SIF"            "CPAT image"
need_file "$CPAT_HEXAMER"        "CPAT hexamer table"
need_file "$CPAT_LOGIT"          "CPAT logit model"
need_file "$LNCBOOST_SIF"        "PlantLncBoost image"
need_dir  "$LNCBOOST_SCRIPT_DIR" "PlantLncBoost scripts"
need_file "$LNCBOOST_MODEL"      "PlantLncBoost model"
need_file "$RBIO_SIF"            "R tools image (LncFinder)"
need_file "$DIAMOND_DB"          "DIAMOND UniProt db"
if (( STANDALONE )); then
    need_file "$CANDIDATES_GTF"   "candidate GTF from step 07"
    need_file "$CANDIDATES_FASTA" "candidate fasta from step 07"
fi
preflight_done

# ---- submit ----
DEP=$(dep "afterok:${MERGE_JOB:-}")
FEELNC_JOB=$(submit    10_FEELnc_filter.sb feelnc    ${DEP:+"$DEP"})
CPAT_JOB=$(submit      11_cpat.sb          cpat      ${DEP:+"$DEP"})
LNCBOOST_JOB=$(submit  12_lncBoost.sb      lncboost  ${DEP:+"$DEP"})
LNCFINDER_JOB=$(submit 13_lncFinder.sb     lncfinder ${DEP:+"$DEP"})
DIAMOND_JOB=$(submit   14_DIAMOND.sb       diamond   ${DEP:+"$DEP"})

save_ids phase3 FEELNC_JOB CPAT_JOB LNCBOOST_JOB LNCFINDER_JOB DIAMOND_JOB
info "Phase 3 submitted (5 parallel jobs)."
monitor_hint
