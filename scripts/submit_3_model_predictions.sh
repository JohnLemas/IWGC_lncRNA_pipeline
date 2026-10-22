#!/bin/bash
# Phase 3 - coding-potential / lncRNA prediction models. All five need only the candidate set from step 07,
# so one run submits all of them.
#
#   10 FEELnc filter   11 CPAT-plant   12 PlantLncBoost   13 LncFinder   14 DIAMOND blastx vs UniProt
#
# Usage: bash submit_3_model_predictions.sh [--dry-run]
#   --dry-run  show what would be submitted without submitting
#
# No job-ID dependencies. Re-run it after a failure: finished models are skipped, only missing ones resubmit.

set -eo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"

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
preflight_done

step_single feelnc    10_FEELnc_filter.sb feelnc    --done "$FEELNC_GTF"      --need "$CANDIDATES_GTF" "$REF_ANN"
step_single cpat      11_cpat.sb          cpat      --done "$CPAT_OUT"        --need "$CANDIDATES_FASTA"
step_single lncboost  12_lncBoost.sb      lncboost  --done "$LNCBOOST_PREDS"  --need "$CANDIDATES_FASTA"
step_single lncfinder 13_lncFinder.sb     lncfinder --done "$LNCFINDER_CSV"   --need "$CANDIDATES_FASTA"
step_single diamond   14_DIAMOND.sb       diamond   --done "$DIAMOND_OUT"     --need "$CANDIDATES_FASTA" "$DIAMOND_DB"

print_status
