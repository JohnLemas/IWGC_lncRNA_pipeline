#!/bin/bash
# Phase 4 - intersect model calls, pull sequences/GTF, classify, join TPM + tissue tracking, write the final table.
#
#   16 intersect + classify + finalize
#      needs: all five model outputs (10-14), the TPM summary (09), the tracking summary (15), the TE annotation
#
# Usage: bash submit_4_final_outputs.sh [--dry-run]
#   --dry-run  show what would be submitted without submitting
#
# No job-ID dependencies: it submits only when every input exists and is non-empty, and says what is missing otherwise.

set -eo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"

need_file "$TE_GFF3"    "EDTA TE annotation (config: TE_GFF3)"
need_file "$REF_ANN"    "mRNA GTF (config: REF_ANN)"
need_file "$RBIO_SIF"   "R tools image"
need_file "$FEELNC_SIF" "FEELnc image (classifier)"
need_dir  "$FEELNC_HOME" "FEELnc install (config: FEELNC_HOME)"
preflight_done

step_single final 16_intersect_classify.sb final --done "$FINAL_ANNOTATION" \
    --need "$FEELNC_GTF" "$CPAT_OUT" "$LNCBOOST_PREDS" "$LNCFINDER_CSV" "$DIAMOND_OUT" \
           "$TPM_SUMMARY" "$TRACKING_SUMMARY" "$CANDIDATES_GTF" "$CANDIDATES_FASTA" "$TE_GFF3" "$REF_ANN"

print_status
