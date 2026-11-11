#!/bin/bash
# Phase 4 - intersect model calls, then build the final deliverables (sequences, GTFs, classification, TE overlap, annotation tables).
#
#   16 intersect  (intersect)  needs: all five model outputs (10-14)
#   17 finalize   (finalize)   needs: the passing IDs from 16, the candidates GTF/FASTA (07), the TPM summary (09),
#                              the tracking summary (15), the TE annotation, the mRNA GTF
#
# Usage: bash submit_4_final_outputs.sh [--dry-run]
#   --dry-run  show what would be submitted without submitting
#
# No job-ID dependencies: each run submits only the steps whose inputs are ready and whose outputs
# are not done, then prints a status table. Typical order: run 1 -> 16;  run 2 -> 17.

set -eo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"

need_file "$TE_GFF3"    "EDTA TE annotation (config: TE_GFF3)"
need_file "$REF_ANN"    "mRNA GTF (config: REF_ANN)"
need_file "$RBIO_SIF"   "R tools image"
need_file "$FEELNC_SIF" "FEELnc image (classifier)"
need_dir  "$FEELNC_HOME" "FEELnc install (config: FEELNC_HOME)"
preflight_done

step_single intersect 16_intersect_classify.sb intersect --done "$FINAL_IDS" \
    --need "$FEELNC_GTF" "$CPAT_OUT" "$LNCBOOST_PREDS" "$LNCFINDER_CSV" "$DIAMOND_OUT"

step_single finalize 17_finalize.sb finalize \
    --done "$FINAL_ALL_FASTA" "$FINAL_ALL_GTF" "$FINAL_ALL_ANNOTATION" \
           "$FINAL_COL_FASTA" "$FINAL_COL_GTF" "$FINAL_COL_ANNOTATION" "$FINAL_ISOFORM_TABLE" \
    --need "$FINAL_IDS" "$CANDIDATES_GTF" "$CANDIDATES_FASTA" "$TPM_SUMMARY" "$TRACKING_SUMMARY" "$TE_GFF3" "$REF_ANN"

print_status