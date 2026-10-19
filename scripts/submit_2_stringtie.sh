#!/bin/bash
# Phase 2 - StringTie assembly, merge, quantification, TPM matrix, tissue tracking.
#
#   05 align (phase 1) ─> 06 assemble (array) ─> 07 merge + candidate fasta ─┬─> 08 quantify (array) ─> 09 TPM matrix
#                                                                              └─> 15 gffcompare tracking
#
# Usage: ./submit_2_stringtie.sh [--dry-run] [--standalone]
#   --dry-run     run the checks and print the sbatch commands without submitting
#   --standalone  do not wait on phase 1 (use when the BAMs already exist)
#
# Dependencies: 06 = aftercorr/afterok(05);  07 = afterok(06);  08, 15 = afterok(07);  09 = afterok(08)
# 15 (tracking) lives here because it only needs the sample assemblies + candidates from 06/07;
# phase 4 waits on it. 09 and 15 do not depend on the model-prediction phase.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"
load_upstream phase1

# ---- preflight ----
count_samples
need_file "$REF_ANN"  "mRNA GTF (config: REF_ANN)"
need_file "$GENOME"   "genome fasta (config: GENOME)"
need_file "$RBIO_SIF" "R tools image (config: RBIO_SIF)"
if (( STANDALONE )); then need_dir "$STAR_DIR" "STAR BAM directory (config: STAR_DIR)"; fi
preflight_done

# ---- submit ----
DEP=$(dep "${ARRAY_DEP_TYPE}:${ALIGN_JOB:-}")
ASSEMBLE_JOB=$(submit 06_stringtie_assemble.sb stringtie --array="1-${N_SAMPLES}" ${DEP:+"$DEP"})

DEP=$(dep "afterok:${ASSEMBLE_JOB}")
MERGE_JOB=$(submit    07_stringtie_merge.sb    stringtie ${DEP:+"$DEP"})

DEP=$(dep "afterok:${MERGE_JOB}")
QUANT_JOB=$(submit    08_stringtie_quantify.sb stringtie --array="1-${N_SAMPLES}" ${DEP:+"$DEP"})
TRACK_JOB=$(submit    15_track_preds.sb        track     ${DEP:+"$DEP"})

DEP=$(dep "afterok:${QUANT_JOB}")
TPM_JOB=$(submit      09_build_tpm_matrix.sb   stringtie ${DEP:+"$DEP"})

save_ids phase2 ASSEMBLE_JOB MERGE_JOB QUANT_JOB TPM_JOB TRACK_JOB
info "Phase 2 submitted. Candidates: ${CANDIDATES_GTF} / ${CANDIDATES_FASTA}"
monitor_hint
