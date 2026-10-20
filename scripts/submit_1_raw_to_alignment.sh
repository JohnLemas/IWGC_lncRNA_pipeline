#!/bin/bash
# Phase 1 - raw reads through alignment.
#
#   00 STAR index ────────────────────────────────────────────┐
#   01 fastp (array) ─┐                                       ├─> 05 STAR align (array)
#   02 barrnap/index ─┴─> 03 SortMeRNA (array) ───────────────┘
#                     │                    └─> 04 QC clean
#                     └─> 04 QC trimmed          04 QC raw runs immediately
#
# Usage: ./submit_1_raw_to_alignment.sh [--dry-run] [--standalone]
#   --dry-run     run the checks and print the sbatch commands without submitting
#   --standalone  accepted for consistency; phase 1 has no upstream phase
#
# Dependencies: 03 = aftercorr/afterok(01) + afterok(02);  05 = aftercorr/afterok(03) + afterok(00)
#               04 stages use afterany so QC reports still appear if a step partly fails.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"

# ---- preflight ----
count_samples
need_file "$GENOME"          "genome fasta (config: GENOME)"
need_file "$BARRNAP_GENOME"  "genome for barrnap (config: BARRNAP_GENOME)"
need_file "$REF_ANN"         "mRNA GTF (config: REF_ANN)"
need_dir  "$RAW_FASTQ"       "raw fastq directory (config: RAW_FASTQ)"
need_file "$SMR_SIF"         "sortmerna/barrnap image (config: SMR_SIF)"
need_file "$SMR_DEFAULT_DB"  "SortMeRNA default db (config: SMR_DEFAULT_DB)"
check_raw_fastq
preflight_done

# ---- submit ----
STAR_IDX_JOB=$(submit 00_star_genomegenerate.sb star_index)
FASTP_JOB=$(submit    01_fastp.sb               fastp    --array="1-${N_SAMPLES}")
BARRNAP_JOB=$(submit  02_barrnap.sb             barrnap)

DEP=$(dep "${ARRAY_DEP_TYPE}:${FASTP_JOB}" "afterok:${BARRNAP_JOB}")
SMR_JOB=$(submit      03_sortmerna.sb           sortmerna --array="1-${N_SAMPLES}" ${DEP:+"$DEP"})

# QC: one job per stage so each waits only on what it needs
QC_RAW_JOB=$(submit   04_QC.sb qc --array=0)
DEP=$(dep "afterany:${FASTP_JOB}")
QC_TRIM_JOB=$(submit  04_QC.sb qc --array=1 ${DEP:+"$DEP"})
DEP=$(dep "afterany:${SMR_JOB}")
QC_CLEAN_JOB=$(submit 04_QC.sb qc --array=2 ${DEP:+"$DEP"})

DEP=$(dep "${ARRAY_DEP_TYPE}:${SMR_JOB}" "afterok:${STAR_IDX_JOB}")
ALIGN_JOB=$(submit    05_alignment_array.sb     align    --array="1-${N_SAMPLES}" ${DEP:+"$DEP"})

save_ids phase1 STAR_IDX_JOB FASTP_JOB BARRNAP_JOB SMR_JOB QC_RAW_JOB QC_TRIM_JOB QC_CLEAN_JOB ALIGN_JOB
info "Phase 1 submitted. BAMs will appear in ${STAR_DIR}"
monitor_hint
