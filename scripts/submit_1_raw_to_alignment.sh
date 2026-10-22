#!/bin/bash
# Phase 1 - raw reads through alignment.
#
#   00 STAR index                     needs: genome + mRNA GTF
#   01 fastp (per sample)             needs: raw fastq
#   02 barrnap + SortMeRNA index      needs: BARRNAP_GENOME
#   03 SortMeRNA (per sample)         needs: fastp output + the barrnap index
#   05 STAR align (per sample)        needs: SortMeRNA output + the STAR index
#   04 QC (raw / trimmed / clean)     runs when its reads exist
#
# Usage: bash submit_1_raw_to_alignment.sh [--dry-run]
#   --dry-run  show what would be submitted without submitting
#
# No job-ID dependencies: each run submits only the steps whose inputs are ready and whose outputs
# are not done, then prints a status table. Run it, wait for the jobs, check the logs, run it again.
# Typical order: run 1 -> 00, 01, 02 and raw QC;  run 2 -> 03;  run 3 -> 05.

set -eo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"
load_samples

need_file "$GENOME"          "genome fasta (config: GENOME)"
need_file "$REF_ANN"         "mRNA GTF (config: REF_ANN)"
check_gtf_ids "$REF_ANN"
need_dir  "$RAW_FASTQ"       "raw fastq directory (config: RAW_FASTQ)"
need_file "$BARRNAP_GENOME"  "genome for barrnap (config: BARRNAP_GENOME)"
need_file "$SMR_SIF"         "sortmerna/barrnap image (config: SMR_SIF)"
need_file "$SMR_DEFAULT_DB"  "SortMeRNA default db (config: SMR_DEFAULT_DB)"
preflight_done

# ---- per-sample tests: i = 0-based index into OLD_IDS / NEW_IDS ----
raw_ok()     { local o=${OLD_IDS[$1]}; [[ -s "${RAW_FASTQ}/${o}_1.fastq.gz" || -s "${RAW_FASTQ}/${o}.fastq.gz" ]]; }
fastp_done() { local n=${NEW_IDS[$1]}; [[ -s "${FASTP_DIR}/${n}_fastp.json" ]] &&
               { [[ -s "${FASTP_DIR}/${n}_1_trim.fastq.gz" && -s "${FASTP_DIR}/${n}_2_trim.fastq.gz" ]] || [[ -s "${FASTP_DIR}/${n}_trim.fastq.gz" ]]; }; }
smr_done()   { local n=${NEW_IDS[$1]}
               { [[ -s "${SMR_CLEAN_DIR}/${n}_1_sortmerna.fastq.gz" && -s "${SMR_CLEAN_DIR}/${n}_2_sortmerna.fastq.gz" ]] || [[ -s "${SMR_CLEAN_DIR}/${n}_singles_sortmerna.fastq.gz" ]]; }; }
bam_done()   { local n=${NEW_IDS[$1]}; [[ -s "${STAR_DIR}/${n}_Aligned.sortedByCoord.out.bam" && -s "${STAR_DIR}/${n}_Aligned.sortedByCoord.out.bam.bai" ]]; }

# ---- steps ----
step_single star_index 00_star_genomegenerate.sb star_index --done "$INTRON_LIMITS" "${STAR_INDEX}/SAindex" --need "$GENOME" "$REF_ANN"
step_array  fastp      01_fastp.sb               fastp      fastp_done raw_ok
step_single barrnap    02_barrnap.sb             barrnap    --done "$SMR_CUSTOM_DB" "$SMR_IDX_DIR" --need "$BARRNAP_GENOME" "$SMR_DEFAULT_DB" "$SMR_SIF"
step_array  sortmerna  03_sortmerna.sb           sortmerna  smr_done fastp_done "$SMR_CUSTOM_DB" "$SMR_IDX_DIR"
step_array  align      05_alignment_array.sb     align      bam_done smr_done "${STAR_INDEX}/SAindex" "$INTRON_LIMITS"

step_qc 0 qc_raw   "${QC_RAW_DIR}/raw_report.html"      "$RAW_FASTQ"     ""
step_qc 1 qc_trim  "${QC_TRIM_DIR}/trimmed_report.html" "$FASTP_DIR"     fastp
step_qc 2 qc_clean "${QC_CLEAN_DIR}/clean_reads.html"   "$SMR_CLEAN_DIR" sortmerna

print_status
