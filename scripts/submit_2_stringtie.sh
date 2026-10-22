#!/bin/bash
# Phase 2 - StringTie assembly, merge, quantification, TPM matrix, tissue tracking.
#
#   06 assemble (per sample)   needs: BAMs
#   07 merge + candidates      needs: an assembly for EVERY sample
#   08 quantify (per sample)   needs: BAMs + merged GTF
#   15 tissue tracking         needs: candidates + every sample's assembly
#   09 TPM matrix              needs: every sample's t_data.ctab
#
# Usage: bash submit_2_stringtie.sh [--dry-run]
#   --dry-run  show what would be submitted without submitting
#
# No job-ID dependencies: each run submits only the steps whose inputs are ready and whose outputs
# are not done, then prints a status table. Typical order: run 1 -> 06;  run 2 -> 07;  run 3 -> 08 and 15;
# run 4 -> 09. If a sample can't be assembled, remove it from the metadata file or fix it first:
# merge waits for all of them.

set -eo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/submit_lib.sh"
parse_flags "$@"
load_samples

need_file "$REF_ANN"  "mRNA GTF (config: REF_ANN)"
check_gtf_ids "$REF_ANN"
need_file "$GENOME"   "genome fasta (config: GENOME)"
need_file "$RBIO_SIF" "R tools image (config: RBIO_SIF)"
preflight_done

bam_ok()        { local n=${NEW_IDS[$1]}; [[ -s "${STAR_DIR}/${n}_Aligned.sortedByCoord.out.bam" && -s "${STAR_DIR}/${n}_Aligned.sortedByCoord.out.bam.bai" ]]; }
assembly_done() { local n=${NEW_IDS[$1]}; [[ -s "${ASSEMBLY_DIR}/${n}_assembled.gtf" ]]; }
quant_done()    { local n=${NEW_IDS[$1]}; [[ -s "${QUANT_DIR}/${n}/${n}_quant.gtf" && -s "${QUANT_DIR}/${n}/t_data.ctab" ]]; }

ASSEMBLIES=(); CTABS=()
for n in "${NEW_IDS[@]}"; do
    ASSEMBLIES+=("${ASSEMBLY_DIR}/${n}_assembled.gtf")
    CTABS+=("${QUANT_DIR}/${n}/t_data.ctab")
done

step_array  assemble 06_stringtie_assemble.sb stringtie assembly_done bam_ok
step_single merge    07_stringtie_merge.sb    stringtie --done "$MERGED_GTF" "$CANDIDATES_GTF" "$CANDIDATES_FASTA" --need "$REF_ANN" "$GENOME" "${ASSEMBLIES[@]}"
step_array  quant    08_stringtie_quantify.sb stringtie quant_done bam_ok "$MERGED_GTF"
step_single track    15_track_preds.sb        track     --done "$TRACKING_SUMMARY" --need "$CANDIDATES_GTF" "$CANDIDATES_FASTA" "${ASSEMBLIES[@]}"
step_single tpm      09_build_tpm_matrix.sb   stringtie --done "$TPM_MATRIX" --need "$CANDIDATES_GTF" "${CTABS[@]}"

print_status
