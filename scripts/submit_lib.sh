#!/bin/bash
# submit_lib.sh - shared helpers for the submit_*.sh scripts. Source it; don't run it.

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR
export CONFIG="${CONFIG:-${PIPELINE_DIR}/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck source=config.sh
source "$CONFIG"

STANDALONE=0
DRY_RUN=0
MISSING=0

die()  { echo -e "ERROR: $*" >&2; exit 1; }
info() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# parse_flags "$@"  ->  --standalone | --dry-run | -h
parse_flags() {
    local a
    for a in "$@"; do
        case "$a" in
            --standalone) STANDALONE=1 ;;
            --dry-run)    DRY_RUN=1 ;;
            -h|--help)    sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
            *)            die "unknown option: $a (try --help)" ;;
        esac
    done
    if (( DRY_RUN )); then JOBID_DIR="${LOG_DIR}/job_ids_dryrun"; else JOBID_DIR="${LOG_DIR}/job_ids"; fi
    mkdir -p "$JOBID_DIR"
}

# ---- preflight checks (collect every problem, then stop once) ----------------
need_file() { if [[ ! -s "$1" ]]; then echo "  MISSING file: $1   ($2)" >&2; MISSING=1; fi; }
need_dir()  { if [[ ! -d "$1" ]]; then echo "  MISSING dir:  $1   ($2)" >&2; MISSING=1; fi; }
preflight_done() {
    if (( MISSING )); then die "preflight failed - fix the paths above (edit $CONFIG)"; fi
    return 0
}

# count_samples -> sets N_SAMPLES; checks every metadata row has tab-delimited columns 1 and 2
count_samples() {
    [[ -s "$METADATA" ]] || die "metadata not found or empty: $METADATA"
    N_SAMPLES=$(awk 'END{print NR}' "$METADATA")
    local bad
    bad=$(tr -d '\r' < "$METADATA" | awk -F'\t' '$1=="" || $2=="" {print "  row " NR ": " $0}')
    [[ -z "$bad" ]] || die "metadata rows missing column 1 or 2 (file must be TAB-delimited):\n${bad}"
    info "${N_SAMPLES} samples in ${METADATA}"
}

# check_raw_fastq -> every sequencing ID (col 1) has <id>_1.fastq.gz or <id>.fastq.gz
check_raw_fastq() {
    local id
    while IFS=$'\t' read -r id _; do
        id="${id%$'\r'}"
        if [[ ! -f "${RAW_FASTQ}/${id}_1.fastq.gz" && ! -f "${RAW_FASTQ}/${id}.fastq.gz" ]]; then
            echo "  MISSING raw fastq for sequencing ID '${id}' in ${RAW_FASTQ}" >&2
            MISSING=1
        fi
    done < "$METADATA"
}

# ---- dependency helpers ------------------------------------------------------
# dep type:id [type:id ...] -> prints --dependency=type:id1:id2,type2:id3 (entries with empty ids are skipped)
dep() {
    local -A by_type=()
    local order=() p t id out=""
    for p in "$@"; do
        [[ "$p" == *:* ]] || continue
        t="${p%%:*}"; id="${p#*:}"
        [[ -n "$id" ]] || continue
        [[ -n "${by_type[$t]+x}" ]] || order+=("$t")
        by_type[$t]+=":${id}"
    done
    (( ${#order[@]} )) || return 0
    for t in "${order[@]}"; do out+="${out:+,}${t}${by_type[$t]}"; done
    echo "--dependency=${out}"
}

# submit <script.sb> <log-subdir> [sbatch args...]  -> prints the job id on stdout
# --output is always set from LOG_DIR (config.sh); --kill-on-invalid-dep cancels the job if an upstream job fails.
submit() {
    local script="$1" logsub="$2"; shift 2
    local sb="${PIPELINE_DIR}/${script}" pat="%x-%j.SLURMout" a
    [[ -f "$sb" ]] || die "step script not found: $sb"
    for a in "$@"; do [[ "$a" == --array=* ]] && pat="%x-%A_%a.SLURMout"; done
    local cmd=(sbatch --parsable --kill-on-invalid-dep=yes --export="ALL,CONFIG=${CONFIG}"
               --output="${LOG_DIR}/${logsub}/${pat}" "$@" "$sb")
    if (( DRY_RUN )); then
        info "DRY RUN: ${cmd[*]}"
        echo "DRY_${script%%_*}"
        return 0
    fi
    mkdir -p "${LOG_DIR}/${logsub}"
    local jid
    jid=$("${cmd[@]}") || die "sbatch failed for ${script}"
    jid="${jid%%;*}"
    info "submitted ${script} -> job ${jid}"
    echo "$jid"
}

# ---- job-id hand-off between phases -----------------------------------------
# save_ids <phase> VAR1 VAR2 ...   writes VAR=value lines to ${JOBID_DIR}/<phase>.env
save_ids() {
    local phase="$1" v; shift
    local f="${JOBID_DIR}/${phase}.env"
    : > "$f"
    for v in "$@"; do printf '%s=%q\n' "$v" "${!v}" >> "$f"; done
    info "job IDs saved -> $f"
}

# load_upstream <phase> [<phase> ...]   sources earlier phases' job ids (skipped with --standalone)
load_upstream() {
    local p f
    if (( STANDALONE )); then info "--standalone: not waiting on any earlier phase"; return 0; fi
    for p in "$@"; do
        f="${JOBID_DIR}/${p}.env"
        if [[ -f "$f" ]]; then
            # shellcheck disable=SC1090
            source "$f"
            info "loaded upstream job ids from ${f}"
        else
            info "NOTE: ${f} not found - ${p} will not be waited on (run it first, or pass --standalone if its outputs already exist)"
        fi
    done
    return 0
}

monitor_hint() {
    echo >&2
    echo "Monitor:  squeue -u \$USER -o '%.10i %.22j %.8T %.10M %.12r %E'" >&2
    echo "Logs:     ${LOG_DIR}/" >&2
}
