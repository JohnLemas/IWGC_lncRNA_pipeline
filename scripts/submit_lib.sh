#!/bin/bash
# submit_lib.sh - helpers for the submit_*.sh scripts. Source it; don't run it.
#
# How the submit scripts work (there are NO job-ID dependencies anywhere):
#   * every step has a "done" test: its output files exist and are non-empty
#   * a step is submitted only if it is not done, its inputs exist and are non-empty,
#     and no jobs of that step are still in the queue
#   * run a script, wait for the jobs, check the logs, then run it again or run the next one.
#     Each run submits whatever is ready, skips what is done, and prints a status table.

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR
export CONFIG="${CONFIG:-${PIPELINE_DIR}/config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck source=config.sh
source "$CONFIG"

DRY_RUN=0
MISSING=0
STATUS_ROWS=()

die()  { echo -e "ERROR: $*" >&2; exit 1; }
info() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# parse_flags "$@"  ->  --dry-run | -h
parse_flags() {
    local a
    for a in "$@"; do
        case "$a" in
            --dry-run)   DRY_RUN=1 ;;
            -h|--help)   sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
            *)           die "unknown option: $a (try --help)" ;;
        esac
    done
}

# ---- static checks (things that must exist before anything runs) -------------
need_file() { if [[ ! -s "$1" ]]; then echo "  MISSING file: $1   ($2)" >&2; MISSING=1; fi; }
need_dir()  { if [[ ! -d "$1" ]]; then echo "  MISSING dir:  $1   ($2)" >&2; MISSING=1; fi; }
preflight_done() {
    if (( MISSING )); then die "preflight failed - fix the paths above (edit $CONFIG)"; fi
    return 0
}

# check_gtf_ids <gtf>: ';' inside a quoted attribute value breaks StringTie's GTF parser
check_gtf_ids() {
    [[ -s "$1" ]] || return 0
    if ! perl -ne 'while (/"([^"]*)"/g) { exit 1 if $1 =~ /;/ }' "$1"; then
        echo "  PROBLEM: ';' inside a quoted attribute value in $1 - StringTie will fail. Fix with: bash ${PIPELINE_DIR}/prep_gtf.sh" >&2
        MISSING=1
    fi
}

# ---- samples -----------------------------------------------------------------
# load_samples -> OLD_IDS[] (column 1), NEW_IDS[] (column 2), N_SAMPLES. Array task N = metadata line N.
load_samples() {
    [[ -s "$METADATA" ]] || die "metadata not found or empty: $METADATA"
    mapfile -t OLD_IDS < <(tr -d '\r' < "$METADATA" | awk -F'\t' '{print $1}')
    mapfile -t NEW_IDS < <(tr -d '\r' < "$METADATA" | awk -F'\t' '{print $2}')
    N_SAMPLES=${#NEW_IDS[@]}
    local i bad="" dups
    for ((i = 0; i < N_SAMPLES; i++)); do
        if [[ -z "${OLD_IDS[$i]}" || -z "${NEW_IDS[$i]}" ]]; then bad+="  line $((i + 1))\n"; fi
    done
    [[ -z "$bad" ]] || die "metadata lines missing column 1 or 2 (file must be TAB-delimited):\n${bad}"
    dups=$(printf '%s\n' "${NEW_IDS[@]}" | sort | uniq -d | head -5 | tr '\n' ' ')
    [[ -z "$dups" ]] || die "duplicate informative IDs in metadata column 2 (their outputs would overwrite each other): $dups"
    info "${N_SAMPLES} samples in ${METADATA}"
}

# ---- tests ---------------------------------------------------------------------
# nonempty <path>: a file with content, or a directory containing at least one non-empty file
nonempty() {
    if [[ -d "$1" ]]; then [[ -n "$(find "$1" -type f -size +0c -print -quit 2>/dev/null)" ]]
    else [[ -s "$1" ]]
    fi
}
rel() { echo "${1#"${PROJ}"/}"; }

# queued_jobs <job-name>: how many of your jobs with that name are pending/running
queued_jobs() {
    local n
    n=$(squeue -h -u "${USER:-$(id -un)}" -n "$1" -o %i 2>/dev/null | wc -l || true)
    echo "${n// /}"
}

# ranges 1 2 3 7 9 10 -> 1-3,7,9-10
ranges() {
    printf '%s\n' "$@" | sort -n | awk '
        NR==1 { s=p=$1; next }
        $1==p+1 { p=$1; next }
        { printf "%s%s", sep, (s==p ? s : s "-" p); sep=","; s=p=$1 }
        END { printf "%s%s\n", sep, (s==p ? s : s "-" p) }'
}

# jobtxt <id>: "job 123: " (blank in a dry run)
jobtxt() { if (( DRY_RUN )); then echo "dry run"; else echo "job $1"; fi; }

# ---- reporting -----------------------------------------------------------------
report() { STATUS_ROWS+=("$(printf '%-12s %-13s %s' "$1" "$2" "$3")"); }
print_status() {
    echo
    printf '%-12s %-13s %s\n' STEP STATUS DETAIL
    printf '%s\n' "${STATUS_ROWS[@]}"
    echo
    echo "DONE = outputs present and non-empty.  RUNNING = jobs still queued (wait).  WAITING = inputs not ready."
    echo "Check logs in ${LOG_DIR}/<step>/ when jobs finish, then run the next script (or this one again)."
    echo "To redo a step, delete its output files first."
    if (( DRY_RUN )); then echo "(dry run: nothing was submitted)"; fi
}

# ---- submitting ----------------------------------------------------------------
# submit <script.sb> <log-subdir> <job-name> [sbatch args...]  -> prints the job id
submit() {
    local script="$1" logsub="$2" name="$3"; shift 3
    local sb="${PIPELINE_DIR}/${script}" pat="%x-%j.SLURMout" a
    [[ -f "$sb" ]] || die "step script not found: $sb"
    for a in "$@"; do [[ "$a" == --array=* ]] && pat="%x-%A_%a.SLURMout"; done
    local cmd=(sbatch --parsable --job-name="$name" --export="ALL,CONFIG=${CONFIG}"
               --output="${LOG_DIR}/${logsub}/${pat}" "$@" "$sb")
    if (( DRY_RUN )); then
        info "DRY RUN: ${cmd[*]}"
        echo "-"
        return 0
    fi
    mkdir -p "${LOG_DIR}/${logsub}"
    local jid
    jid=$("${cmd[@]}") || die "sbatch failed for ${script}"
    echo "${jid%%;*}"
}

# step_single <label> <script> <log-subdir> [--args "<sbatch args>"] --done <files/dirs...> --need <files/dirs...>
#   done: every --done path is non-empty        need: every --need path is non-empty
step_single() {
    local label="$1" script="$2" logsub="$3"; shift 3
    local mode="" extra=() done_files=() need_files=() f missing="" n jid
    while (( $# )); do
        case "$1" in
            --args) read -r -a extra <<< "$2"; shift ;;
            --done) mode=done ;;
            --need) mode=need ;;
            *) if [[ "$mode" == done ]]; then done_files+=("$1"); else need_files+=("$1"); fi ;;
        esac
        shift
    done
    local name="${SPEC}_${label}" all_done=1
    (( ${#done_files[@]} )) || all_done=0
    for f in "${done_files[@]}"; do nonempty "$f" || { all_done=0; break; }; done
    if (( all_done )); then report "$label" DONE "$(rel "${done_files[0]}")"; return 0; fi
    n=$(queued_jobs "$name")
    if (( n > 0 )); then report "$label" RUNNING "${n} job(s) still queued/running"; return 0; fi
    local miss=()
    for f in "${need_files[@]}"; do nonempty "$f" || miss+=("$(rel "$f")"); done
    if (( ${#miss[@]} )); then
        missing=" ${miss[*]:0:3}"
        if (( ${#miss[@]} > 3 )); then missing+=" ... (+$(( ${#miss[@]} - 3 )) more)"; fi
        report "$label" WAITING "missing or empty:${missing}"; return 0
    fi
    jid=$(submit "$script" "$logsub" "$name" "${extra[@]}")
    report "$label" "$( ((DRY_RUN)) && echo WOULD-SUBMIT || echo SUBMITTED )" "$(jobtxt "$jid")"
}

# step_array <label> <script> <log-subdir> <done_fn> <need_fn> [shared inputs...]
#   done_fn i / need_fn i: test sample i (0-based index into NEW_IDS/OLD_IDS); return 0 for yes.
#   Only samples that are not done AND have their inputs are put in the array, so re-running
#   a script only does the missing work.
step_array() {
    local label="$1" script="$2" logsub="$3" done_fn="$4" need_fn="$5"; shift 5
    local name="${SPEC}_${label}" i f missing="" n jid
    local todo=() ready=() held=()
    for ((i = 0; i < N_SAMPLES; i++)); do
        "$done_fn" "$i" || todo+=("$i")
    done
    if (( ${#todo[@]} == 0 )); then report "$label" DONE "all ${N_SAMPLES} samples"; return 0; fi
    n=$(queued_jobs "$name")
    if (( n > 0 )); then report "$label" RUNNING "${n} task(s) still queued/running; ${#todo[@]} sample(s) not done yet"; return 0; fi
    for f in "$@"; do nonempty "$f" || missing+=" $(rel "$f")"; done
    if [[ -n "$missing" ]]; then report "$label" WAITING "missing or empty:${missing}"; return 0; fi
    for i in "${todo[@]}"; do
        if "$need_fn" "$i"; then ready+=("$((i + 1))"); else held+=("${NEW_IDS[$i]}"); fi
    done
    if (( ${#ready[@]} == 0 )); then
        report "$label" WAITING "inputs missing for all ${#todo[@]} unfinished samples (e.g. ${held[*]:0:3})"
        return 0
    fi
    local spec; spec=$(ranges "${ready[@]}")
    jid=$(submit "$script" "$logsub" "$name" --array="$spec")
    report "$label" "$( ((DRY_RUN)) && echo WOULD-SUBMIT || echo SUBMITTED )" "$(jobtxt "$jid"): ${#ready[@]} task(s) [${spec}]"
    if (( ${#held[@]} )); then
        report "" "  held back" "${#held[@]} sample(s) have no input yet: ${held[*]:0:6}$( ((${#held[@]} > 6)) && echo ' ...')"
    fi
}

# step_qc <array-index> <label> <report.html> <reads-dir> <upstream-label|"">
#   QC is "done" when its report exists and no read file is newer than it (QC itself resumes, so a
#   re-run after more samples finish only does the new files).
step_qc() {
    local idx="$1" label="$2" report_file="$3" indir="$4" upstream="$5"
    local name="${SPEC}_${label}" jid
    if ! compgen -G "${indir}/*.gz" > /dev/null; then report "$label" WAITING "no reads in $(rel "$indir")"; return 0; fi
    if [[ -n "$upstream" ]] && (( $(queued_jobs "${SPEC}_${upstream}") > 0 )); then
        report "$label" WAITING "${upstream} still running"; return 0
    fi
    if (( $(queued_jobs "$name") > 0 )); then report "$label" RUNNING "still queued/running"; return 0; fi
    if [[ -s "$report_file" && -z "$(find "$indir" -maxdepth 1 -name '*.gz' -newer "$report_file" -print -quit)" ]]; then
        report "$label" DONE "$(rel "$report_file")"; return 0
    fi
    jid=$(submit 04_QC.sb qc "$name" --array="$idx")
    report "$label" "$( ((DRY_RUN)) && echo WOULD-SUBMIT || echo SUBMITTED )" "$(jobtxt "$jid"): reads in $(rel "$indir")"
}
