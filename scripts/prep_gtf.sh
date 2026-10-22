#!/bin/bash
# prep_gtf.sh - make REF_ANN safe for StringTie / gffcompare / FEELnc.
#
# Problem: GFF3->GTF converters (AGAT etc.) decode URL-escaped characters, so IDs like
# gene-PHT4%3B1 become gene_id "gene-PHT4;1". StringTie splits attributes on ';' and
# dies with "Error parsing attribute gene_id ('"' required for GTF)".
#
# Fix: replace ';' with '_' inside every double-quoted attribute value (columns 1-8 and
# every other attribute are left untouched). Applied consistently, so a gene keeps one
# ID across all its lines. Original is kept next to it as <name>.orig.gtf.
#
# Usage: bash prep_gtf.sh        (reads REF_ANN from config.sh; safe to re-run)

set -euo pipefail
source "${CONFIG:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh}"

[[ -s "$REF_ANN" ]] || { echo "ERROR: REF_ANN not found: $REF_ANN" >&2; exit 1; }

ORIG="${REF_ANN%.gtf}.orig.gtf"
TMP="${REF_ANN}.tmp"

# Always rebuild from the untouched original so re-running is idempotent
if [[ ! -e "$ORIG" ]]; then cp -p "$REF_ANN" "$ORIG"; echo "Backed up original -> $ORIG"; fi

# shellcheck disable=SC2016
perl -pe 's/"([^"]*)"/my $x=$1; $x =~ s{;}{_}g; "\"$x\""/ge' "$ORIG" > "$TMP"

n_in=$(wc -l < "$ORIG"); n_out=$(wc -l < "$TMP")
[[ "$n_in" -eq "$n_out" ]] || { echo "ERROR: line count changed ($n_in -> $n_out); leaving $REF_ANN untouched" >&2; rm -f "$TMP"; exit 1; }

changed=$(diff <(cat "$ORIG") "$TMP" | grep -c '^>' || true)
mv "$TMP" "$REF_ANN"
echo "Done: $n_out lines, $changed lines had ';' inside quoted values and were fixed."
echo "Cleaned GTF: $REF_ANN"
