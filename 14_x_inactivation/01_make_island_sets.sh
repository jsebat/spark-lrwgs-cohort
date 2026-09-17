#!/usr/bin/env bash
# Build the two CpG island sets the skew estimator needs.
#
#   chrX_cgi.bed        candidate X-linked islands; 04_xci_skew.py narrows these to the informative subset
#                       empirically, so no escape annotation is needed here.
#   autosomal_cgi.bed   islands used as the noise floor for the folded estimator. Autosomes are diploid and
#                       not subject to inactivation, so their haplotype asymmetry in the same reads is assay
#                       noise with no true signal.
#
# Input is any CpG island BED (for example the UCSC cpgIslandExt track for your build). Islands are filtered
# on length so that very short ones, where a handful of CpGs drive the mean, do not enter either set.
set -euo pipefail

CGI=""; OUT="islands"; MIN_LEN=200; MAX_LEN=5000; AUTO_N=2000; SEED=1

usage() {
  cat <<'EOF'
usage: 01_make_island_sets.sh --cgi <cpg_islands.bed> --out <dir> [options]

  --cgi       CpG island BED for the reference build (required)
  --out       output directory (default: islands)
  --min-len   minimum island length in bp (default: 200)
  --max-len   maximum island length in bp (default: 5000)
  --auto-n    number of autosomal islands to sample for the noise floor (default: 2000; 0 = keep all)
  --seed      sampling seed (default: 1)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cgi) CGI="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --min-len) MIN_LEN="$2"; shift 2;;
    --max-len) MAX_LEN="$2"; shift 2;;
    --auto-n) AUTO_N="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage; exit 1;;
  esac
done

[[ -n "$CGI" ]] || { echo "--cgi is required" >&2; usage; exit 1; }
[[ -f "$CGI" ]] || { echo "no such file: $CGI" >&2; exit 1; }
mkdir -p "$OUT"

# chrX islands. Keep the whole chromosome: PAR islands are harmless here because the informative-island
# filter in 04 removes anything that is not intermediate in females and low in males.
awk -v OFS='\t' -v lo="$MIN_LEN" -v hi="$MAX_LEN" \
  '!/^(#|track|browser)/ && $1=="chrX" && ($3-$2)>=lo && ($3-$2)<=hi {print $1,$2,$3}' \
  "$CGI" | sort -k1,1 -k2,2n > "$OUT/chrX_cgi.bed"

# Autosomal islands, restricted to the primary assembly so that alt and random contigs do not contribute
# spurious haplotype asymmetry.
awk -v OFS='\t' -v lo="$MIN_LEN" -v hi="$MAX_LEN" \
  '!/^(#|track|browser)/ && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ && ($3-$2)>=lo && ($3-$2)<=hi {print $1,$2,$3}' \
  "$CGI" | sort -k1,1 -k2,2n > "$OUT/autosomal_cgi.all.bed"

if [[ "$AUTO_N" -gt 0 ]]; then
  awk -v seed="$SEED" 'BEGIN{srand(seed)} {print rand()"\t"$0}' "$OUT/autosomal_cgi.all.bed" \
    | sort -k1,1g | head -n "$AUTO_N" | cut -f2- \
    | sort -k1,1 -k2,2n > "$OUT/autosomal_cgi.bed"
else
  cp "$OUT/autosomal_cgi.all.bed" "$OUT/autosomal_cgi.bed"
fi

printf 'chrX islands:      %s\n' "$(wc -l < "$OUT/chrX_cgi.bed")"
printf 'autosomal islands: %s (of %s)\n' \
  "$(wc -l < "$OUT/autosomal_cgi.bed")" "$(wc -l < "$OUT/autosomal_cgi.all.bed")"
printf 'written to %s/\n' "$OUT"
