#!/usr/bin/env bash
# Refuse to release anything that looks like PHI, a person's name, a sample/family id,
# a consulting-family reference, or a personal path. Exit 1 on any hit.
set -uo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || dirname "$(dirname "$0")")"
NAMES='gendron|quad_vip|\bvip\b|\btran\b|wgs_stephen_family|\bstephen\b|\bpeter\b|\bjoyce\b|\bvinh\b|long_read_nygc|consultants'
IDS='SP0[0-9]{6}|SF0[0-9]{6}|REACH0009[0-9]{2}|\bF0323\b|\bF0176\b'
PATHS='/home/[a-z]+|sebatlab dropbox|spark phenotype data|sparkdatarelease'
hits=0
while IFS= read -r f; do
  m=$(grep -inoE "$NAMES|$IDS|$PATHS" "$f" 2>/dev/null | head -5)
  if [ -n "$m" ]; then echo "HIT $f"; echo "$m" | sed 's/^/    /'; hits=$((hits+1)); fi
done < <( (git ls-files 2>/dev/null; git ls-files --others --exclude-standard 2>/dev/null) | grep -v '^scripts/phi_scan.sh$' )
if [ "$hits" -gt 0 ]; then echo "phi_scan: $hits file(s) with hits - DO NOT COMMIT"; exit 1; fi
echo "phi_scan: clean"
