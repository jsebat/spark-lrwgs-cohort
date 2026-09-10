#!/usr/bin/env bash
# Zhou-lab InfiniumAnnotation hg38 manifests + masks (release v8.1, 2026-07) and the UCSC hg19->hg38 chain.
# Run from the repo root on any machine that needs resources/manifests (they are gitignored):  tools/fetch_manifests.sh
# Checksums are those recorded at first download (2026-09-09); a mismatch aborts.
set -euo pipefail
R="${1:-resources}"; mkdir -p "$R/manifests" "$R/liftover"
BASE=https://github.com/zhou-lab/InfiniumAnnotationData/raw/main/Anno
declare -A SHA=(
  [HM450/HM450.hg38.manifest.tsv.gz]=b14f51ee6a80fb98d2d8485807ffcb2af1f7471d156c99c192211dc98b29817b
  [EPIC/EPIC.hg38.manifest.tsv.gz]=93fc4de126c9ad8bcfba3d16f2fddb8a2e2a6695a49c7cd1ce1de7198ac7ee56
  [EPICv2/EPICv2.hg38.manifest.tsv.gz]=85f52e02fbc13505ca62aaab7d7c13a73a9f8df13dfbdf988d335a81df3ad955
  [HM450/HM450.hg38.mask.tsv.gz]=5de40616d0f6b400204663583bdadfc2cdbe8c57b9c619195f385f258e48c416
  [EPIC/EPIC.hg38.mask.tsv.gz]=4095b0a64f6dcb1795278792244ba572280a0639248423f616121876bd372f56
  [EPICv2/EPICv2.hg38.mask.tsv.gz]=c19947ba0fb264ddf696a88d3be30a0df46858c5342c5b376521e5b00d1d587c
)
{ echo "# Infinium manifests (hg38), fetched $(date -u +%F) by tools/fetch_manifests.sh"; echo; echo "| file | URL | sha256 |"; echo "|---|---|---|"; } > "$R/manifests/SOURCE.md"
for k in "${!SHA[@]}"; do
  f="$R/manifests/$(basename "$k")"
  [ -s "$f" ] || curl -sfL -o "$f" "$BASE/$k"
  got=$(sha256sum "$f" | cut -d" " -f1)
  [ "$got" = "${SHA[$k]}" ] || { echo "CHECKSUM MISMATCH for $f: $got"; exit 1; }
  echo "| $(basename "$k") | $BASE/$k | $got |" >> "$R/manifests/SOURCE.md"
  echo "ok $(basename "$k")"
done
C="$R/liftover/hg19ToHg38.over.chain.gz"
[ -s "$C" ] || curl -sfL -o "$C" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
echo "| hg19ToHg38.over.chain.gz | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | $(sha256sum "$C" | cut -d" " -f1) |" > "$R/liftover/SOURCE.md"
echo "ok hg19ToHg38.over.chain.gz"
