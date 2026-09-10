#!/usr/bin/env bash
# Layer 2 annotation resources (hg38) -> resources/layer2/. Run on Expanse from the repo root.
# All public. Records URLs + date in resources/layer2/SOURCE.md. Re-running skips files that already exist.
set -euo pipefail
R="${1:-resources/layer2}"; mkdir -p "$R"; cd "$R"
declare -A URL=(
  [GRCh38-cCREs.bed]=https://downloads.wenglab.org/V3/GRCh38-cCREs.bed                                  # ENCODE SCREEN v3 cCREs (PLS/pELS/dELS/CA-*)
  [cpgIslandExt.hg38.txt.gz]=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz   # UCSC CpG islands
  [rmsk.hg38.txt.gz]=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz                   # RepeatMasker (satellites, LINE-1)
  [ncbiRefSeqCurated.hg38.txt.gz]=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz  # gene bodies
  [hg38.chrom.sizes]=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
  [solo_WCGW_inCommonPMDs_hg38.bed.gz]=https://zhouserver.research.chop.edu/GenomeAnnotation/hg38/solo_WCGW_inCommonPMDs_hg38.bed.gz  # Zhou 2018 Nat Genet
  [PMD_coordinates_hg38.bed.gz]=https://zhouserver.research.chop.edu/GenomeAnnotation/hg38/PMD_coordinates_hg38.bed.gz              # Zhou 2018 PMD/HMD
)
{
  echo "# Layer 2 resources (hg38)"; echo; echo "Fetched $(date -u +%F) by tools/fetch_layer2_resources.sh."; echo
  echo "| file | URL | bytes | sha256 |"; echo "|---|---|---|---|"
} > SOURCE.md.tmp
for f in "${!URL[@]}"; do
  [ -s "$f" ] || curl -sfL -o "$f" "${URL[$f]}"
  echo "| $f | ${URL[$f]} | $(stat -c %s "$f") | $(sha256sum "$f" | cut -c1-16)... |" >> SOURCE.md.tmp
done
mv SOURCE.md.tmp SOURCE.md
echo "done: $(ls | wc -l) files in $R"
