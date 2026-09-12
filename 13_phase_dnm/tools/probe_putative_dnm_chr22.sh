#!/usr/bin/env bash
# probe_c.sh — putative-DNM count on chr22 for the first callset family.
# Light (one chromosome of a 3-4 sample VCF); login-node safe. Prints NO sample identifiers.
# Purpose: DESIGN P13 — how many child-het / parents-hom-ref sites exist per child, by min GQ,
# to estimate the true-DNM contamination fraction of SynthDNM's negative class on HiFi.
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
BC="${BCFTOOLS:-bcftools}"   # set BCFTOOLS in config/phase_dnm.env (e.g. a micromamba env binary)
L="${LRC:?set LRC in config/phase_dnm.env}"
M=$L/MANIFEST_lrWGS_105.tsv
FAM=$(ls "$L/callsets" | head -1)
V=$(find "$L/callsets/$FAM/out/joint_small_variants_vcf" -maxdepth 2 -name '*.vcf.gz' | head -1)
I=$(find "$L/callsets/$FAM/out/joint_small_variants_vcf_index" -maxdepth 2 -name '*.tbi' | head -1)
[ -s "$V" ] && [ -s "$I" ] || { echo "vcf or index not found"; exit 1; }
echo "roles (role affected) in first family:"
awk -F'\t' -v f="$FAM" 'NR>1 && $2==f {print $7, $6}' "$M" | sort | uniq -c
echo "vcf samples: $($BC query -l "$V##idx##$I" | wc -l)   chr22 records: $($BC view -H -r chr22 "$V##idx##$I" | wc -l)"
awk -F'\t' -v f="$FAM" 'NR>1 && $2==f && $3!="0" && $3!="" && $4!="0" && $4!="" {print $1"\t"$3"\t"$4}' "$M" |
while IFS=$'\t' read -r c fa mo; do
  echo "== one child, chr22 (child het, both parents hom-ref) =="
  $BC view -r chr22 -s "$c,$fa,$mo" "$V##idx##$I" 2>/dev/null |
  $BC query -f '%POS\t%REF\t%ALT[\t%GT:%GQ:%DP]\n' |
  awk -F'\t' '{ n++; split($4,c,":"); split($5,f,":"); split($6,m,":");
     het = (c[1]=="0/1"||c[1]=="1/0"||c[1]=="0|1"||c[1]=="1|0");
     if (het && f[1] ~ /^0[\/|]0$/ && m[1] ~ /^0[\/|]0$/) {
        d++; g=c[2]+0; if (f[2]+0<g) g=f[2]+0; if (m[2]+0<g) g=m[2]+0;
        if (g>=20) d20++; if (g>=30) d30++;
        if (g>=20 && c[3]+0>=10 && f[3]+0>=10 && m[3]+0>=10) dq++;
        if (length($2)!=length($3)) di++ } }
     END { printf "sites_chr22=%d  putative_DNM=%d  minGQ>=20:%d  minGQ>=30:%d  GQ20&DP10_all3:%d  indel_frac=%.2f\n",
                  n, d+0, d20+0, d30+0, dq+0, (d?di/d:0) }'
done
echo "scale note: chr22 is ~1.3% of autosomal variant sites; multiply by ~75 for a genome-wide order of magnitude"
