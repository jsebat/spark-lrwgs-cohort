#!/usr/bin/env bash
# Read-level review of candidate de novo SVs (required before a candidate is called validation-ready).
# A joint-caller genotype on a handful of reads is not enough: a constitutional heterozygous deletion must
# (1) have junction reads at BOTH breakpoints, (2) deplete one haplotype's reads across the interval, and
# (3) leave no heterozygous SNVs inside the interval on the deleted haplotype (het fraction drops towards 0;
#     in a female X or autosome, het SNVs persisting at the flanking rate on both haplotypes exclude a
#     constitutional deletion and leave only mosaicism or chimeric reads).
# Usage: sv_read_review.sh <candidates.tsv> <manifest.tsv> <small_variant_vcf_pattern> <out.tsv>
#   candidates.tsv : family  proband  chrom  start  end  svtype  ... (header; from denovo_sv_priority)
#   manifest.tsv   : sample_id family_id father_id mother_id sex affected role bam_path (header)
#   vcf pattern    : path with {FAMILY} placeholder to the family's phased small-variant VCF (+ index alongside)
# Requires: samtools, bcftools (via the container in $SIF), haplotagged BAMs with HP tags and .bai alongside.
set -uo pipefail
CAND="${1:?candidates.tsv}"; MAN="${2:?manifest.tsv}"; VCFPAT="${3:?vcf pattern with {FAMILY}}"; OUT="${4:?out.tsv}"
SIF="${SIF:?set SIF to the container with samtools/bcftools}"
sx(){ singularity exec -B /expanse:/expanse -B "$HOME:$HOME" "$SIF" "$@"; }
FLANK=${FLANK:-20000}; BP=${BP:-300}; MINBIG=${MINBIG:-200}
hetfrac(){ # vcf sample region -> "sites het_frac"
  sx bcftools view -r "$3" -s "$2" "$1" 2>/dev/null | sx bcftools query -f '[%GT]\n' | awk '{n++; if($1~/^(0\/1|1\/0|0\|1|1\|0)$/)h++} END{printf "%d %.2f", n, (n?h/n:0)}'
}
hpcount(){ # bam region -> "total HP1 HP2"
  echo "$(sx samtools view -c "$1" "$2") $(sx samtools view -c -d HP:1 "$1" "$2") $(sx samtools view -c -d HP:2 "$1" "$2")"
}
junction(){ # bam region -> "total signature HP1 HP2 noHP"
  sx samtools view "$1" "$2" | awk -v M=$MINBIG '{sig=0; if($0 ~ /SA:Z/) sig=1; c=$6; while(match(c,/[0-9]+D/)){ if(substr(c,RSTART,RLENGTH-1)+0>=M) sig=1; c=substr(c,RSTART+RLENGTH)}; n++; if(sig){d++; hp="0"; if(match($0,/HP:i:[0-9]/)) hp=substr($0,RSTART+5,1); t[hp]++}} END{printf "%d %d %d %d %d", n, d+0, t["1"]+0, t["2"]+0, t["0"]+0}'
}
printf 'family\tproband\tchrom\tstart\tend\tsvtype\tjunction_left(total,sig,HP1,HP2,noHP)\tjunction_right\treads_inside(total,HP1,HP2)\treads_flankL\treads_flankR\thet_inside(sites,frac)\thet_flankL\thet_flankR\treview\n' > "$OUT"
tail -n +2 "$CAND" | while IFS=$'\t' read -r fam pro chrom start end svtype rest; do
  bam=$(awk -F'\t' -v s="$pro" 'NR>1 && $1==s{print $8}' "$MAN")
  vcf="${VCFPAT//\{FAMILY\}/$fam}"
  [ -s "$bam" ] || { echo -e "$fam\t$pro\t$chrom\t$start\t$end\t$svtype\tBAM_MISSING" >> "$OUT"; continue; }
  jl=$(junction "$bam" "$chrom:$((start-BP))-$((start+BP))"); jr=$(junction "$bam" "$chrom:$((end-BP))-$((end+BP))")
  ri=$(hpcount "$bam" "$chrom:$start-$end"); rl=$(hpcount "$bam" "$chrom:$((start-FLANK))-$((start-BP))"); rr=$(hpcount "$bam" "$chrom:$((end+BP))-$((end+FLANK))")
  hi="NA NA"; hl="NA NA"; hr="NA NA"
  if [ -s "$vcf" ] && [ $((end-start)) -ge 5000 ]; then hi=$(hetfrac "$vcf" "$pro" "$chrom:$start-$end"); hl=$(hetfrac "$vcf" "$pro" "$chrom:$((start-FLANK))-$((start-BP))"); hr=$(hetfrac "$vcf" "$pro" "$chrom:$((end+BP))-$((end+FLANK))"); fi
  # verdict: constitutional-consistent if junction reads at both ends AND (interval < 5 kb OR one haplotype depleted >= 40% inside OR het fraction inside <= half of flank mean)
  review=$(awk -v jl="$jl" -v jr="$jr" -v ri="$ri" -v rl="$rl" -v rr="$rr" -v hi="$hi" -v hl="$hl" -v hr="$hr" -v len=$((end-start)) 'BEGIN{
    split(jl,a," "); split(jr,b," "); split(ri,c," "); split(rl,d," "); split(rr,e," "); split(hi,f," "); split(hl,g," "); split(hr,h," ");
    if(a[2]<2 || b[2]<2){print "FAIL_junction_reads(<2 at a breakpoint)"; exit}
    if(len<5000){print "PASS_small_event(junction reads both ends)"; exit}
    fl1=(d[2]+e[2])/2; fl2=(d[3]+e[3])/2; dep=0; if(fl1>0 && c[2]<=0.6*fl1) dep=1; if(fl2>0 && c[3]<=0.6*fl2) dep=1;
    hf=0; if(f[1]!="NA"){fm=(g[2]+h[2])/2; if(f[1]>=20 && f[2]<=0.5*fm) hf=1; if(f[1]>=20 && fm>0.05 && f[2]>0.8*fm) hf=-1}
    if(dep && hf>=0){print "PASS_constitutional_consistent"; exit}
    if(hf<0){print "FAIL_het_SNVs_persist_inside(both haplotypes present: mosaic or artifact)"; exit}
    if(!dep){print "FAIL_no_haplotype_depletion"; exit}
    print "REVIEW_manually"}')
  echo -e "$fam\t$pro\t$chrom\t$start\t$end\t$svtype\t$jl\t$jr\t$ri\t$rl\t$rr\t$hi\t$hl\t$hr\t$review" | tr ' ' ',' >> "$OUT"
done
echo "wrote $OUT"; column -t -s$'\t' "$OUT" | cut -c1-220
