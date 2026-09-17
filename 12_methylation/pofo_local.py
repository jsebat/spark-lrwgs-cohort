#!/usr/bin/env python3
"""Parent-of-origin at ICRs using the NEAREST informative SNVs (window +/- 20 kb, then +/- 50 kb), with a block-level
consistency flag. Fixes the phase-switch artifact seen in D74 stage 1 (two 'discordant' children)."""
import os, csv, subprocess, collections
ME=os.environ["DATA_ROOT"] + "/meth/"; MAN=os.environ["MANIFEST"]
SIF=os.environ["SIF"]; sx=["singularity","exec","-B","/expanse:/expanse","-B",os.path.expanduser("~")+":"+os.path.expanduser("~"),SIF]
man={r["sample_id"]:r for r in csv.DictReader(open(MAN),delimiter="\t")}
rows=list(csv.DictReader(open(ME+"icr_parent_of_origin_check.tsv"),delimiter="\t"))
icr={}
for ln in open(ME+"icr_regions.bed"):
    c,a,b,n=ln.rstrip("\n").split("\t"); icr[n]=(c,int(a),int(b))
def run(cmd): return subprocess.run(cmd,stdout=subprocess.PIPE,universal_newlines=True).stdout
def informative(child,c,a,b,win):
    fam=man[child]["family_id"]; vcf=os.path.join(os.environ["PHASED_VCF_DIR"], "%s.phased.vcf.gz" % fam)
    txt=run(sx+["bcftools","query","-s",",".join([child,man[child]["father_id"],man[child]["mother_id"]]),"-r","%s:%d-%d"%(c,max(1,a-win),b+win),"-i",'TYPE="snp"',"-f","%POS[\t%GT:%PS]\n",vcf])
    calls=[]
    for ln in txt.splitlines():
        f=ln.split("\t")
        if len(f)<4: continue
        pos=int(f[0]); cg=f[1].split(":")[0]; fg=f[2].split(":")[0].replace("|","/"); mg=f[3].split(":")[0].replace("|","/"); ps=f[1].split(":")[1] if ":" in f[1] else "."
        if "|" not in cg or cg[0]==cg[2] or "." in fg or "." in mg: continue
        h1,h2=cg[0],cg[2]; fa=set(fg.split("/")); mo=set(mg.split("/"))
        if h1 in mo and h1 not in fa and h2 in fa: calls.append((pos,"mat",ps))
        elif h2 in mo and h2 not in fa and h1 in fa: calls.append((pos,"pat",ps))
        elif h1 in fa and h1 not in mo and h2 in mo: calls.append((pos,"pat",ps))
        elif h2 in fa and h2 not in mo and h1 in mo: calls.append((pos,"mat",ps))
    return calls
out=open(ME+"icr_parent_of_origin_local.tsv","w"); out.write("child\tfamily\ticr\texpected\tmethylated_hap\thap1_meth\thap2_meth\tn_local\tlocal_hap1_maternal_frac\tn_block\tblock_hap1_maternal_frac\tswitch_flag\tmethylated_parent\tconcordant\n")
tally=collections.Counter(); switches=0
for r in rows:
    if r["methylated_hap"]=="none": tally["no_ASM"]+=1; continue
    c,a,b=icr[r["icr"]]; ch=r["child"]
    loc=informative(ch,c,a,b,20000)
    if len(loc)<2: loc=informative(ch,c,a,b,50000)
    blk=informative(ch,c,a,b,300000)
    def frac(calls):
        return (sum(1 for _,k,_ in calls if k=="mat")/len(calls)) if calls else None
    fl=frac(loc); fb=frac(blk)
    par=""; 
    if fl is not None and len(loc)>=2 and (fl>=0.9 or fl<=0.1):
        h1mat = fl>=0.9; par=("maternal" if h1mat else "paternal") if r["methylated_hap"]=="hap1" else ("paternal" if h1mat else "maternal")
    sw = "SWITCH" if (fl is not None and fb is not None and len(loc)>=2 and len(blk)>=10 and abs(fl-fb)>0.5) else ""
    if sw: switches+=1
    conc = "" if not par else ("yes" if par==r["expected_methylated"] else "NO")
    tally[conc or "no_local_phase"]+=1
    out.write("\t".join(map(str,[ch,r["family"],r["icr"],r["expected_methylated"],r["methylated_hap"],r["hap1_meth"],r["hap2_meth"],len(loc),"" if fl is None else round(fl,2),len(blk),"" if fb is None else round(fb,2),sw,par,conc]))+"\n")
out.close()
print("local-window parent-of-origin:",dict(tally)," phase-switch flags:",switches)
for ln in open(ME+"icr_parent_of_origin_local.tsv"):
    f=ln.rstrip("\n").split("\t")
    if f[-1]=="NO" or f[11]=="SWITCH": print("  ",f[0],f[2],"meth=",f[4],"local n=%s frac=%s | block n=%s frac=%s"%(f[7],f[8],f[9],f[10]),f[11],"->",f[12],f[13])
