# -*- coding: utf-8 -*-
"""Syndrome-match scoring: how well does each lrWGS proband match each published
syndrome phenotype, and does the true SV carrier rank at the top of their own?

Published feature sets (union):
  FBRSL1  Ummat et al. 2020, PMC7519918, Table 1 (3 de novo truncating patients)
  MECP2   Neul et al. 2010 revised Rett criteria (4 main, 11 supportive, regression)
  DNMT3A  Tatton-Brown et al. 2018 TBRS, 55 individuals (+ GeneReviews)
  CELSR1  DDG2P: hereditary lymphoedema (monoallelic); NDD only BIALLELIC (2025)

Coding of a SPARK field for one proband:
  present  value in {1, True, yes}, or an ordinal/numeric threshold met
  absent   value in {0, False, no}, OR a blank checkbox under a gate that IS positive
           (item was shown and left unticked)
  unknown  'na_survey_logic' (= NOT ASKED per the data dictionary), blank with no positive
           gate, or no row in that instrument. Never inferred as absent.
A feature with no SPARK field at all is 'not in SPARK' and excluded from scoring.

Score_S(p) = present / (present + absent) over features of syndrome S assessable in p.
Rank the carrier among all lrWGS probands on their own syndrome; empirical
p = fraction of probands scoring >= the carrier (carrier included).
"""
import csv
import sys, os
from collections import OrderedDict, defaultdict

P = os.environ.get("PHENO_DIR") or sys.exit("set PHENO_DIR to the SPARK phenotype release directory")
L = os.path.dirname(os.path.abspath(__file__))
ours = sorted(set(x.strip() for x in open(os.environ.get("PROBANDS", "probands.txt")) if x.strip() and not x.startswith("#")))
CARRIER = OrderedDict(l.split()[:2] for l in open(os.environ.get("CARRIERS", "carriers.tsv")) if l.strip() and not l.startswith("#"))
NA = "na_survey_logic"


def load(table, keep):
    out = {}
    f = os.path.join(P, table + "-" + os.environ.get("PHENO_RELEASE", "2026-06-25") + ".csv")
    with open(f, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            s = r.get("subject_sp_id")
            if s in keep:
                out[s] = r
    return out


keep = set(ours)
med = load("basic_medical_screening", keep)
core = load("core_descriptive_variables", keep)
rbsr = load("rbsr", keep)
bhc = load("background_history_child", keep)
dcdq = load("dcdq", keep)
print("lrWGS probands: %d | with medical screening %d | core %d | rbsr %d | bhc %d | dcdq %d"
      % (len(ours), len(med), len(core), len(rbsr), len(bhc), len(dcdq)))

# gate -> sub-item checkbox families (blank under a positive gate = absent)
GATES = {"birth_def_": "med_cond_birth_def", "growth_": "med_cond_growth",
         "visaud_": "med_cond_visaud", "neuro_": "med_cond_neuro",
         "mood_": "mood_or_anx", "behav_": "attn_behav", "gen_test_": "gen_test"}


def truthy(v):
    return (v or "").strip().lower() in ("1", "true", "yes")


def falsy(v):
    return (v or "").strip().lower() in ("0", "false", "no")


def code_flag(row, field):
    """present/absent/unknown for a yes/no medical-screening style field."""
    if row is None or field not in row:
        return "unknown"
    v = (row.get(field) or "").strip()
    if truthy(v):
        return "present"
    if falsy(v):
        return "absent"
    if v == NA:
        return "unknown"   # dictionary: na_survey_logic = NOT ASKED (survey logic); never infer absence
    if v == "":
        for pre, gate in GATES.items():
            if field.startswith(pre) and field != gate:
                g = (row.get(gate) or "").strip()
                if truthy(g):
                    return "absent"   # shown under a positive gate and left unticked
                return "unknown"      # gate blank / not asked -> item never shown
        return "unknown"
    return "unknown"


def any_present(row, fields):
    codes = [code_flag(row, f) for f in fields]
    if "present" in codes:
        return "present"
    if all(c == "absent" for c in codes):
        return "absent"
    return "unknown" if "unknown" in codes else "absent"


def num(row, field):
    try:
        return float((row or {}).get(field, ""))
    except (TypeError, ValueError):
        return None


def code_num(row, field, op, thr):
    v = num(row, field)
    if v is None:
        return "unknown"
    return "present" if (v >= thr if op == ">=" else v <= thr) else "absent"


LANG = {"0_no_words": 0, "1_single_words": 1, "2_short_sentences": 2, "3_long_sentences": 3}


def code_lang_le1(s):
    v = (core.get(s, {}).get("language_level_latest") or "").strip()
    if v not in LANG:
        return "unknown"
    return "present" if LANG[v] <= 1 else "absent"


def F(s):
    m, c, r, b, d = med.get(s), core.get(s), rbsr.get(s), bhc.get(s), dcdq.get(s)
    return {
        # ---- shared
        "Intellectual disability / GDD": any_present(m, ["dev_id"]) if m else (
            "present" if truthy(c.get("cognitive_impairment_latest")) else
            "absent" if falsy(c.get("cognitive_impairment_latest")) else "unknown") if c else "unknown",
        "Delayed speech / language": any_present(m, ["dev_lang", "dev_speech", "dev_lang_dis"]),
        "Minimally verbal (<= single words)": code_lang_le1(s),
        "Autistic behaviour": "present",   # cohort-wide by ascertainment; shown, excluded from score
        "Seizures / epilepsy": code_flag(m, "neuro_sz"),
        "Motor delay / abnormal gait": any_present(m, ["dev_motor"]),
        "Developmental regression": ("present" if (truthy((c or {}).get("regress_lang_y_n")) or truthy((c or {}).get("regress_other_y_n")))
                                     else "absent" if (falsy((c or {}).get("regress_lang_y_n")) and falsy((c or {}).get("regress_other_y_n")))
                                     else "unknown"),
        "Sleep disturbance": any_present(m, ["sleep_probs", "sleep_dx"]),
        "Feeding / swallowing difficulty": any_present(m, ["eating_probs", "feeding_dx"]),
        "Facial dysmorphism": code_flag(m, "birth_def_fac"),
        "Hearing impairment": code_flag(m, "visaud_deaf"),
        "Heart defect": code_flag(m, "birth_def_thorac_heart"),
        "Cleft palate": code_flag(m, "birth_def_cleft_palate"),
        "Skeletal anomaly (any)": code_flag(m, "birth_def_bone"),
        "Kyphoscoliosis / spine": code_flag(m, "birth_def_bone_spine"),
        "Microcephaly": code_flag(m, "growth_microceph"),
        "Macrocephaly": code_flag(m, "growth_macroceph"),
        "Postnatal growth retardation / short": any_present(m, ["growth_short", "growth_low_wt"]),
        "Obesity": code_flag(m, "growth_obes"),
        "Neonatal respiratory support": code_flag(m, "birth_oxygen"),
        "Hand stereotypies (RBS-R hand/finger)": code_num(r, "q03_hand_finger", ">=", 1),
        "Brain malformation": code_flag(m, "birth_def_cns_brain"),
        "Behavioural / psychiatric (ADHD, anxiety)": any_present(m, ["behav_adhd", "attn_behav", "mood_or_anx"]),
        "Strabismus": code_flag(m, "visaud_strab"),
        "Polydactyly": code_flag(m, "birth_def_bone_polydact"),
    }


# which published features belong to which syndrome, and what is NOT in SPARK
SYN = OrderedDict()
SYN["FBRSL1"] = (["Intellectual disability / GDD", "Delayed speech / language", "Microcephaly",
                  "Feeding / swallowing difficulty", "Postnatal growth retardation / short",
                  "Skeletal anomaly (any)", "Heart defect", "Cleft palate", "Neonatal respiratory support",
                  "Hearing impairment", "Facial dysmorphism"],
                 ["Camptodactyly / contractures", "Asplenia", "Anal anomaly", "Skin creases"])
SYN["MECP2"] = (["Developmental regression", "Minimally verbal (<= single words)",
                 "Motor delay / abnormal gait", "Hand stereotypies (RBS-R hand/finger)",
                 "Sleep disturbance", "Kyphoscoliosis / spine", "Postnatal growth retardation / short",
                 "Seizures / epilepsy", "Feeding / swallowing difficulty"],
                ["Breathing disturbance awake", "Bruxism", "Abnormal tone", "Vasomotor disturbance",
                 "Small cold hands/feet", "Laughing/screaming spells", "Diminished pain response",
                 "Eye pointing"])
SYN["DNMT3A"] = (["Macrocephaly", "Intellectual disability / GDD", "Facial dysmorphism", "Obesity",
                  "Motor delay / abnormal gait", "Kyphoscoliosis / spine", "Seizures / epilepsy",
                  "Behavioural / psychiatric (ADHD, anxiety)", "Skeletal anomaly (any)"],
                 ["Tall stature", "Joint hypermobility", "Hypotonia"])
SYN["CELSR1"] = (["Brain malformation", "Intellectual disability / GDD", "Seizures / epilepsy",
                  "Behavioural / psychiatric (ADHD, anxiety)"],
                 ["Lymphoedema (monoallelic G2P phenotype)", "Neural tube defect"])

feat = {s: F(s) for s in ours}


def score(s, syn):
    fs, _ = SYN[syn]
    pres = sum(1 for f in fs if feat[s][f] == "present")
    abse = sum(1 for f in fs if feat[s][f] == "absent")
    n = pres + abse
    return (pres / float(n) if n else None, pres, n)


print("")
print("=" * 96)
print("CARRIER MATCH TO OWN SYNDROME, ranked among all %d lrWGS probands" % len(ours))
print("=" * 96)
print("%-8s %-10s %-11s %-6s %-14s %-9s %s" % ("gene", "carrier", "score", "feats", "rank", "emp_p", "cohort median score"))
for s, g in CARRIER.items():
    sc = {p: score(p, g) for p in ours}
    mine = sc[s]
    vals = [v[0] for v in sc.values() if v[0] is not None]
    if mine[0] is None:
        print("%-8s %-10s  not assessable (no features codable)" % (g, s))
        continue
    ge = sum(1 for v in vals if v >= mine[0] - 1e-9)
    vals_s = sorted(vals, reverse=True)
    rank = vals_s.index(mine[0]) + 1
    med_ = sorted(vals)[len(vals) // 2]
    print("%-8s %-10s %5.2f       %d/%-3d  %2d of %-9d %-9.3f %.2f"
          % (g, s, mine[0], mine[1], mine[2], rank, len(vals), ge / float(len(vals)), med_))

print("")
print("=" * 96)
print("SPECIFICITY: each carrier's score on ALL four syndromes (row = carrier, col = syndrome)")
print("=" * 96)
print("%-10s %-8s " % ("carrier", "gene") + " ".join("%8s" % g for g in SYN))
for s, g in CARRIER.items():
    print("%-10s %-8s " % (s, g) + " ".join(("%8.2f" % score(s, k)[0]) if score(s, k)[0] is not None else "%8s" % "-" for k in SYN))

print("")
print("=" * 96)
print("UNION FEATURE TABLE  (+ present, - absent, U unknown/not assessed; * = published core feature of that column's syndrome)")
print("=" * 96)
allf = []
for g, (fs, _) in SYN.items():
    for f in fs:
        if f not in allf:
            allf.append(f)
sym = {"present": "+", "absent": "-", "unknown": "U"}   # U = not assessed (distinct from absent)
print("%-42s " % "feature" + " ".join("%-8s" % g for g in CARRIER.values()))
for f in allf:
    cells = []
    for s, g in CARRIER.items():
        star = "*" if f in SYN[g][0] else " "
        cells.append("%s%s" % (sym[feat[s][f]], star))
    print("%-42s " % f[:42] + " ".join("%-8s" % c for c in cells))
print("")
print("published features with NO SPARK field (excluded from scoring):")
for g, (_, nofield) in SYN.items():
    print("  %-8s %s" % (g, "; ".join(nofield)))

out = os.path.join(L, "syndrome_match_table.tsv")
with open(out, "w") as fh:
    fh.write("feature\t" + "\t".join("%s_%s" % (g, s) for s, g in CARRIER.items()) + "\n")
    for f in allf:
        fh.write(f + "\t" + "\t".join(feat[s][f] for s in CARRIER) + "\n")
print("wrote " + out)
