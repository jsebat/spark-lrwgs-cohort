from phase_dnm.classify import likelihood as L

P = L.LikParams()


def rows(A, O, T, U, N1=(0, 10), N2=(0, 10), UNT=(0, 0)):
    return {"A": A, "O": O, "T": T, "U": U, "N1": N1, "N2": N2, "UNT": UNT}


def best(post):
    return max(L.HYPS, key=lambda h: post["lik_post_" + h])


def test_germline_pattern_wins_even_against_a_strong_artefact_prior():
    r = rows(A=(10, 0), O=(0, 10), T=(0, 10), U=(0, 10))
    p = L.posterior(r, P, po_resolved=True)
    assert best(p) == "germline" and p["phase_score"] > 0.7, p
    # the nearest alternative at ~10 reads on T is a LOW-FRACTION PARENTAL MOSAIC, not an artefact: the reported
    # likelihood ratio is against that alternative and is depth-limited (P10)
    assert p["lik_best_alternative"] == "parental_mosaic" and 0 < p["lik_log10lr_germline"] < 3
    assert p["lik_post_artefact"] < 1e-3 and p["lik_post_inherited"] < 1e-3
    # more depth per haplotype -> more confidence (a low-fraction parental mosaic on T becomes excludable)
    r20 = rows(A=(20, 0), O=(0, 20), T=(0, 20), U=(0, 20), N1=(0, 20), N2=(0, 20))
    p20 = L.posterior(r20, P, True)
    assert p20["phase_score"] > 0.9 and p20["phase_score"] >= p["phase_score"]
    assert p20["lik_log10lr_germline"] > p["lik_log10lr_germline"]
    # and the parental-mosaic posterior shrinks with depth on T
    assert p20["lik_post_parental_mosaic"] < p["lik_post_parental_mosaic"]


def test_inherited_when_transmitted_haplotype_carries_alt():
    r = rows(A=(10, 0), O=(0, 10), T=(9, 1), U=(0, 10))
    p = L.posterior(r, P, po_resolved=True)
    assert best(p) == "inherited" and p["phase_score"] < 0.01


def test_parental_mosaic_intermediate_fraction_on_T():
    r = rows(A=(10, 0), O=(0, 10), T=(4, 16), U=(0, 20))
    p = L.posterior(r, P, po_resolved=True)
    assert best(p) == "parental_mosaic", p
    # one alt read of twenty on T is within error: still germline
    r1 = rows(A=(10, 0), O=(0, 10), T=(1, 19), U=(0, 20))
    assert best(L.posterior(r1, P, True)) == "germline"


def test_child_mosaic_subclonal_on_A():
    r = rows(A=(6, 14), O=(0, 20), T=(0, 12), U=(0, 12), N1=(0, 12), N2=(0, 12))
    p = L.posterior(r, P, po_resolved=True)
    assert best(p) == "child_mosaic", p


def test_artefact_when_alt_is_everywhere():
    r = rows(A=(4, 6), O=(3, 7), T=(3, 7), U=(4, 6), N1=(3, 7), N2=(4, 6))
    p = L.posterior(r, P, po_resolved=True)
    assert best(p) == "artefact" and p["phase_score"] < 0.01


def test_unresolved_transmission_averages_the_pair_and_lowers_confidence():
    r = rows(A=(10, 0), O=(0, 10), T=(0, 10), U=(0, 10))
    resolved = L.posterior(r, P, True)["phase_score"]
    unresolved = L.posterior(r, P, False)["phase_score"]
    assert unresolved <= resolved + 1e-9
    # with the alt on ONE of the two unlabelled parental haplotypes, unresolved cannot tell germline from nothing:
    # inherited must dominate germline just as when resolved
    r2 = rows(A=(10, 0), O=(0, 10), T=(0, 10), U=(9, 1))
    assert best(L.posterior(r2, P, False)) == "inherited"


def test_low_depth_gives_low_confidence_not_a_call():
    # two alt reads on an otherwise clean matrix: the posterior must stay well below the 10-read case; the hard
    # evidence floor (min_alt_reads, P10) lives in the rule layer, and P15 uses both layers
    r = rows(A=(2, 0), O=(0, 2), T=(0, 2), U=(0, 2), N1=(0, 2), N2=(0, 2))
    p2 = L.posterior(r, P, True)["phase_score"]
    p10 = L.posterior(rows(A=(10, 0), O=(0, 10), T=(0, 10), U=(0, 10)), P, True)["phase_score"]
    assert p2 < 0.9 and p2 < p10


def test_rows_from_evidence_uses_orientation_and_transmission():
    ev = {"C1_alt": "10", "C1_ref": "0", "C2_alt": "0", "C2_ref": "10", "F1_alt": "0", "F1_ref": "10", "F2_alt": "9", "F2_ref": "1",
          "M1_alt": "0", "M1_ref": "10", "M2_alt": "0", "M2_ref": "10", "C_untagged_alt": "1", "C_untagged_dp": "2",
          "child_hap1_is": "P", "F_transmitted_hap": "1", "M_transmitted_hap": "2"}
    rows_, po = L.rows_from_evidence(ev)
    assert po and rows_["A"] == (10, 0) and rows_["T"] == (0, 10) and rows_["U"] == (9, 1) and rows_["N1"] == (0, 10)
    s = L.score_evidence_row(ev, P)
    # alt on the father's UNTRANSMITTED haplotype AND on the child's paternal haplotype fits neither germline nor
    # inheritance (the child did not receive that haplotype): the model must not call germline
    assert best(s) == "artefact" and s["phase_score"] < 0.1
    ev["F_transmitted_hap"] = "2"           # now the alt-bearing paternal haplotype IS the transmitted one
    assert best(L.score_evidence_row(ev, P)) == "inherited"
    ev["F2_alt"], ev["F2_ref"] = "0", "10"  # father clean on both haplotypes: germline
    assert best(L.score_evidence_row(ev, P)) == "germline"
    ev["child_hap1_is"] = "."; ev["F_transmitted_hap"] = "."
    s = L.score_evidence_row(ev, P)
    assert s["phase_score"] is not None       # undetermined origin still scores (unordered pair)
