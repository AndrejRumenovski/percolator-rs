#!/usr/bin/env python3
"""Two controls for the real-data PEP calibration.

A. RAW SEARCH SCORE.  The same estimator, same entrapment labels, same
   spectrum-level competition, but the score is Comet's own XCorr or ln(E-value)
   instead of the rescored value.  The learner is out of the loop entirely.

B. REFERENCE IMPLEMENTATION.  The C++ Percolator 3.09 default PEP, produced by
   its score-aware monotone I-spline fit (not standalone QVALITY), on the
   identical PIN files under --post-processing-tdc so the reported set is one
   PSM per spectrum, as percolator-rs reports.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps

csv.field_size_limit(10 ** 9)
BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])
NAMES = ["[0,1e-4)", "[1e-4,1e-3)", "[1e-3,5e-3)", "[5e-3,.01)", "[.01,.02)",
         "[.02,.05)", "[.05,.10)", "[.10,.20)", "[.20,.50)", "[.50,1]"]


def classify(members, decoy):
    if decoy:
        members = [m.removeprefix("DECOY_").removeprefix("decoy_") for m in members]
    pure = bool(members) and all(m.startswith("ENT_") for m in members)
    mixed = any(m.startswith("ENT_") for m in members) and not pure
    return pure, mixed


def read_pin(path, score_col, sign=1.0):
    """Spectrum-level TDC competition on one raw PIN feature.

    `sign` orients the column so that larger is always better; it is applied
    before the competition, not after, or the winner would be the worst
    candidate for a column like lnExpect where small is good."""
    best = {}
    with open(path, newline="") as h:
        head = h.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(head)}
        si, mi, li, ni = ci["ScanNr"], ci["ExpMass"], ci["Label"], ci[score_col]
        pi = ci["Peptide"]
        for line in h:
            f = line.rstrip("\n").split("\t")
            # Match Dataset::spectrum_key in production: a scan can contain
            # multiple precursor masses, and those are separate competitions.
            key = (f[si], f[mi])
            v = sign * float(f[ni])
            cur = best.get(key)
            if cur is None or v > cur[0]:
                best[key] = (v, int(f[li]), f[pi + 1:])
    score = np.empty(len(best)); label = np.empty(len(best), np.int8)
    pure = np.empty(len(best), bool); mixed = np.empty(len(best), bool)
    for i, (v, lab, prot) in enumerate(best.values()):
        mem = [x for p in prot for x in p.replace(";", "\t").split("\t") if x]
        pu, mx = classify(mem, lab < 0)
        score[i] = v; label[i] = 1 if lab > 0 else -1; pure[i] = pu; mixed[i] = mx
    return score, label, pure, mixed


def read_results(path, decoy):
    score, q, pep, pure, mixed = [], [], [], [], []
    with open(path, newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            vals = []
            for k, v in row.items():
                if k == "proteinIds" or k is None:
                    if isinstance(v, list):
                        vals.extend(v)
                    elif v:
                        vals.append(v)
            mem = [x for v in vals for x in v.replace(";", "\t").split("\t") if x]
            pu, mx = classify(mem, decoy)
            score.append(float(row["score"])); q.append(float(row["q-value"]))
            p = row["posterior_error_prob"]
            pep.append(float("nan") if p == "NA" else float(p))
            pure.append(pu); mixed.append(mx)
    return (np.array(score), np.array(q), np.array(pep),
            np.array(pure, bool), np.array(mixed, bool))


def table(pep, pure, mixed, label, f_global):
    """Calibration by PEP bin.  Predicted vs entrapment-implied observed, plus
    the adjustment-free internal null (entrapment targets per entrapment decoy)."""
    t = label > 0
    idx = np.digitize(pep, BINS) - 1
    rows = []
    for b in range(len(BINS) - 1):
        mt = t & (idx == b); md = (~t) & (idx == b)
        n = int(mt.sum())
        if n == 0:
            continue
        et = int(pure[mt].sum()); ed = int(pure[md].sum())
        nd = int((~pure[md] & ~mixed[md]).sum())
        fb = ed / (ed + nd) if (ed + nd) >= 20 else f_global
        mp = float(pep[mt].mean()); raw = et / n
        rows.append({"bin": NAMES[b], "n": n, "mean_pep": mp,
                     "n_ent_target": et, "n_ent_decoy": ed, "n_nat_decoy": nd,
                     "obs_f1": raw, "obs_adj": raw / fb if fb > 0 else None,
                     "gap_f1": raw - mp,
                     "gap_adj": (raw / fb - mp) if fb > 0 else None,
                     "ratio_adj": (raw / fb / mp) if (fb > 0 and mp > 0) else None,
                     "ent_t_over_ent_d": (et / ed) if ed > 0 else None})
    n = sum(r["n"] for r in rows)
    summary = {"n_targets": n,
               "weighted_signed_adj": sum(r["n"] * (r["gap_adj"] or 0) for r in rows) / n,
               "weighted_signed_f1": sum(r["n"] * r["gap_f1"] for r in rows) / n,
               "weighted_abs_adj": sum(r["n"] * abs(r["gap_adj"] or 0) for r in rows) / n}
    return rows, summary


def main():
    pin_root = Path(sys.argv[1])
    cpp_root = Path(sys.argv[2])
    out = {}
    stems = sorted(d.name for d in pin_root.iterdir()
                   if d.is_dir() and d.name.startswith("comet-") and d.name != "comet-out")

    # ---- A. raw search score ----
    for col, sign in (("Xcorr", 1.0), ("lnExpect", -1.0)):
        S, L, P, M = [], [], [], []
        for st in stems:
            s, l, p, m = read_pin(pin_root / st / "comet.pin", col, sign)
            S.append(s); L.append(l); P.append(p); M.append(m)
        # per-file estimation, then pool (matches how percolator-rs is run)
        peps, labels, pures, mixeds = [], [], [], []
        for s, l, p, m in zip(S, L, P, M):
            _, pep = qvalues_and_peps(s, l)
            peps.append(pep); labels.append(l); pures.append(p); mixeds.append(m)
        pep = np.concatenate(peps); lab = np.concatenate(labels)
        pure = np.concatenate(pures); mixed = np.concatenate(mixeds)
        d = lab < 0
        fg = float(pure[d & ~mixed].sum()) / max(int((d & ~mixed).sum()), 1)
        rows, summary = table(pep, pure, mixed, lab, fg)
        out[f"raw_{col}"] = {"f_global": fg, "summary": summary, "bins": rows}

    # ---- B. reference C++ PEPs under matched competition ----
    peps, labels, pures, mixeds = [], [], [], []
    for st in stems:
        for decoy, suffix in ((0, "target"), (1, "decoy")):
            f = cpp_root / f"{st}.{suffix}.tsv"
            if not f.exists():
                continue
            _, _, pep, pu, mx = read_results(f, decoy)
            peps.append(pep); pures.append(pu); mixeds.append(mx)
            labels.append(np.full(pep.size, -1 if decoy else 1, np.int8))
    pep = np.concatenate(peps); lab = np.concatenate(labels)
    pure = np.concatenate(pures); mixed = np.concatenate(mixeds)
    d = lab < 0
    fg = float(pure[d & ~mixed].sum()) / max(int((d & ~mixed).sum()), 1)
    rows, summary = table(pep, pure, mixed, lab, fg)
    out["cpp_percolator_default_ispline_tdc"] = {"f_global": fg, "summary": summary, "bins": rows,
                              "n_pep_at_floor_1e-10": int((pep[lab > 0] <= 1e-10).sum()),
                              "min_pep_target": float(np.nanmin(pep[lab > 0]))}
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
