#!/usr/bin/env python3
"""Matched-depth internal null, and the training dose-response, in PEP space.

The internal null needs no ground-truth model and no entrapment extrapolation:
walking the pooled target+decoy list best-first, count accepted entrapment
targets against accepted entrapment decoys.  Both populations are certainly
false, so under target-decoy exchangeability the ratio is 1 at every depth.
Comparing scorings at a matched number of accepted entrapment decoys removes
the sensitivity difference between them.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps
from pep_rootcause_controls import read_pin, read_results

DEPTHS = [2, 5, 10, 25, 50, 100, 133, 250, 500, 1000, 2500]


def curve(score, label, pure):
    """Cumulative entrapment targets at each accepted-entrapment-decoy count."""
    o = np.argsort(-score, kind="stable")
    t = label[o] > 0
    p = pure[o]
    ct = np.cumsum(t & p)
    cd = np.cumsum((~t) & p)
    out = {}
    for d in DEPTHS:
        j = np.searchsorted(cd, d)
        out[d] = float(ct[j] / d) if j < cd.size else None
    return out


def pep_summary(pep, label, pure, mixed):
    t = label > 0
    d = ~t
    fg = float(pure[d & ~mixed].sum()) / max(int((d & ~mixed).sum()), 1)
    out = {"f_global": fg, "n_targets": int(t.sum()),
           "min_pep": float(np.nanmin(pep[t])),
           "n_pep_lt_1e3": int((pep[t] < 1e-3).sum())}
    for lo, hi, name in ((0, 1e-3, "lt_1e-3"), (0, 1e-2, "lt_1e-2"), (0, 5e-2, "lt_5e-2")):
        mt = t & (pep >= lo) & (pep < hi)
        md = d & (pep >= lo) & (pep < hi)
        et, ed = int(pure[mt].sum()), int(pure[md].sum())
        out[name] = {"n": int(mt.sum()),
                     "mean_pep": float(pep[mt].mean()) if mt.any() else None,
                     "ent_target": et, "ent_decoy": ed,
                     "ratio": (et / ed) if ed else None,
                     "obs_adj": (et / mt.sum() / fg) if mt.any() else None}
        if out[name]["obs_adj"] is not None and out[name]["mean_pep"]:
            out[name]["calib_ratio"] = out[name]["obs_adj"] / out[name]["mean_pep"]
    return out


def load_run(dirs):
    S, P, L, M, PE = [], [], [], [], []
    for d in dirs:
        for decoy, fn in ((0, "target.tsv"), (1, "decoy.tsv")):
            f = Path(d) / fn
            if not f.exists():
                continue
            s, _, pep, pu, mx = read_results(f, decoy)
            S.append(s); PE.append(pep); P.append(pu); M.append(mx)
            L.append(np.full(s.size, -1 if decoy else 1, np.int8))
    return (np.concatenate(S), np.concatenate(L), np.concatenate(P),
            np.concatenate(M), np.concatenate(PE))


def main():
    pin_root = Path(sys.argv[1]); maxiter_root = Path(sys.argv[2])
    cur_root = Path(sys.argv[3]); cpp_root = Path(sys.argv[4])
    stems = sorted(d.name for d in pin_root.iterdir()
                   if d.is_dir() and d.name.startswith("comet-") and d.name != "comet-out")
    out = {"depths": DEPTHS, "scorings": []}

    # raw search scores
    for col, sign, name in (("Xcorr", 1.0, "raw Comet XCorr"),
                            ("lnExpect", -1.0, "raw Comet -ln(E-value)")):
        S, L, P, M, PEP = [], [], [], [], []
        for st in stems:
            s, l, p, m = read_pin(pin_root / st / "comet.pin", col, sign)
            _, pep = qvalues_and_peps(s, l)
            S.append(s); L.append(l); P.append(p); M.append(m); PEP.append(pep)
        s = np.concatenate(S); l = np.concatenate(L)
        p = np.concatenate(P); m = np.concatenate(M); pep = np.concatenate(PEP)
        out["scorings"].append({"name": name, "curve": curve(s, l, p),
                                "pep": pep_summary(pep, l, p, m)})

    # percolator-rs training dose-response (seeds 1-3 pooled per maxiter)
    for mi in (0, 1, 2, 3, 5, 10):
        dirs = sorted(str(x) for x in (maxiter_root / f"mi-{mi}").glob("seed-*/comet-*"))
        if not dirs:
            continue
        s, l, p, m, pep = load_run(dirs)
        # curve must be per-scoring-list; pooling three seeds of the same six
        # files triples every count, which the ratio is invariant to, but the
        # matched depths are not -- so scale the depths by the seed count.
        seeds = len({Path(d).parent.name for d in dirs})
        c = curve(s, l, p)
        out["scorings"].append({"name": f"percolator-rs maxiter={mi}", "seeds": seeds,
                                "curve": {k: v for k, v in c.items()},
                                "curve_depths_are_pooled_over_seeds": True,
                                "pep": pep_summary(pep, l, p, m)})

    # current default run, seed 1 only, for a like-for-like curve
    dirs = sorted(str(x) for x in (cur_root / "seed-1").glob("comet-*"))
    s, l, p, m, pep = load_run(dirs)
    out["scorings"].append({"name": "percolator-rs default, seed 1",
                            "curve": curve(s, l, p), "pep": pep_summary(pep, l, p, m)})

    # reference C++ under matched competition
    dirs = []
    for st in stems:
        f = cpp_root / f"{st}.target.tsv"
        if f.exists():
            dirs.append(st)
    S, L, P, M, PEP = [], [], [], [], []
    for st in dirs:
        for decoy, suffix in ((0, "target"), (1, "decoy")):
            s, _, pep, pu, mx = read_results(cpp_root / f"{st}.{suffix}.tsv", decoy)
            S.append(s); PEP.append(pep); P.append(pu); M.append(mx)
            L.append(np.full(s.size, -1 if decoy else 1, np.int8))
    s = np.concatenate(S); l = np.concatenate(L); p = np.concatenate(P)
    m = np.concatenate(M); pep = np.concatenate(PEP)
    out["scorings"].append({"name": "C++ Percolator 3.09 (qvality), TDC",
                            "curve": curve(s, l, p), "pep": pep_summary(pep, l, p, m)})
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
