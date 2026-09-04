#!/usr/bin/env python3
"""Identical score/label inputs, two PEP estimators.

Feeds the *same* percolator-rs entrapment scores to QVALITY 3.09 (the intended
reference methodology: a monotone logistic/spline fit to the target-vs-null
score distributions) and to the production percolator-rs estimator, then
compares the two on the same rows.  Any difference is attributable to the
estimator alone: the scores, the labels, the competition and the entrapment
ground truth are byte-identical between the arms.
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])
NAMES = ["[0,1e-4)", "[1e-4,1e-3)", "[1e-3,5e-3)", "[5e-3,.01)", "[.01,.02)",
         "[.02,.05)", "[.05,.10)", "[.10,.20)", "[.20,.50)", "[.50,1]"]
QVALITY = os.path.expanduser("~/opt/percolator-root/usr/bin/qvality")


def run_qvality(target, null, evaluate_at, workdir):
    """Fit QVALITY on (target scores, null scores) and read the fitted PEP off
    its output grid at `evaluate_at`.  The mixed model must receive the target
    scores only, exactly as Percolator calls it internally."""
    t = Path(workdir) / "t.txt"; n = Path(workdir) / "n.txt"
    o = Path(workdir) / "o.txt"
    np.savetxt(t, target, fmt="%.10f")
    np.savetxt(n, null, fmt="%.10f")
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = os.path.expanduser("~/opt/perc-libs") + ":" + env.get("LD_LIBRARY_PATH", "")
    r = subprocess.run([QVALITY, "-o", str(o), str(t), str(n)],
                       capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    grid = np.loadtxt(o, skiprows=1)
    # qvality reports (score, PEP, q) on its own grid; map each target score to
    # the PEP of the nearest grid score at or below it (the fit is monotone).
    gs, gp = grid[:, 0], grid[:, 1]
    order = np.argsort(gs)
    gs, gp = gs[order], gp[order]
    idx = np.clip(np.searchsorted(gs, evaluate_at, side="right") - 1, 0, gs.size - 1)
    return gp[idx], grid, r.stderr


def calib(pep, pure, mixed, is_target, f_global):
    idx = np.digitize(pep, BINS) - 1
    rows = []
    for b in range(len(BINS) - 1):
        mt = is_target & (idx == b); md = (~is_target) & (idx == b)
        n = int(mt.sum())
        if n == 0:
            continue
        et, ed = int(pure[mt].sum()), int(pure[md].sum())
        nd = int((~pure[md] & ~mixed[md]).sum())
        fb = ed / (ed + nd) if (ed + nd) >= 20 else f_global
        mp = float(pep[mt].mean()); raw = et / n
        rows.append({"bin": NAMES[b], "n": n, "mean_pep": mp, "ent_target": et,
                     "ent_decoy": ed, "obs_adj": raw / fb, "gap": raw / fb - mp,
                     "ratio": (raw / fb / mp) if mp > 0 else None,
                     "ent_t_over_ent_d": (et / ed) if ed else None})
    n = sum(r["n"] for r in rows)
    return rows, {"n": n,
                  "weighted_signed": sum(r["n"] * r["gap"] for r in rows) / n,
                  "weighted_abs": sum(r["n"] * abs(r["gap"]) for r in rows) / n,
                  "sum_pep": float(pep[is_target].sum())}


def main():
    npz = dict(np.load(sys.argv[1], allow_pickle=True))
    seed = 1
    out = {"datasets": []}
    pooled = {k: [] for k in ("pep_rs", "pep_qv", "pure", "mixed", "is_t")}
    with tempfile.TemporaryDirectory(dir=sys.argv[2]) as wd:
        for di, name in enumerate(npz["datasets"]):
            m = (npz["seed"] == seed) & (npz["ds"] == di)
            score = npz["score"][m]; dec = npz["decoy"][m] == 1
            pure = npz["pure"][m]; mixed = npz["mixed"][m]
            pep_rs_reported = npz["pep"][m]
            is_t = ~dec
            # recompute with the production estimator to confirm the reported
            # column is reproduced from (score,label) alone
            lab = np.where(is_t, 1, -1).astype(np.int8)
            _, pep_rs = qvalues_and_peps(score, lab)
            agree = float(np.max(np.abs(pep_rs[is_t] - pep_rs_reported[is_t])))
            pep_qv, grid, err = run_qvality(score[is_t], score[dec], score, wd)
            fg = float(pure[dec & ~mixed].sum()) / max(int((dec & ~mixed).sum()), 1)
            rows_rs, sum_rs = calib(pep_rs, pure, mixed, is_t, fg)
            rows_qv, sum_qv = calib(pep_qv, pure, mixed, is_t, fg)
            out["datasets"].append({
                "dataset": str(name), "n_rows": int(m.sum()),
                "max_abs_reported_minus_recomputed_target_pep": agree,
                "n_decoys": int(dec.sum()), "n_targets": int(is_t.sum()),
                "f_global": fg,
                "qvality_grid_points": int(grid.shape[0]),
                "qvality_min_pep": float(grid[:, 1].min()),
                "rs": {"summary": sum_rs, "bins": rows_rs},
                "qvality": {"summary": sum_qv, "bins": rows_qv},
            })
            for k, v in (("pep_rs", pep_rs), ("pep_qv", pep_qv), ("pure", pure),
                         ("mixed", mixed), ("is_t", is_t)):
                pooled[k].append(v)
            print(f"done {name}", file=sys.stderr)
    P = {k: np.concatenate(v) for k, v in pooled.items()}
    dec = ~P["is_t"]
    fg = float(P["pure"][dec & ~P["mixed"]].sum()) / max(int((dec & ~P["mixed"]).sum()), 1)
    for tag, key in (("rs", "pep_rs"), ("qvality", "pep_qv")):
        rows, summary = calib(P[key], P["pure"], P["mixed"], P["is_t"], fg)
        out[f"pooled_{tag}"] = {"summary": summary, "bins": rows}
    out["pooled_f_global"] = fg
    out["n_targets_total"] = int(P["is_t"].sum())
    out["spearman_like_rank_agreement"] = float(
        np.corrcoef(np.argsort(np.argsort(P["pep_rs"][P["is_t"]])),
                    np.argsort(np.argsort(P["pep_qv"][P["is_t"]])))[0, 1])
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
