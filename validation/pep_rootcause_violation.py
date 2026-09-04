#!/usr/bin/env python3
"""Which assumption violation reproduces the observed PEP optimism?

Each scenario keeps the estimator untouched and breaks exactly one assumption of
the target-decoy contract in src/stats.rs.  Ground truth is known per row, so no
entrapment extrapolation enters anywhere.

  exchangeable   nothing broken (control)
  foreign        a class of INCORRECT targets drawn from a shifted distribution
                 with no decoy counterpart -- the homology channel, i.e. decoys
                 under-count incorrect targets in the tail
  decoy_shift    decoys drawn from a distribution shifted below the incorrect
                 targets they are supposed to mirror -- provenance leakage
  decoy_narrow   decoys share the incorrect-target mean but have a lighter tail
  clustered      exchangeable scores, but PSMs arrive in correlated clusters
                 (one peptide contributing many spectra) -- dependence only
  ties           exchangeable scores rounded onto a coarse grid -- discreteness
  trained        scores produced by a discriminant fitted on the same data
                 through 3-fold cross-validation -- semi-supervised selection
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])
PHI = lambda s, m=0.0, sd=1.0: np.exp(-0.5 * ((s - m) / sd) ** 2) / sd


def competition(t, d):
    win = t >= d
    return np.where(win, t, d), np.where(win, 1, -1).astype(np.int8), win


def scenario(rng, name, n, pi0, mu, **kw):
    """Return score, label, false_target_flag, true_lfdr (or None)."""
    if name == "exchangeable":
        null = rng.random(n) < pi0
        t = np.where(null, rng.normal(0, 1, n), rng.normal(mu, 1, n))
        d = rng.normal(0, 1, n)
        s, lab, win = competition(t, d)
        false_t = win & null
        truth = pi0 * PHI(s) / (pi0 * PHI(s) + (1 - pi0) * PHI(s, mu))
        return s, lab, false_t, truth

    if name == "foreign":
        ph, muh = kw["p_foreign"], kw["mu_foreign"]
        u = rng.random(n)
        cls = np.where(u < pi0, 0, np.where(u < pi0 + ph, 1, 2))  # 0 null,1 foreign,2 correct
        t = np.where(cls == 0, rng.normal(0, 1, n),
             np.where(cls == 1, rng.normal(muh, 1, n), rng.normal(mu, 1, n)))
        d = rng.normal(0, 1, n)
        s, lab, win = competition(t, d)
        false_t = win & (cls != 2)
        pc = 1 - pi0 - ph
        num = pi0 * PHI(s) + ph * PHI(s, muh)
        truth = num / (num + pc * PHI(s, mu))
        return s, lab, false_t, truth

    if name in ("decoy_shift", "decoy_narrow"):
        null = rng.random(n) < pi0
        t = np.where(null, rng.normal(0, 1, n), rng.normal(mu, 1, n))
        if name == "decoy_shift":
            d = rng.normal(-kw["delta"], 1, n)
        else:
            d = rng.normal(0, kw["sd"], n)
        s, lab, win = competition(t, d)
        false_t = win & null
        return s, lab, false_t, None

    if name == "clustered":
        m = kw["cluster"]
        g = n // m
        null_g = rng.random(g) < pi0
        base_g = np.where(null_g, rng.normal(0, 1, g), rng.normal(mu, 1, g))
        null = np.repeat(null_g, m)
        t = np.repeat(base_g, m) + rng.normal(0, kw["jitter"], g * m)
        dg = rng.normal(0, 1, g)
        d = np.repeat(dg, m) + rng.normal(0, kw["jitter"], g * m)
        s, lab, win = competition(t, d)
        return s, lab, win & null, None

    if name == "ties":
        null = rng.random(n) < pi0
        t = np.where(null, rng.normal(0, 1, n), rng.normal(mu, 1, n))
        d = rng.normal(0, 1, n)
        s, lab, win = competition(t, d)
        return np.round(s / kw["grid"]) * kw["grid"], lab, win & null, None

    if name == "trained":
        # p informative features plus q pure-noise features; the decoy and the
        # incorrect target are exchangeable in every feature, so any optimism
        # here comes from fitting the discriminant, not from the data.
        p, qn, folds = kw["informative"], kw["noise"], 3
        null = rng.random(n) < pi0
        xt = rng.normal(0, 1, (n, p + qn))
        xt[~null, :p] += mu / np.sqrt(p)
        xd = rng.normal(0, 1, (n, p + qn))
        yt, yd = _cv_scores(rng, xt, xd, folds, kw.get("iters", 3))
        s, lab, win = competition(yt, yd)
        return s, lab, win & null, None

    raise ValueError(name)


def _cv_scores(rng, xt, xd, folds, iters):
    """Percolator-shaped semi-supervised loop, cross-validated by fold."""
    n = xt.shape[0]
    fold = rng.integers(0, folds, n)
    yt = np.zeros(n)
    yd = np.zeros(n)
    for f in range(folds):
        tr = fold != f
        te = ~tr
        x = np.vstack([xt[tr], xd[tr]])
        lab = np.concatenate([np.ones(tr.sum()), -np.ones(tr.sum())])
        w = x[: tr.sum()].mean(0) - x[tr.sum():].mean(0)   # initial direction
        for _ in range(iters):
            sc = x @ w
            qv, _ = qvalues_and_peps(sc, lab.astype(np.int8))
            pos = (lab > 0) & (qv < 0.01)
            neg = lab < 0
            if pos.sum() < 10:
                break
            a, b = x[pos].mean(0), x[neg].mean(0)
            cov = np.cov(np.vstack([x[pos], x[neg]]).T) + 1e-6 * np.eye(x.shape[1])
            w = np.linalg.solve(cov, a - b)                 # regularized LDA step
        yt[te] = xt[te] @ w
        yd[te] = xd[te] @ w
    return yt, yd


def table(pep, false_flag, truth=None):
    idx = np.digitize(pep, BINS) - 1
    rows = []
    for b in range(len(BINS) - 1):
        m = idx == b
        if not m.any():
            continue
        r = {"lo": float(BINS[b]), "hi": float(BINS[b + 1]), "n": int(m.sum()),
             "mean_pep": float(pep[m].mean()), "obs_false": float(false_flag[m].mean())}
        r["gap"] = r["obs_false"] - r["mean_pep"]
        r["ratio"] = r["obs_false"] / r["mean_pep"] if r["mean_pep"] > 0 else None
        if truth is not None:
            r["mean_true_lfdr"] = float(truth[m].mean())
        rows.append(r)
    n = sum(x["n"] for x in rows)
    return rows, (sum(x["n"] * x["gap"] for x in rows) / n,
                  sum(x["n"] * abs(x["gap"]) for x in rows) / n)


def run(rng, name, reps=6, n=200000, pi0=0.5, mu=3.0, **kw):
    peps, fl, tr = [], [], []
    for _ in range(reps):
        s, lab, false_t, truth = scenario(rng, name, n, pi0, mu, **kw)
        _, pep = qvalues_and_peps(s, lab)
        t = lab > 0
        peps.append(pep[t]); fl.append(false_t[t].astype(float))
        if truth is not None:
            tr.append(truth[t])
    pep = np.concatenate(peps); fl = np.concatenate(fl)
    truth = np.concatenate(tr) if tr else None
    rows, (signed, absolute) = table(pep, fl, truth)
    # internal decoy-based null: at the same score threshold, how many false
    # targets per decoy?  1.0 under exchangeability, no ground truth needed.
    return {"scenario": name, "reps": reps, "n": n, "pi0": pi0, "mu": mu, **kw,
            "weighted_signed": signed, "weighted_abs": absolute,
            "n_targets": int(pep.size), "sum_pep": float(pep.sum()),
            "n_false": int(fl.sum()),
            "sum_pep_over_false": float(pep.sum() / max(fl.sum(), 1)),
            "bins": rows}


def main():
    rng = np.random.default_rng(90210)
    out = []
    out.append(run(rng, "exchangeable"))
    for ph, muh in ((0.02, 2.5), (0.05, 2.5), (0.02, 3.0), (0.05, 3.0), (0.10, 3.0)):
        out.append(run(rng, "foreign", p_foreign=ph, mu_foreign=muh, pi0=0.45))
    for delta in (0.05, 0.1, 0.2, 0.4):
        out.append(run(rng, "decoy_shift", delta=delta))
    for sd in (0.9, 0.95):
        out.append(run(rng, "decoy_narrow", sd=sd))
    for m, j in ((5, 0.3), (20, 0.3), (20, 0.05)):
        out.append(run(rng, "clustered", cluster=m, jitter=j))
    for grid in (0.01, 0.1, 0.5):
        out.append(run(rng, "ties", grid=grid))
    for it in (0, 1, 3, 10):
        out.append(run(rng, "trained", reps=3, n=100000, informative=4, noise=16, iters=it))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
