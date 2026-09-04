#!/usr/bin/env python3
"""Does semi-supervised training alone break PEP calibration, or does it need a
feature that carries sequence provenance?

Both arms use the identical learner, identical fold structure and identical
number of iterations.  They differ in one feature only:

  clean      every feature is exchangeable between a decoy and an incorrect
             target -- the learner has nothing to exploit
  provenance one extra feature is shifted in decoys relative to incorrect
             targets by `delta`, and carries NO information about whether a
             match is correct (correct and incorrect targets share its mean).
             This is the enzN/enzC channel of a semi-tryptic search with
             reversed decoys, in its simplest possible form.

Ground truth is known per row.  Nothing is fitted to the calibration outcome.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])


def simulate(rng, n, pi0, mu, iters, delta, informative=4, noise=12, folds=3):
    null = rng.random(n) < pi0
    p = informative + noise + 1
    xt = rng.normal(0, 1, (n, p))
    xd = rng.normal(0, 1, (n, p))
    xt[~null, :informative] += mu / np.sqrt(informative)   # correctness signal
    xd[:, -1] -= delta                                     # provenance channel
    yt = np.zeros(n); yd = np.zeros(n)
    fold = rng.integers(0, folds, n)
    for f in range(folds):
        tr, te = fold != f, fold == f
        x = np.vstack([xt[tr], xd[tr]])
        lab = np.concatenate([np.ones(tr.sum()), -np.ones(tr.sum())]).astype(np.int8)
        w = x[:tr.sum()].mean(0) - x[tr.sum():].mean(0)
        for _ in range(iters):
            qv, _ = qvalues_and_peps(x @ w, lab)
            pos = (lab > 0) & (qv < 0.01)
            if pos.sum() < 20:
                break
            neg = lab < 0
            cov = np.cov(np.vstack([x[pos], x[neg]]).T) + 1e-6 * np.eye(p)
            w = np.linalg.solve(cov, x[pos].mean(0) - x[neg].mean(0))
        yt[te] = xt[te] @ w
        yd[te] = xd[te] @ w
    win = yt >= yd
    return (np.where(win, yt, yd), np.where(win, 1, -1).astype(np.int8), win & null)


def run(rng, iters, delta, reps=4, n=120000, pi0=0.5, mu=3.0):
    peps, fl = [], []
    for _ in range(reps):
        s, lab, false_t = simulate(rng, n, pi0, mu, iters, delta)
        _, pep = qvalues_and_peps(s, lab)
        t = lab > 0
        peps.append(pep[t]); fl.append(false_t[t].astype(float))
    pep = np.concatenate(peps); f = np.concatenate(fl)
    idx = np.digitize(pep, BINS) - 1
    rows, num = [], 0.0
    for b in range(len(BINS) - 1):
        m = idx == b
        if not m.any():
            continue
        mp, ob = float(pep[m].mean()), float(f[m].mean())
        rows.append({"lo": float(BINS[b]), "n": int(m.sum()), "mean_pep": mp,
                     "obs": ob, "gap": ob - mp,
                     "ratio": ob / mp if mp > 0 else None})
        num += m.sum() * (ob - mp)
    low = pep < 1e-2
    return {"iters": iters, "delta": delta, "reps": reps, "n_targets": int(pep.size),
            "weighted_signed": num / pep.size,
            "low_n": int(low.sum()),
            "low_mean_pep": float(pep[low].mean()) if low.any() else None,
            "low_obs": float(f[low].mean()) if low.any() else None,
            "low_ratio": (float(f[low].mean()) / float(pep[low].mean()))
                         if low.any() and pep[low].mean() > 0 else None,
            "bins": rows}


def main():
    rng = np.random.default_rng(555)
    out = []
    for delta in (0.0, 0.15, 0.3):
        for iters in (0, 1, 2, 3, 5, 10):
            out.append(run(rng, iters, delta))
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()


# --------------------------------------------------------------------------
# Second arm: homolog-shaped incorrect targets.
#
# A false match to a homolog of an abundant true peptide is not null on the
# fragment-match features -- IAPEEHPVLLTEAPLNPK differs from its native actin
# counterpart at one of eighteen residues, so almost every fragment ion still
# matches.  It is also a real tryptic peptide, so it is not null on the
# sequence-plausibility features either.  A reversed decoy is null on both.
# Nothing here is fitted to the observed calibration.
# --------------------------------------------------------------------------

def simulate_homolog(rng, n, pi_correct, pi_homolog, iters, a=3.0, b=1.0, rho=0.55,
                     n_match=4, n_seq=2, n_noise=10, folds=3):
    u = rng.random(n)
    cls = np.where(u < pi_correct, 2, np.where(u < pi_correct + pi_homolog, 1, 0))
    p = n_match + n_seq + n_noise
    xt = rng.normal(0, 1, (n, p))
    xd = rng.normal(0, 1, (n, p))
    m = slice(0, n_match)
    s = slice(n_match, n_match + n_seq)
    xt[cls == 2, m] += a / np.sqrt(n_match)          # correct: full match signal
    xt[cls == 1, m] += rho * a / np.sqrt(n_match)    # homolog: partial match signal
    xt[:, s] += b / np.sqrt(n_seq)                   # every target is a real peptide
    #  decoys are reversed sequences: null on match AND on sequence plausibility
    yt = np.zeros(n); yd = np.zeros(n)
    fold = rng.integers(0, folds, n)
    for f in range(folds):
        tr, te = fold != f, fold == f
        x = np.vstack([xt[tr], xd[tr]])
        lab = np.concatenate([np.ones(tr.sum()), -np.ones(tr.sum())]).astype(np.int8)
        w = np.zeros(p); w[0] = 1.0                  # start from one match feature
        for _ in range(iters):
            qv, _ = qvalues_and_peps(x @ w, lab)
            pos = (lab > 0) & (qv < 0.01)
            if pos.sum() < 20:
                break
            neg = lab < 0
            cov = np.cov(np.vstack([x[pos], x[neg]]).T) + 1e-6 * np.eye(p)
            w = np.linalg.solve(cov, x[pos].mean(0) - x[neg].mean(0))
        yt[te] = xt[te] @ w
        yd[te] = xd[te] @ w
    win = yt >= yd
    return (np.where(win, yt, yd), np.where(win, 1, -1).astype(np.int8),
            win & (cls != 2), win & (cls == 1))


def run_homolog(rng, iters, pi_homolog, reps=3, n=120000, pi_correct=0.45):
    peps, fl, hm, lab_all = [], [], [], []
    for _ in range(reps):
        s, lab, false_t, homo = simulate_homolog(rng, n, pi_correct, pi_homolog, iters)
        _, pep = qvalues_and_peps(s, lab)
        t = lab > 0
        peps.append(pep[t]); fl.append(false_t[t].astype(float)); hm.append(homo[t])
    pep = np.concatenate(peps); f = np.concatenate(fl)
    idx = np.digitize(pep, BINS) - 1
    rows, num = [], 0.0
    for b in range(len(BINS) - 1):
        mm = idx == b
        if not mm.any():
            continue
        mp, ob = float(pep[mm].mean()), float(f[mm].mean())
        rows.append({"lo": float(BINS[b]), "n": int(mm.sum()), "mean_pep": mp,
                     "obs": ob, "gap": ob - mp, "ratio": ob / mp if mp > 0 else None})
        num += mm.sum() * (ob - mp)
    low = pep < 1e-2
    return {"arm": "homolog", "iters": iters, "pi_homolog": pi_homolog,
            "n_targets": int(pep.size), "weighted_signed": num / pep.size,
            "low_n": int(low.sum()),
            "low_mean_pep": float(pep[low].mean()) if low.any() else None,
            "low_obs": float(f[low].mean()) if low.any() else None,
            "low_ratio": (float(f[low].mean()) / float(pep[low].mean()))
                         if low.any() and pep[low].mean() > 0 else None,
            "bins": rows}
