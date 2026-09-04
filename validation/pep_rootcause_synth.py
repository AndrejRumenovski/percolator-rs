#!/usr/bin/env python3
"""Synthetic PEP calibration where every TDC assumption holds by construction.

Generative model, one row per spectrum, exactly the contract in src/stats.rs:

  * a spectrum is null with probability pi0, otherwise it carries a true peptide;
  * a null spectrum draws target and decoy scores i.i.d. from f0 -- so incorrect
    targets and decoys are exchangeable by construction and an incorrect target
    beats its decoy with probability exactly 1/2;
  * a correct spectrum draws its target from f1 and its decoy from f0;
  * the reported row is the winner of that competition;
  * spectra are independent.

Ground truth is therefore known per row, with no entrapment extrapolation, and
the true local error probability among reported target winners has the closed
form  lfdr(s) = pi0 f0(s) / (pi0 f0(s) + (1-pi0) f1(s))  -- the competition term
F0(s) cancels between numerator and denominator.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps, pava_non_decreasing

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])


def draw(rng, n, pi0, mu, family):
    """Return (score, label, is_false_target, true_lfdr) for n spectra."""
    null = rng.random(n) < pi0
    if family == "normal":
        f0 = lambda k: rng.normal(0.0, 1.0, k)
        f1 = lambda k: rng.normal(mu, 1.0, k)
        dens0 = lambda s: np.exp(-0.5 * s * s)
        dens1 = lambda s: np.exp(-0.5 * (s - mu) ** 2)
    elif family == "gumbel":                      # search-score-like right tail
        f0 = lambda k: rng.gumbel(0.0, 1.0, k)
        f1 = lambda k: rng.gumbel(mu, 1.0, k)
        dens0 = lambda s: np.exp(-(s + np.exp(-s)))
        dens1 = lambda s: np.exp(-((s - mu) + np.exp(-(s - mu))))
    elif family == "scaled":                      # alternative is also wider
        f0 = lambda k: rng.normal(0.0, 1.0, k)
        f1 = lambda k: rng.normal(mu, 1.6, k)
        dens0 = lambda s: np.exp(-0.5 * s * s)
        dens1 = lambda s: np.exp(-0.5 * ((s - mu) / 1.6) ** 2) / 1.6
    else:
        raise ValueError(family)

    t = np.where(null, f0(n), f1(n))
    d = f0(n)
    target_wins = t >= d
    score = np.where(target_wins, t, d)
    label = np.where(target_wins, 1, -1).astype(np.int8)
    is_false_target = target_wins & null
    with np.errstate(under="ignore"):
        a = pi0 * dens0(score)
        b = (1.0 - pi0) * dens1(score)
    true_lfdr = np.where(a + b > 0, a / np.maximum(a + b, 1e-300), 1.0)
    return score, label, is_false_target, true_lfdr


def bin_table(pep, false_flag, truth=None):
    idx = np.digitize(pep, BINS) - 1
    rows = []
    for b in range(len(BINS) - 1):
        m = idx == b
        k = int(m.sum())
        if k == 0:
            continue
        row = {
            "lo": float(BINS[b]), "hi": float(BINS[b + 1]), "n": k,
            "mean_pep": float(pep[m].mean()),
            "obs_false": float(false_flag[m].mean()),
        }
        row["gap"] = row["obs_false"] - row["mean_pep"]
        if truth is not None:
            row["mean_true_lfdr"] = float(truth[m].mean())
        rows.append(row)
    return rows


def weighted_errors(pep, false_flag, truth=None):
    tab = bin_table(pep, false_flag, truth)
    n = sum(r["n"] for r in tab)
    signed = sum(r["n"] * r["gap"] for r in tab) / n
    absolute = sum(r["n"] * abs(r["gap"]) for r in tab) / n
    out = {"weighted_signed": signed, "weighted_abs": absolute,
           "sum_pep": float(pep.sum()), "n_false": int(false_flag.sum())}
    if truth is not None:
        out["sum_true_lfdr"] = float(truth.sum())
    return out, tab


def one_run(rng, n, pi0, mu, family):
    score, label, false_target, true_lfdr = draw(rng, n, pi0, mu, family)
    _, pep = qvalues_and_peps(score, label)
    t = label > 0
    return {
        "pep": pep[t], "false": false_target[t].astype(float),
        "truth": true_lfdr[t], "score": score[t],
    }


def main():
    out = {}
    rng = np.random.default_rng(31337)

    # --- primary grid ---------------------------------------------------
    grid = []
    for family in ("normal", "gumbel", "scaled"):
        for pi0 in (0.3, 0.5, 0.7, 0.9):
            for mu in (2.0, 3.0, 4.0):
                for n in (20000, 200000):
                    r = one_run(rng, n, pi0, mu, family)
                    err, tab = weighted_errors(r["pep"], r["false"], r["truth"])
                    grid.append({
                        "family": family, "pi0": pi0, "mu": mu, "n": n,
                        "n_targets": int(r["pep"].size),
                        **err,
                        "sum_pep_over_n_false": (float(r["pep"].sum()) /
                                                 max(int(r["false"].sum()), 1)),
                        "bins": tab,
                    })
    out["grid"] = grid

    # --- replicate study at one setting, for sampling error -------------
    reps = []
    for _ in range(40):
        r = one_run(rng, 200000, 0.5, 3.0, "normal")
        err, tab = weighted_errors(r["pep"], r["false"], r["truth"])
        low = r["pep"] < 1e-3
        reps.append({
            **err,
            "n_low": int(low.sum()),
            "mean_pep_low": float(r["pep"][low].mean()) if low.any() else None,
            "obs_false_low": float(r["false"][low].mean()) if low.any() else None,
            "true_lfdr_low": float(r["truth"][low].mean()) if low.any() else None,
        })
    out["replicates_normal_pi0.5_mu3_n200k"] = reps

    # --- size scaling of the low-PEP bin --------------------------------
    scaling = []
    for n in (2000, 10000, 50000, 200000, 1000000):
        acc = []
        for _ in range(8):
            r = one_run(rng, n, 0.5, 3.0, "normal")
            low = r["pep"] < 1e-3
            if low.sum() > 0:
                acc.append((int(low.sum()), float(r["pep"][low].mean()),
                            float(r["false"][low].mean()),
                            float(r["truth"][low].mean())))
        if acc:
            a = np.array(acc)
            scaling.append({"n": n, "reps": len(acc),
                            "mean_n_low": float(a[:, 0].mean()),
                            "mean_pep_low": float(a[:, 1].mean()),
                            "mean_obs_false_low": float(a[:, 2].mean()),
                            "mean_true_lfdr_low": float(a[:, 3].mean())})
    out["low_bin_size_scaling"] = scaling
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
