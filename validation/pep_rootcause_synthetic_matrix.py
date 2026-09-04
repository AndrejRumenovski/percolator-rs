#!/usr/bin/env python3
"""Compact assumption-holding matrix covering calibration edge cases."""

import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from pep_rootcause_lib import qvalues_and_peps

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2,
                 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])


def normal_density(x, mean):
    return np.exp(-0.5 * (x - mean) ** 2)


def normal_mass(x, mean, grid):
    hi = (x + grid / 2 - mean) / math.sqrt(2)
    lo = (x - grid / 2 - mean) / math.sqrt(2)
    erf = np.vectorize(math.erf)
    return 0.5 * (erf(hi) - erf(lo))


def one(rng, n, pi0, separation, grid=None):
    null = rng.random(n) < pi0
    target = np.where(null, rng.normal(0, 1, n), rng.normal(separation, 1, n))
    decoy = rng.normal(0, 1, n)
    if grid is not None:
        target = np.round(target / grid) * grid
        decoy = np.round(decoy / grid) * grid
    tie = target == decoy
    target_wins = (target > decoy) | (tie & (rng.random(n) < 0.5))
    score = np.where(target_wins, target, decoy)
    label = np.where(target_wins, 1, -1).astype(np.int8)
    false = target_wins & null
    _, pep = qvalues_and_peps(score, label)
    selected = label > 0
    score_t = score[selected]
    pep_t = pep[selected]
    false_t = false[selected]
    if pi0 == 1:
        truth = np.ones(score_t.size)
    else:
        f0 = normal_density(score_t, 0) if grid is None else normal_mass(score_t, 0, grid)
        f1 = normal_density(score_t, separation) if grid is None else normal_mass(score_t, separation, grid)
        truth = pi0 * f0 / np.maximum(pi0 * f0 + (1 - pi0) * f1, 1e-300)
    return pep_t, false_t, truth


def aggregate(name, parts, **parameters):
    pep = np.concatenate([x[0] for x in parts])
    false = np.concatenate([x[1] for x in parts])
    truth = np.concatenate([x[2] for x in parts])
    index = np.digitize(pep, BINS) - 1
    bins = []
    for bi in range(len(BINS) - 1):
        mask = index == bi
        if mask.any():
            bins.append({
                "lo": float(BINS[bi]), "hi": float(BINS[bi + 1]), "n": int(mask.sum()),
                "mean_predicted": float(pep[mask].mean()),
                "observed_false": float(false[mask].mean()),
                "mean_true_posterior": float(truth[mask].mean()),
            })
    return {
        "case": name, **parameters, "targets": int(pep.size),
        "mean_predicted": float(pep.mean()), "observed_false": float(false.mean()),
        "mean_true_posterior": float(truth.mean()),
        "observed_minus_predicted": float(false.mean() - pep.mean()),
        "true_posterior_minus_predicted": float(truth.mean() - pep.mean()),
        "bins": bins,
    }


def main():
    rng = np.random.default_rng(20260830)
    specs = [
        ("complete_null", 200_000, 1.0, 0.0, None, 4),
        ("strong_separation", 200_000, 0.5, 4.0, None, 4),
        ("moderate_overlap", 200_000, 0.5, 2.0, None, 4),
        ("heavy_overlap", 200_000, 0.5, 0.5, None, 4),
        ("target_heavy_class_imbalance", 200_000, 0.1, 3.0, None, 4),
        ("near_balanced_class_imbalance", 200_000, 0.9, 3.0, None, 4),
        ("exact_repeated_scores", 200_000, 0.5, 2.0, 0.5, 4),
        ("dense_score_ties", 200_000, 0.5, 2.0, 0.1, 4),
        ("small_sample", 50, 0.5, 3.0, None, 1000),
        ("sparse_extreme_tail", 2_000, 0.5, 4.0, None, 200),
        ("large_sample", 1_000_000, 0.5, 3.0, None, 1),
        ("dense_extreme_tail", 1_000_000, 0.5, 4.0, None, 1),
    ]
    out = []
    for name, n, pi0, separation, grid, reps in specs:
        parts = [one(rng, n, pi0, separation, grid) for _ in range(reps)]
        out.append(aggregate(name, parts, n_per_replicate=n, replicates=reps,
                             pi0=pi0, separation=separation, score_grid=grid))

    # Every score identical and labels exactly balanced: target PEP is exactly 1,
    # matching the complete-null truth.  This also exercises one enormous tie.
    n = 200_000
    labels = np.tile(np.array([1, -1], np.int8), n // 2)
    _, pep = qvalues_and_peps(np.zeros(n), labels)
    t = labels > 0
    out.append({"case": "all_scores_tied_complete_null", "rows": n,
                "targets": int(t.sum()), "mean_predicted": float(pep[t].mean()),
                "observed_false": 1.0,
                "observed_minus_predicted": float(1.0 - pep[t].mean())})
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
