#!/usr/bin/env python3
"""Significance and uncertainty for the adjustment-free internal null.

Statistic: among entrapment PSMs falling in a PEP region, the share that are
targets rather than decoys.  Both populations are certainly false, so under
target-decoy exchangeability that share is 1/2 whatever the database sizes,
opportunity ratios or entrapment fraction happen to be.  Reported at PSM level
and collapsed to distinct peptides, which removes the clustering that PSM-level
counting hides.
"""
import json
import math
import sys

import numpy as np

REGIONS = [("PEP < 1e-3", 0.0, 1e-3), ("PEP < 0.01", 0.0, 1e-2),
           ("PEP < 0.05", 0.0, 5e-2), ("PEP < 0.5", 0.0, 0.5),
           ("PEP >= 0.5", 0.5, 1.0 + 1e-12)]


def z_binom(k, n):
    if n == 0:
        return None
    return (k - n / 2.0) / math.sqrt(n * 0.25)


def main():
    d = dict(np.load(sys.argv[1], allow_pickle=True))
    pep, pure, dec, pid = d["pep"], d["pure"], d["decoy"] == 1, d["pepid"]
    seed, ds = d["seed"], d["ds"]
    out = {"psm_level": [], "peptide_level": [], "by_seed": [], "bootstrap": {}}

    for name, lo, hi in REGIONS:
        m = pure & (pep >= lo) & (pep < hi)
        t = int((m & ~dec).sum()); k = int((m & dec).sum())
        out["psm_level"].append({"region": name, "ent_targets": t, "ent_decoys": k,
                                 "ratio": (t / k) if k else None,
                                 "z": z_binom(t, t + k)})
        # collapse to distinct peptide sequences within each seed, then pool
        tp = len(np.unique(pid[m & ~dec])); dp = len(np.unique(pid[m & dec]))
        out["peptide_level"].append({"region": name, "ent_target_peptides": tp,
                                     "ent_decoy_peptides": dp,
                                     "ratio": (tp / dp) if dp else None,
                                     "z": z_binom(tp, tp + dp)})

    for s in np.unique(seed):
        row = {"seed": int(s)}
        for name, lo, hi in REGIONS[:3]:
            m = pure & (pep >= lo) & (pep < hi) & (seed == s)
            t = int((m & ~dec).sum()); k = int((m & dec).sum())
            row[name] = {"t": t, "d": k, "ratio": (t / k) if k else None,
                         "z": z_binom(t, t + k)}
        out["by_seed"].append(row)

    # cluster bootstrap of the PEP<0.01 target share, resampling whole runs and
    # whole peptide sequences
    rng = np.random.default_rng(20260827)
    m = pure & (pep < 1e-2)
    for unit_name, unit in (("run", ds[m]), ("peptide", pid[m])):
        istgt = (~dec)[m]
        codes, inv = np.unique(unit, return_inverse=True)
        groups = [np.nonzero(inv == i)[0] for i in range(codes.size)]
        vals = []
        for _ in range(2000):
            pick = rng.integers(0, len(groups), len(groups))
            idx = np.concatenate([groups[i] for i in pick])
            if idx.size == 0:
                continue
            share = istgt[idx].mean()
            vals.append(share / (1 - share) if share < 1 else float("inf"))
        v = np.array([x for x in vals if np.isfinite(x)])
        out["bootstrap"][unit_name] = {
            "n_units": int(codes.size),
            "point_ratio": float(istgt.mean() / (1 - istgt.mean())),
            "ci95": [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))],
        }
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
