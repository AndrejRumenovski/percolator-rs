#!/usr/bin/env python3
"""Load the entrapment PSM outputs into one cached array set.

Columns kept per row: dataset index, seed, decoy flag, score, q, PEP, a
pure-entrapment flag, a mixed flag, and a peptide key for clustered resampling.
"""
import csv
import sys
from pathlib import Path

import numpy as np

csv.field_size_limit(10 ** 9)


def strip_flanks(p):
    return p[2:-2] if len(p) > 4 and p[1] == "." and p[-2] == "." else p


def load_file(path, decoy):
    scores, qs, peps, pure, mixed, peptides = [], [], [], [], [], []
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
            if decoy:
                mem = [m.removeprefix("DECOY_").removeprefix("decoy_") for m in mem]
            p = bool(mem) and all(m.startswith("ENT_") for m in mem)
            scores.append(float(row["score"]))
            qs.append(float(row["q-value"]))
            pv = row["posterior_error_prob"]
            peps.append(float("nan") if pv == "NA" else float(pv))
            pure.append(p)
            mixed.append(any(m.startswith("ENT_") for m in mem) and not p)
            peptides.append(strip_flanks(row["peptide"]))
    return scores, qs, peps, pure, mixed, peptides


def build(root, seeds, cache):
    if cache.exists():
        return dict(np.load(cache, allow_pickle=True))
    datasets = sorted(d.name for d in (root / f"seed-{seeds[0]}").iterdir() if d.is_dir())
    cols = {k: [] for k in
            ("ds", "seed", "decoy", "score", "q", "pep", "pure", "mixed", "pepid")}
    keys = {}
    for seed in seeds:
        for di, ds in enumerate(datasets):
            for decoy, fn in ((0, "target.tsv"), (1, "decoy.tsv")):
                s, q, p, pu, mx, pe = load_file(root / f"seed-{seed}" / ds / fn, decoy)
                n = len(s)
                cols["ds"].append(np.full(n, di, np.int16))
                cols["seed"].append(np.full(n, seed, np.int16))
                cols["decoy"].append(np.full(n, decoy, np.int8))
                cols["score"].append(np.asarray(s)); cols["q"].append(np.asarray(q))
                cols["pep"].append(np.asarray(p))
                cols["pure"].append(np.asarray(pu, bool)); cols["mixed"].append(np.asarray(mx, bool))
                ids = np.empty(n, np.int64)
                for i, key in enumerate(pe):
                    k = (decoy, key)
                    ids[i] = keys.setdefault(k, len(keys))
                cols["pepid"].append(ids)
            print(f"loaded seed {seed} {ds}", file=sys.stderr)
    out = {k: np.concatenate(v) for k, v in cols.items()}
    out["datasets"] = np.array(datasets)
    np.savez_compressed(cache, **out)
    return out


if __name__ == "__main__":
    root = Path(sys.argv[1])
    cache = Path(sys.argv[2])
    d = build(root, [1, 2, 3, 4, 5], cache)
    print({k: (v.shape, v.dtype) for k, v in d.items()})
