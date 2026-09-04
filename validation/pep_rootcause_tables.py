#!/usr/bin/env python3
"""Real-data PEP calibration tables, stratified and with uncertainty."""
import json
import sys

import numpy as np

BINS = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0 + 1e-12])
NAMES = ["[0,1e-4)", "[1e-4,1e-3)", "[1e-3,5e-3)", "[5e-3,.01)", "[.01,.02)",
         "[.02,.05)", "[.05,.10)", "[.10,.20)", "[.20,.50)", "[.50,1]"]


def bin_index(pep):
    return np.digitize(pep, BINS) - 1


def f_global(d, mask=None):
    """Fraction of pure-vs-native decoys that are entrapment decoys."""
    m = d["decoy"] == 1
    if mask is not None:
        m &= mask
    pure = d["pure"][m]
    mixed = d["mixed"][m]
    keep = ~mixed
    return float(pure[keep].sum()) / max(int(keep.sum()), 1)


def calib(d, sel_target, sel_decoy, f_mode="bin"):
    """Per-bin calibration table for the target rows selected by sel_target."""
    t = (d["decoy"] == 0) & sel_target
    dc = (d["decoy"] == 1) & sel_decoy
    bt = bin_index(d["pep"])
    fglob = f_global(d, sel_decoy)
    rows = []
    for b in range(len(BINS) - 1):
        mt = t & (bt == b)
        n = int(mt.sum())
        if n == 0:
            continue
        md = dc & (bt == b)
        dpure = int((d["pure"] & md).sum())
        dnat = int((~d["pure"] & ~d["mixed"] & md).sum())
        fbin = dpure / (dpure + dnat) if (dpure + dnat) > 0 else float("nan")
        f = fbin if (f_mode == "bin" and dpure + dnat >= 20) else fglob
        tpure = int((d["pure"] & mt).sum())
        raw = tpure / n
        rows.append({
            "bin": NAMES[b], "n": n, "mean_pep": float(d["pep"][mt].mean()),
            "n_pure_target": tpure, "n_pure_decoy": dpure, "n_native_decoy": dnat,
            "f_used": f, "f_bin": fbin,
            "obs_raw": raw,                       # entrapment targets / all targets
            "obs_f1": raw,                        # strict lower bound (f = 1)
            "obs_adj": raw / f if f > 0 else float("nan"),
            "ent_target_over_ent_decoy": (tpure / dpure) if dpure > 0 else None,
        })
    for r in rows:
        r["gap_adj"] = r["obs_adj"] - r["mean_pep"]
        r["gap_f1"] = r["obs_f1"] - r["mean_pep"]
    n = sum(r["n"] for r in rows)
    summary = {
        "n_targets": n,
        "weighted_signed_adj": sum(r["n"] * r["gap_adj"] for r in rows) / n,
        "weighted_abs_adj": sum(r["n"] * abs(r["gap_adj"]) for r in rows) / n,
        "weighted_signed_f1": sum(r["n"] * r["gap_f1"] for r in rows) / n,
        "weighted_abs_f1": sum(r["n"] * abs(r["gap_f1"]) for r in rows) / n,
        "f_global": fglob,
    }
    return rows, summary


def main():
    d = dict(np.load(sys.argv[1], allow_pickle=True))
    out = {"datasets": [str(x) for x in d["datasets"]]}
    allsel = np.ones(d["pep"].size, bool)

    # audit-style pooled table, both f conventions and a bin-local f
    for mode, tag in (("global", "pooled_fglobal"), ("bin", "pooled_fbin")):
        rows, summary = calib(d, allsel, allsel, f_mode=mode)
        out[tag] = {"bins": rows, "summary": summary}

    # contribution decomposition, on the audit's own convention
    rows = out["pooled_fglobal"]["bins"]
    n = sum(r["n"] for r in rows)
    out["contribution"] = [{"bin": r["bin"], "weight": r["n"] / n,
                            "gap_adj": r["gap_adj"],
                            "contrib_adj": r["n"] * r["gap_adj"] / n,
                            "gap_f1": r["gap_f1"],
                            "contrib_f1": r["n"] * r["gap_f1"] / n} for r in rows]

    # per seed
    out["by_seed"] = []
    for s in np.unique(d["seed"]):
        m = d["seed"] == s
        rows, summary = calib(d, m, m, f_mode="global")
        out["by_seed"].append({"seed": int(s), "summary": summary, "bins": rows})

    # per dataset (pooled over seeds)
    out["by_dataset"] = []
    for i, name in enumerate(out["datasets"]):
        m = d["ds"] == i
        rows, summary = calib(d, m, m, f_mode="global")
        out["by_dataset"].append({"dataset": name, "summary": summary, "bins": rows})

    # cluster bootstrap of the pooled weighted signed error
    rng = np.random.default_rng(4242)
    def stat(mask):
        _, s = calib(d, mask, mask, f_mode="global")
        return s["weighted_signed_adj"], s["weighted_signed_f1"]
    boots = {}
    for unit_name, unit in (("dataset", d["ds"]), ("seed", d["seed"])):
        codes = np.unique(unit)
        vals = []
        for _ in range(400):
            pick = rng.choice(codes, size=codes.size, replace=True)
            mask = np.zeros(unit.size, bool)
            for c in pick:                      # with replacement, so union-of-draws
                mask |= unit == c
            vals.append(stat(mask))
        v = np.array(vals)
        boots[unit_name] = {
            "adj": [float(np.percentile(v[:, 0], 2.5)), float(np.percentile(v[:, 0], 97.5))],
            "f1": [float(np.percentile(v[:, 1], 2.5)), float(np.percentile(v[:, 1], 97.5))],
        }
    out["bootstrap_ci"] = boots
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
