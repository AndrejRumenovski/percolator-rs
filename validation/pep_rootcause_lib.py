"""Reimplementation of the reported percolator-rs PEP estimator, in numpy.

Validated bit-for-bit against `validation/pep_rootcause_probe.rs --score-file`,
which calls the production `stats::qvalues_and_peps` directly.  Used only so the
simulations below can run at a scale the process-per-replicate route cannot.
"""
import numpy as np


def tie_groups(scores_desc):
    """Boundaries of exact-score tie groups in a descending-sorted array."""
    n = scores_desc.size
    if n == 0:
        return np.empty(0, dtype=np.int64)
    ends = np.nonzero(np.diff(scores_desc) != 0)[0]
    return np.concatenate([ends, [n - 1]])


def qvalues_and_peps(scores, labels, p=0.5, plus_one=True):
    """Return (q, pep) aligned to input order.  labels: +1 target, -1 decoy."""
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int8)
    n = scores.size
    lam = p / (1.0 - p)
    order = np.argsort(-scores, kind="stable")
    s = scores[order]
    is_t = labels[order] > 0
    tcum = np.cumsum(is_t).astype(np.float64)
    dcum = np.cumsum(~is_t).astype(np.float64) + (1.0 if plus_one else 0.0)
    ends = tie_groups(s)

    # ---- q-values: one raw FDP per tie group, then reverse cumulative min ----
    fdp_g = np.clip(lam * dcum[ends] / np.maximum(tcum[ends], 1.0), 0.0, 1.0)
    starts = np.concatenate([[0], ends[:-1] + 1])
    fdp = np.repeat(fdp_g, ends - starts + 1)
    q_sorted = np.minimum.accumulate(fdp[::-1])[::-1]
    q = np.empty(n)
    q[order] = q_sorted

    # ---- PEP: increments of F, shared inside a group, then isotonic ----
    tg = tcum[ends] - np.concatenate([[0.0], tcum[ends[:-1]]])  # targets per group
    keep = tg > 0
    d_at = dcum[ends][keep]
    t_at = tcum[ends][keep]
    ng = tg[keep]
    # F_g = min(lam*D_g, T_g, F_{g-1} + n_g).  The third bound is sequential, so
    # it is applied in a loop over groups that hold targets.
    cap = np.minimum(lam * d_at, t_at)
    f = np.empty(cap.size)
    assigned = 0.0
    for i in range(cap.size):
        assigned = min(cap[i], assigned + ng[i])
        f[i] = assigned
    inc = np.diff(np.concatenate([[0.0], f]))
    share = np.maximum(inc / ng, 0.0)
    values = np.repeat(share, ng.astype(np.int64))
    values = pava_non_decreasing(values)
    values = np.clip(values, 1e-12, 1.0)

    pep = np.ones(n)
    t_rows = order[is_t]
    pep[t_rows] = values
    # Decoys take the value of the nearest target at or above them.
    filled = np.empty(n)
    idx = np.cumsum(is_t) - 1
    idx = np.maximum(idx, 0)
    filled = values[idx]
    d_rows = order[~is_t]
    pep[d_rows] = np.clip(filled[~is_t], 1e-12, 1.0)
    return q, pep


def pava_non_decreasing(y):
    """Unit-weight isotonic regression under a non-decreasing constraint."""
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    if n == 0:
        return y.copy()
    val = np.empty(n)
    wt = np.empty(n)
    ln = np.empty(n, dtype=np.int64)
    top = -1
    for i in range(n):
        top += 1
        val[top] = y[i]
        wt[top] = 1.0
        ln[top] = 1
        while top > 0 and val[top - 1] > val[top]:
            w = wt[top - 1] + wt[top]
            val[top - 1] = (val[top - 1] * wt[top - 1] + val[top] * wt[top]) / w
            wt[top - 1] = w
            ln[top - 1] += ln[top]
            top -= 1
    return np.repeat(val[: top + 1], ln[: top + 1])
