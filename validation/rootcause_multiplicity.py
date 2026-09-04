#!/usr/bin/env python3
"""Does within-spectrum re-competition over an unbalanced candidate pool drive it?

For each spectrum the PIN offers n_T target and n_D decoy candidates.  Percolator
re-ranks with the learned score and keeps the max.  Under exchangeability of
null candidates the target should win with probability n_T/(n_T+n_D), not 1/2.
Stratify the accepted set by pool composition and see whether the entrapment
excess tracks it.
"""
import csv
from pathlib import Path
from collections import defaultdict, Counter

csv.field_size_limit(10**9)
PIN = Path("/home/andrej-rumenovski/percolator_rs_out/entrapment")
OUT = Path("/run/media/andrej-rumenovski/New Volume/percolator_rs_final_audit_20260827/empirical-current/entrapment/seed-1")
DATASETS = sorted(p.name for p in OUT.iterdir() if p.is_dir())

# pool composition per spectrum key, and Comet's own rank-1 label
pool = {}
rank1 = {}
for ds in DATASETS:
    with open(PIN/ds/"comet.pin") as h:
        r = csv.reader(h, delimiter="\t"); hdr = next(r)
        iL,iS,iX = hdr.index("Label"), hdr.index("SpecId"), hdr.index("Xcorr")
        best = {}
        for row in r:
            key = row[iS].rsplit("_",1)[0]
            t,d = pool.get(key,(0,0))
            pool[key] = (t+(row[iL]=="1"), d+(row[iL]=="-1"))
            x = float(row[iX])
            if key not in best or x > best[key][0]: best[key] = (x, row[iL])
        for k,(x,l) in best.items(): rank1[k] = (x, l)

def load(path, decoy):
    out=[]
    with open(path, newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            vals=[]
            for k,v in row.items():
                if k=="proteinIds" or k is None:
                    if isinstance(v,list): vals.extend(v)
                    elif v: vals.append(v)
            mem=[p for v in vals for p in v.replace(";","\t").split("\t") if p]
            if decoy: mem=[m.removeprefix("DECOY_") for m in mem]
            pure = bool(mem) and all(m.startswith("ENT_") for m in mem)
            mixed = any(m.startswith("ENT_") for m in mem) and not pure
            key = row["PSMId"].rsplit("_",1)[0]
            out.append((key, float(row["q-value"]), pure, mixed, decoy))
    return out

rows=[]
for ds in DATASETS:
    rows += load(OUT/ds/"target.tsv", False)
    rows += load(OUT/ds/"decoy.tsv",  True)
print("output PSMs:", len(rows), " with pool info:", sum(1 for r in rows if r[0] in pool))

# --- Global check: does the target-win rate track the pool composition? ---
print("\n=== A. Win rate vs pool composition (ALL spectra, percolator winners) ===")
by = defaultdict(lambda:[0,0])
for key,q,pure,mixed,dec in rows:
    p = pool.get(key)
    if not p: continue
    by[p][1 if dec else 0]+= 1
print(f"{'(nT,nD)':<10}{'spectra':>9}{'T wins':>9}{'D wins':>9}{'P(T win)':>10}{'nT/(nT+nD)':>12}")
tot_obs=tot_exp=0
for k in sorted(by, key=lambda k:-sum(by[k])):
    t,d = by[k]; n=t+d
    if n < 50: continue
    exp = k[0]/(k[0]+k[1]) if sum(k) else float('nan')
    print(f"{str(k):<10}{n:>9}{t:>9}{d:>9}{t/n:>10.4f}{exp:>12.4f}")

# --- B. Restrict to spectra Comet's own rank-1 called a decoy (null-enriched) ---
print("\n=== B. Null-enriched stratum: spectra whose Comet rank-1 was a DECOY ===")
by2 = defaultdict(lambda:[0,0])
for key,q,pure,mixed,dec in rows:
    p = pool.get(key)
    if not p or rank1.get(key,(0,"1"))[1] != "-1": continue
    by2[p][1 if dec else 0]+=1
print(f"{'(nT,nD)':<10}{'spectra':>9}{'T wins':>9}{'D wins':>9}{'P(T win)':>10}{'nT/(nT+nD)':>12}")
for k in sorted(by2, key=lambda k:-sum(by2[k])):
    t,d = by2[k]; n=t+d
    if n < 50: continue
    print(f"{str(k):<10}{n:>9}{t:>9}{d:>9}{t/n:>10.4f}{k[0]/(k[0]+k[1]):>12.4f}")
T=sum(v[0] for v in by2.values()); D=sum(v[1] for v in by2.values())
EXP=sum((v[0]+v[1])*k[0]/(k[0]+k[1]) for k,v in by2.items() if sum(k))
print(f"pooled: T={T} D={D} P(T win)={T/(T+D):.4f}  pool-predicted={EXP/(T+D):.4f}  (0.5 if pools were balanced)")

# --- C. R_ent/D_ent at q<0.01 stratified by whether the pool favoured targets ---
print("\n=== C. Accepted-set entrapment ratio, stratified by pool composition ===")
for label, pred in (("nT > nD", lambda p: p[0] > p[1]),
                    ("nT = nD", lambda p: p[0] == p[1]),
                    ("nT < nD", lambda p: p[0] < p[1])):
    for t in (0.01, 0.05):
        R=Rent=D=Dent=Dnat=0
        for key,q,pure,mixed,dec in rows:
            p = pool.get(key)
            if not p or not pred(p) or q >= t: continue
            if dec:
                D+=1
                if pure: Dent+=1
                elif not mixed: Dnat+=1
            else:
                R+=1
                if pure: Rent+=1
        f = Dent/(Dent+Dnat) if Dent+Dnat else float('nan')
        adj = (Rent/f)/R if R and f else float('nan')
        print(f"{label:<8} q<{t:<6} R={R:>6} R_ent={Rent:>4} D={D:>5} D_ent={Dent:>4} "
              f"f={f:.4f} Rent/Dent={(Rent/Dent if Dent else float('nan')):>6.3f} adjFDP={adj:.5f}")
