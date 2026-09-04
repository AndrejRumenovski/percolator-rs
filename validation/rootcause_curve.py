#!/usr/bin/env python3
"""Shared curve utilities: entrapment decomposition at q thresholds and at
matched numbers of accepted targets."""
import csv
from pathlib import Path
csv.field_size_limit(10**9)
THRESHOLDS = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)

def classify(mem, decoy):
    if decoy: mem=[m.removeprefix("DECOY_").removeprefix("decoy_") for m in mem]
    pure = bool(mem) and all(m.startswith("ENT_") for m in mem)
    mixed = any(m.startswith("ENT_") for m in mem) and not pure
    return pure, mixed

def load_tsv(path, decoy, qcol="q-value", scol="score", pcol="proteinIds"):
    out=[]
    with open(path, newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            vals=[]
            for k,v in row.items():
                if k==pcol or k is None:
                    if isinstance(v,list): vals.extend(v)
                    elif v: vals.append(v)
            mem=[p for v in vals for p in v.replace(";","\t").split("\t") if p]
            pure,mixed = classify(mem, decoy)
            out.append((float(row[scol]), float(row[qcol]), 1 if not decoy else -1, pure, mixed))
    return out

def q_table(rows, thresholds=THRESHOLDS):
    res=[]
    for t in thresholds:
        sel=[r for r in rows if r[1] < t]
        tg=[r for r in sel if r[2]==1]; dc=[r for r in sel if r[2]==-1]
        Rent=sum(r[3] for r in tg)
        Dent=sum(r[3] for r in dc); Dnat=sum(1 for r in dc if not r[3] and not r[4])
        f=Dent/(Dent+Dnat) if Dent+Dnat else float("nan")
        res.append(dict(t=t,R=len(tg),Rent=Rent,D=len(dc),Dent=Dent,Dnat=Dnat,f=f,
                        ratio=Rent/Dent if Dent else float("nan"),
                        rawfdp=Rent/len(tg) if tg else float("nan"),
                        adj=(Rent/f)/len(tg) if tg and f else float("nan")))
    return res

def r_table(rows, targets_wanted):
    """Walk the score-sorted list; report entrapment stats when R hits each target count."""
    rows=sorted(rows, key=lambda x:-x[0])
    R=D=Rent=Dent=Dnat=0
    want=sorted(targets_wanted); i=0; out=[]
    for s,q,lab,pure,mixed in rows:
        if lab==1:
            R+=1; Rent+=pure
        else:
            D+=1; Dent+=pure; Dnat+= (not pure and not mixed)
        while i<len(want) and R>=want[i]:
            f=Dent/(Dent+Dnat) if Dent+Dnat else float("nan")
            out.append(dict(R=R,D=D,Rent=Rent,Dent=Dent,Dnat=Dnat,f=f,
                            tdc_fdp=D/R if R else float("nan"),
                            ratio=Rent/Dent if Dent else float("nan"),
                            rawfdp=Rent/R if R else float("nan"),
                            adj=(Rent/f)/R if R and f else float("nan")))
            i+=1
        if i>=len(want): break
    return out

def show_q(name, rows):
    print(f"--- {name} (q thresholds) ---")
    print(f"{'t':<8}{'R':>8}{'R_ent':>7}{'D':>7}{'D_ent':>7}{'f':>8}{'Rent/Dent':>11}{'rawEntFDP':>11}{'adjFDP':>9}{'adj/nom':>9}")
    for r in q_table(rows):
        print(f"{r['t']:<8}{r['R']:>8}{r['Rent']:>7}{r['D']:>7}{r['Dent']:>7}{r['f']:>8.4f}"
              f"{r['ratio']:>11.3f}{r['rawfdp']:>11.5f}{r['adj']:>9.5f}{r['adj']/r['t']:>9.2f}")

def show_r(name, rows, wants):
    print(f"--- {name} (matched accepted-target counts) ---")
    print(f"{'R':>8}{'D':>8}{'TDC D/R':>9}{'R_ent':>7}{'D_ent':>7}{'f':>8}{'Rent/Dent':>11}{'rawEntFDP':>11}{'adjFDP':>9}")
    for r in r_table(rows, wants):
        print(f"{r['R']:>8}{r['D']:>8}{r['tdc_fdp']:>9.4f}{r['Rent']:>7}{r['Dent']:>7}{r['f']:>8.4f}"
              f"{r['ratio']:>11.3f}{r['rawfdp']:>11.5f}{r['adj']:>9.5f}")
