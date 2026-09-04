import csv, math, statistics
import numpy as np
from pathlib import Path
csv.field_size_limit(10**9)
RS=Path("/run/media/andrej-rumenovski/New Volume/percolator_rs_final_audit_20260827/empirical-current/entrapment")
DS=[d.name for d in sorted((RS/"seed-1").iterdir()) if d.is_dir()]
def strip(p): return p[2:-2] if len(p)>4 and p[1]=="." and p[-2]=="." else p
def load(seed):
    ds_i=[];dec=[];q=[];pu=[];mx=[];pep=[]
    for i,ds in enumerate(DS):
        for decoy,fn in ((0,"target.tsv"),(1,"decoy.tsv")):
            with open(RS/f"seed-{seed}"/ds/fn,newline="") as h:
                for row in csv.DictReader(h,delimiter="\t"):
                    vals=[]
                    for k,v in row.items():
                        if k=="proteinIds" or k is None:
                            if isinstance(v,list): vals.extend(v)
                            elif v: vals.append(v)
                    mem=[x for v in vals for x in v.replace(";","\t").split("\t") if x]
                    if decoy: mem=[m.removeprefix("DECOY_") for m in mem]
                    p=bool(mem) and all(m.startswith("ENT_") for m in mem)
                    ds_i.append(i);dec.append(decoy);q.append(float(row["q-value"]))
                    pu.append(p);mx.append(any(m.startswith("ENT_") for m in mem) and not p)
                    pep.append((decoy,strip(row["peptide"])))
    return (np.array(ds_i),np.array(dec,bool),np.array(q),np.array(pu,bool),np.array(mx,bool),pep)
dsi,dec,q,pu,mx,pep = load(1)
sel = q<0.01
R=int((sel&~dec).sum()); Rent=int((sel&~dec&pu).sum())
De=int((sel&dec&pu).sum()); Dn=int((sel&dec&~pu&~mx).sum())
fbulk=int((dec&pu).sum())/(int((dec&pu).sum())+int((dec&~pu&~mx).sum())); ftail=De/(De+Dn)
print("=== SENSITIVITY TO f (seed 1, q<0.01) ===")
print(f"  R={R}  R_ent={Rent}  entrapment decoys={De}  native decoys={Dn}")
for name,f in (("f = 1.0000  no adjustment: strict lower bound on FDP",1.0),
               (f"f = {ftail:.4f}  decoys accepted at q<0.01  <- what the audit uses",ftail),
               (f"f = {fbulk:.4f}  all decoys (no tail selection)",fbulk),
               ("f = 0.5000  the design's stated amino-acid balance",0.5)):
    print(f"    {name:<58} adjusted FDP = {(Rent/f)/R:.5f}")
rng=np.random.default_rng(17)
def boot(unit,B=4000):
    codes,inv=np.unique(unit,return_inverse=True)
    order=np.argsort(inv,kind="stable"); inv_s=inv[order]
    starts=np.searchsorted(inv_s,np.arange(len(codes)))
    s_,d_,p_,m_=sel[order],dec[order],pu[order],mx[order]
    cR=np.add.reduceat((s_&~d_).astype(np.int64),starts)
    cRe=np.add.reduceat((s_&~d_&p_).astype(np.int64),starts)
    cDe=np.add.reduceat((s_&d_&p_).astype(np.int64),starts)
    cDn=np.add.reduceat((s_&d_&~p_&~m_).astype(np.int64),starts)
    n=len(codes); idx=rng.integers(0,n,size=(B,n))
    R_=cR[idx].sum(1);Re_=cRe[idx].sum(1);De_=cDe[idx].sum(1);Dn_=cDn[idx].sum(1)
    ok=(R_>0)&(De_>0)
    v=np.sort((Re_[ok]/(De_[ok]/(De_[ok]+Dn_[ok])))/R_[ok])
    return v.mean(),v[int(.025*len(v))],v[int(.975*len(v))],len(v)
print("\n=== CLUSTER BOOTSTRAP of the adjusted FDP at q<0.01 (seed 1) ===")
for nm,u,B in (("resample the 6 LC-MS/MS runs",dsi,4000),
               ("resample peptide sequences",np.array([hash(p) for p in pep],np.int64),2000),
               ("resample PSMs (ignores clustering)",np.arange(len(q)),1200)):
    m,lo,hi,nb=boot(u,B)
    print(f"  {nm:<36} mean={m:.5f}  95% CI [{lo:.5f}, {hi:.5f}]  ({nb} resamples)")
print("\n=== FULL THRESHOLD CURVE, three readings of the same data (mean of 5 seeds) ===")
print(f"{'nominal q':>10}{'adjFDP(audit f)':>17}{'ratio':>8}{'raw entFDP(f=1)':>18}{'ratio':>8}{'Rent/Dent':>11}{'R':>9}")
acc={}
for seed in (1,2,3,4,5):
    d2,de2,q2,pu2,mx2,_=load(seed)
    for t in (0.001,0.005,0.01,0.02,0.05,0.10):
        s=q2<t
        Rr=int((s&~de2).sum()); Re=int((s&~de2&pu2).sum())
        Dee=int((s&de2&pu2).sum()); Dnn=int((s&de2&~pu2&~mx2).sum())
        f=Dee/(Dee+Dnn) if Dee+Dnn else float('nan')
        acc.setdefault(t,[]).append(((Re/f)/Rr, Re/Rr, Re/Dee if Dee else float('nan'), Rr))
for t in (0.001,0.005,0.01,0.02,0.05,0.10):
    a=acc[t]; m=lambda i: statistics.fmean(x[i] for x in a)
    print(f"{t:>10}{m(0):>17.5f}{m(0)/t:>8.2f}{m(1):>18.5f}{m(1)/t:>8.2f}{m(2):>11.3f}{m(3):>9.0f}")
