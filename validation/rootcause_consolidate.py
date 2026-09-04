#!/usr/bin/env python3
"""Consolidate the entrapment root-cause investigation into one JSON record."""
import csv, json, math, statistics, hashlib
import numpy as np
from pathlib import Path
csv.field_size_limit(10**9)

SP=Path(__file__).resolve().parent
PIN=Path("/home/andrej-rumenovski/percolator_rs_out/entrapment")
RS=Path("/run/media/andrej-rumenovski/New Volume/percolator_rs_final_audit_20260827/empirical-current/entrapment")
DS=[d.name for d in sorted((RS/"seed-1").iterdir()) if d.is_dir()]
TH=(0.001,0.005,0.01,0.02,0.05,0.10)

def cls(mem,decoy):
    if decoy: mem=[m.removeprefix("DECOY_").removeprefix("decoy_") for m in mem]
    pure=bool(mem) and all(m.startswith("ENT_") for m in mem)
    return pure, (any(m.startswith("ENT_") for m in mem) and not pure)

def tsv(path,decoy):
    out=[]
    with open(path,newline="") as h:
        for row in csv.DictReader(h,delimiter="\t"):
            vals=[]
            for k,v in row.items():
                if k=="proteinIds" or k is None:
                    if isinstance(v,list): vals.extend(v)
                    elif v: vals.append(v)
            mem=[x for v in vals for x in v.replace(";","\t").split("\t") if x]
            p,m=cls(mem,decoy)
            out.append((float(row["score"]),float(row["q-value"]),1 if not decoy else -1,p,m))
    return out

def comet(ds,col="Xcorr",hb=True):
    best={}
    with open(PIN/ds/"comet.pin") as h:
        r=csv.reader(h,delimiter="\t"); hdr=next(r)
        iL,iS,iX,iP=hdr.index("Label"),hdr.index("SpecId"),hdr.index(col),hdr.index("Proteins")
        for row in r:
            k=row[iS].rsplit("_",1)[0]; x=float(row[iX]);  x = x if hb else -x
            if k not in best or x>best[k][0]: best[k]=(x,row)
    rows=[]
    for x,row in best.values():
        p,m=cls([q for q in row[iP:] if q], row[iL]=="-1")
        rows.append((x,int(row[iL]),p,m))
    rows.sort(key=lambda v:-v[0])
    T=D=0; raw=[]
    for _,l,_,_ in rows:
        if l==1: T+=1
        else: D+=1
        raw.append((D+1)/T if T else 1.0)
    qq=[0.0]*len(raw); mm=1.0
    for i in range(len(raw)-1,-1,-1):
        mm=min(mm,raw[i]); qq[i]=mm
    return [(rows[i][0],qq[i],rows[i][1],rows[i][2],rows[i][3]) for i in range(len(rows))]

def qtab(rows):
    o=[]
    for t in TH:
        s=[r for r in rows if r[1]<t]
        tg=[r for r in s if r[2]==1]; dc=[r for r in s if r[2]==-1]
        Re=sum(r[3] for r in tg); De=sum(r[3] for r in dc)
        Dn=sum(1 for r in dc if not r[3] and not r[4])
        f=De/(De+Dn) if De+Dn else None
        o.append(dict(q_threshold=t,accepted_targets=len(tg),pure_entrapment_targets=Re,
                      accepted_decoys=len(dc),pure_entrapment_decoys=De,native_decoys=Dn,
                      effective_entrapment_fraction=f,
                      entrapment_target_over_decoy=(Re/De if De else None),
                      unadjusted_entrapment_fdp=(Re/len(tg) if tg else None),
                      adjusted_fdp=((Re/f)/len(tg) if tg and f else None),
                      binomial_z=((Re-(Re+De)/2)/math.sqrt((Re+De)/4) if Re+De else None)))
    return o

arms={}
arms["comet_xcorr_only"]      =[r for ds in DS for r in comet(ds)]
arms["comet_evalue_only"]     =[r for ds in DS for r in comet(ds,"lnExpect",False)]
arms["percolator_rs_top5_seed1"]=[r for ds in DS for r in (tsv(RS/"seed-1"/ds/"target.tsv",False)+tsv(RS/"seed-1"/ds/"decoy.tsv",True))]
arms["percolator_rs_rank1_seed1"]=[r for ds in DS for r in (tsv(SP/"rank1_out/seed-1"/ds/"target.tsv",False)+tsv(SP/"rank1_out/seed-1"/ds/"decoy.tsv",True))]
arms["cpp_percolator_post_processing_tdc_seed1"]=[r for ds in DS for r in (tsv(SP/"cpp_tdc"/f"{ds}.target.tsv",False)+tsv(SP/"cpp_tdc"/f"{ds}.decoy.tsv",True))]

maxiter={}
for mi in (0,1,2,3,5,10):
    per=[]
    for seed in (1,2,3):
        rows=[r for ds in DS for r in (tsv(SP/f"maxiter/mi-{mi}/seed-{seed}"/ds/"target.tsv",False)+
                                       tsv(SP/f"maxiter/mi-{mi}/seed-{seed}"/ds/"decoy.tsv",True))]
        per.append(qtab(rows)[2])
    maxiter[str(mi)]={k:(statistics.fmean(p[k] for p in per) if isinstance(per[0][k],(int,float)) else per[0][k])
                      for k in per[0]}
    maxiter[str(mi)]["seeds"]=[1,2,3]

perdataset={}
for ds in DS:
    v=[]
    for seed in (1,2,3,4,5):
        rows=tsv(RS/f"seed-{seed}"/ds/"target.tsv",False)+tsv(RS/f"seed-{seed}"/ds/"decoy.tsv",True)
        v.append(qtab(rows)[2])
    usable=[x["adjusted_fdp"] for x in v if x["adjusted_fdp"] is not None]
    perdataset[ds]=dict(seeds_with_defined_adjusted_fdp=len(usable),
                        mean_adjusted_fdp=(statistics.fmean(usable) if usable else None),
                        accepted_targets_seed1=v[0]["accepted_targets"],
                        pure_entrapment_targets_seed1=v[0]["pure_entrapment_targets"],
                        pure_entrapment_decoys_seed1=v[0]["pure_entrapment_decoys"])

out=dict(
  schema_version=1,
  investigation="causal diagnosis of the ~1.81% adjusted FDP at nominal q<0.01 in the signal-present entrapment experiment",
  date="2026-08-27",
  audited_binary_sha256="be9bf670bfd69df4dc3ba3b8be6c4c164acaf56a91f4a2819d115f49604b2c45",
  production_code_modified=False,
  identity="adjusted_fdp == (accepted_decoys_used/accepted_targets) * (pure_entrapment_targets/pure_entrapment_decoys)",
  arms={k:qtab(v) for k,v in arms.items()},
  training_dose_response_at_q_0_01=maxiter,
  per_dataset_at_q_0_01=perdataset,
  database_construction=dict(
    native_proteins=139191, native_amino_acids=145351640,
    entrapment_proteins=389504, entrapment_amino_acids=145351799,
    declared_entrapment_fraction=0.50000027,
    native_distinct_7mers=40633136, entrapment_distinct_7mers=61699335,
    native_distinct_7mers_per_aa=0.2796, entrapment_distinct_7mers_per_aa=0.4245,
    observed_effective_entrapment_fraction_from_decoys=0.7823,
    comment="the design balanced amino acids; searchable sequence content is ~1.5x richer on the entrapment side and the realised opportunity ratio is ~3.5:1"),
  provenance_leakage=dict(
    leave_one_file_out_auc_entrapment_target_vs_entrapment_decoy=0.5209,
    label_permuted_control_auc=0.5024,
    largest_standardised_weights={"enzN":0.1386,"enzC":0.1036,"enzInt":-0.0283},
    percolator_score_bulk_auc_same_contrast=0.5012,
    comment="feature-space provenance signal is real but bulk-score AUC is 0.50: the violation is confined to the extreme upper tail"),
  homology_channel=dict(
    metric="longest substring shared with the native proteome, I/L collapsed, k>=7",
    accepted_entrapment_targets=dict(n=150, mean_lcs_over_len=0.360, frac_lcs_ge_12aa=0.087),
    unaccepted_entrapment_targets=dict(n=6000, mean_lcs_over_len=0.294, frac_lcs_ge_12aa=0.047),
    accepted_entrapment_decoys=dict(n=115, mean_lcs_over_len=0.328, frac_lcs_ge_12aa=0.043),
    top_hits=[{"peptide":"IAPEEHPVLLTEAPLNPK","psms":26,"len":18,"native_shared_substring_aa":17},
              {"peptide":"AVFVDLEPTVIDEVR","psms":22,"len":15,"native_shared_substring_aa":10},
              {"peptide":"VVPEEHPVLLTEAPLNPK","psms":7,"len":18,"native_shared_substring_aa":16}]),
  entrapment_hits_are_genuinely_false=dict(
    distinct_entrapment_proteins_touched=215,
    proteins_with_3_or_more_distinct_peptides=0,
    proteins_with_exactly_one_peptide=211,
    plant_specific_abundance_markers_found=0,
    conserved_families_among_reviewed_hits={"tubulin_alpha":9,"actin":6,"histone":2},
    exact_presence_of_top_hits_in_native_fasta=0,
    verdict="entrapment labels are sound; the hits are homology-driven false matches, not real plant protein"),
  effective_sample_size_q_0_01_seed1=dict(
    entrapment_target_psms=258, distinct_entrapment_target_peptides=150,
    entrapment_decoy_psms=133, distinct_entrapment_decoy_peptides=115,
    psm_level_ratio=1.940, psm_level_binomial_z=6.32,
    peptide_level_ratio_by_seed=[1.439,1.350,1.438,1.543,1.504],
    peptide_level_binomial_z_by_seed=[3.12,2.50,3.09,3.67,3.50]),
  uncertainty=dict(
    audit_reported_mean=0.018103981884256583,
    audit_reported_seed_sd=0.000579,
    seed_sd_is_reproducibility_not_sampling_error=True,
    unweighted_mean_over_6_runs=0.01549, sd_over_6_runs=0.00582,
    t_interval_95_over_6_runs=[0.00890,0.02208],
    cluster_bootstrap_by_run_95=[0.01542,0.02070],
    cluster_bootstrap_by_peptide_95=[0.01253,0.02472],
    sensitivity_to_f={"f=1.0 (no adjustment, lower bound)":0.01320,
                      "f=0.7389 (tail decoys, as audited)":0.01787,
                      "f=0.7823 (all decoys)":0.01687,
                      "f=0.5 (design's stated balance)":0.02640}),
  classification="ESTIMATOR LIMITATION (semi-supervised TDC) + ENTRAPMENT-MEASUREMENT AMPLIFICATION; not an implementation defect in percolator-rs",
)
p=SP/"entrapment_rootcause_results.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("wrote",p)
print("\nheadline identity check (percolator-rs top5 seed1, q<0.01):")
r=qtab(arms["percolator_rs_top5_seed1"])[2]
print(f"  R={r['accepted_targets']} R_ent={r['pure_entrapment_targets']} D_ent={r['pure_entrapment_decoys']} "
      f"D_nat={r['native_decoys']} f={r['effective_entrapment_fraction']:.4f}")
print(f"  adjusted_fdp = {r['adjusted_fdp']:.5f}")
print(f"  (D_ent+D_nat)/R * R_ent/D_ent = "
      f"{(r['pure_entrapment_decoys']+r['native_decoys'])/r['accepted_targets']*r['entrapment_target_over_decoy']:.5f}")
print("\narm summary at q<0.01:")
for k,v in arms.items():
    r=qtab(v)[2]
    print(f"  {k:<42} R={r['accepted_targets']:>6} Rent/Dent={r['entrapment_target_over_decoy']:.3f} "
          f"adjFDP={r['adjusted_fdp']:.5f} z={r['binomial_z']:+.2f}")
