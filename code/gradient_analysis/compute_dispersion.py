#!/usr/bin/env python3
"""
compute_dispersion.py -- global, within-network, and between-network
dispersion from aligned PC1. Saves CSVs for resi_dispersion.R.

Usage: python compute_dispersion.py <DATASET>

UPDATE: primary outcome is raw dispersion 
"""

import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

DATASET = sys.argv[1]
N_PARCELS = 200

ROOT = "/cbica/projects/network_replication/covariate_analyses/sex_diff"
CACHE_DIR = f"{ROOT}/cache"
OUTPUT_DIR = f"{ROOT}/output/{DATASET}/gradient_dispersion"

ALIGNED_PC1_IN = f"{OUTPUT_DIR}/aligned_pc1_{DATASET}.npy"
IDS_CACHE = f"{CACHE_DIR}/subject_ids_{DATASET}.csv"
MATS_CACHE = f"{CACHE_DIR}/all_matrices_{DATASET}.npy"

DEMO_PATH = {
    "PNC":  f"{ROOT}/input/PNC/sample_info/PNC_demographics_finalsample.csv",
    "NKI":  f"{ROOT}/input/NKI/sample_info/NKI_demographics_finalsample.csv",
    "HCPD": f"{ROOT}/input/HCPD/sample_info/HCPD_demographics_finalsample.csv",
    "HBN":  f"{ROOT}/input/HBN/sample_info/HBN_demographics_finalsample.csv",
}

DISP_OUT = f"{OUTPUT_DIR}/dispersion_global_{DATASET}.csv"
NET_DISP_OUT = f"{OUTPUT_DIR}/dispersion_network_{DATASET}.csv"
PC1_LONG_OUT = f"{OUTPUT_DIR}/pc1_loadings_long_{DATASET}.csv"
PC1_SUBJECT_SUMMARY_OUT = f"{OUTPUT_DIR}/pc1_subject_summary_{DATASET}.csv"
PC1_SEX_BRAIN_OUT = f"{OUTPUT_DIR}/pc1_parcel_by_sex_{DATASET}.csv"

SA_PARCEL_LABELS = "/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv"

NET_MAP_17 = [
    ("VisCent", "visualCentral"), ("VisPeri", "visualPeripheral"),
    ("SomMotA", "somatomotorA"), ("SomMotB", "somatomotorBAuditory"),
    ("DorsAttnA", "dorsalAttentionA"), ("DorsAttnB", "dorsalAttentionB"),
    ("SalVentAttnA", "salienceVentralAttentionA"),
    ("SalVentAttnB", "salienceVentralAttentionB"),
    ("LimbicA_TempPole", "limbicTemporopolar"),
    ("LimbicB_OFC", "limbicOrbitofrontal"),
    ("ContA", "frontoparietalControlA"), ("ContB", "frontoparietalControlB"),
    ("ContC", "frontoparietalControlC"),
    ("DefaultA", "defaultA"), ("DefaultB", "defaultB"), ("DefaultC", "defaultC"),
    ("TempPar", "temporoparietal"),
]


def assign_network(label):
    for pattern, name in NET_MAP_17:
        if pattern in label:
            return name
    raise ValueError(f"unmatched parcel label: {label}")


def load_demographics(dataset):
    if dataset == "all_datasets":
        frames = []
        for ds in ("PNC", "NKI", "HCPD", "HBN"):
            d = pd.read_csv(DEMO_PATH[ds])[["sub", "sex", "age", "meanFD_avgSes"]]
            d["dataset"] = ds
            frames.append(d)
        return pd.concat(frames, ignore_index=True)
    d = pd.read_csv(DEMO_PATH[dataset])[["sub", "sex", "age", "meanFD_avgSes"]]
    d["dataset"] = dataset
    return d


def infer_hemisphere(label):
    """Return LH/RH when hemisphere is encoded in the Schaefer parcel label."""
    if "_LH_" in label or label.startswith("LH_"):
        return "LH"
    if "_RH_" in label or label.startswith("RH_"):
        return "RH"
    return np.nan


def main():
    print(f"{'=' * 60}\nDATASET: {DATASET}\n{'=' * 60}\n")

    aligned_pc1 = np.load(ALIGNED_PC1_IN)
    ids = pd.read_csv(IDS_CACHE)["sub"].tolist()
    assert aligned_pc1.shape[0] == len(ids)

    # ---- global dispersion ----
    centroid = np.median(aligned_pc1, axis=1, keepdims=True)
    dispersion = np.sum((aligned_pc1 - centroid) ** 2, axis=1)

    df = pd.DataFrame({"sub": ids, "dispersion": dispersion})

    demo = load_demographics(DATASET)
    n_before = len(df)
    df = df.merge(demo, on="sub", how="inner")
    assert len(df) == n_before, "ID mismatch merging demographics -- check prefix format"

    mats = np.load(MATS_CACHE, mmap_mode="r")
    iu = np.triu_indices(N_PARCELS, k=1)
    df["mean_fc"] = np.array([m[iu].mean() for m in mats])

    # Global dispersion, with and without mean_fc.
    FORM = "{y} ~ bs(age, df=3) + C(sex, Treatment('Female')) + meanFD_avgSes"
    m0 = smf.ols(FORM.format(y="dispersion"), data=df).fit()
    m1 = smf.ols(FORM.format(y="dispersion") + " + mean_fc", data=df).fit()
    print("Global dispersion preview (see resi_dispersion.R for the real result):")
    print(f"  residual skew  {pd.Series(m0.resid).skew():.2f}")
    t0 = [t for t in m0.params.index if "sex" in t][0]
    t1 = [t for t in m1.params.index if "sex" in t][0]
    print(f"  sex beta (no mean_fc) {m0.params[t0]:.4f}  p {m0.pvalues[t0]:.4g}")
    print(f"  sex beta (+ mean_fc)  {m1.params[t1]:.4f}  p {m1.pvalues[t1]:.4g}")

    df.to_csv(DISP_OUT, index=False)
    print(f"saved -> {DISP_OUT}")

    # ---- PC1 subject-level summaries ----
    # Each subject's PC1 is a 200-parcel vector, so retain several descriptive
    # summaries rather than treating PC1 as if it were a single subject score.
    pc1_subject = pd.DataFrame({
        "sub": ids,
        "pc1_mean": np.mean(aligned_pc1, axis=1),
        "pc1_median": np.median(aligned_pc1, axis=1),
        "pc1_sd": np.std(aligned_pc1, axis=1, ddof=1),
        "pc1_q25": np.quantile(aligned_pc1, 0.25, axis=1),
        "pc1_q75": np.quantile(aligned_pc1, 0.75, axis=1),
        "pc1_min": np.min(aligned_pc1, axis=1),
        "pc1_max": np.max(aligned_pc1, axis=1),
    })
    pc1_subject = pc1_subject.merge(
        df[["sub", "dataset", "sex", "age", "meanFD_avgSes", "mean_fc"]],
        on="sub", how="inner")
    assert len(pc1_subject) == len(ids), "merge lost subjects in PC1 summary export"
    pc1_subject.to_csv(PC1_SUBJECT_SUMMARY_OUT, index=False)
    print(f"saved PC1 subject summaries -> {PC1_SUBJECT_SUMMARY_OUT}")

    # ---- within/between network dispersion ----
    labels = pd.read_csv(SA_PARCEL_LABELS, index_col=0).index.to_series().reset_index(drop=True)
    assert len(labels) == N_PARCELS
    net_labels = labels.apply(assign_network).to_numpy()
    net_names = sorted(set(net_labels))
    net_idx = {n: (net_labels == n) for n in net_names}

    # ---- PC1 long-format export for density/distribution plots ----
    # One row = one subject x one parcel. This preserves the full aligned PC1
    # information and can be used for pooled, sex-stratified, dataset-stratified,
    # or parcel-specific distributions in R.
    pc1_wide = pd.DataFrame(aligned_pc1, columns=np.arange(1, N_PARCELS + 1))
    pc1_wide.insert(0, "sub", ids)
    pc1_long = pc1_wide.melt(id_vars="sub", var_name="parcel", value_name="pc1")
    pc1_long["parcel"] = pc1_long["parcel"].astype(int)

    parcel_info = pd.DataFrame({
        "parcel": np.arange(1, N_PARCELS + 1),
        "parcel_label": labels.to_numpy(),
        "hemisphere": labels.apply(infer_hemisphere).to_numpy(),
        "network": net_labels,
    })
    pc1_long = pc1_long.merge(parcel_info, on="parcel", how="left")
    pc1_long = pc1_long.merge(
        df[["sub", "dataset", "sex", "age", "meanFD_avgSes", "mean_fc"]],
        on="sub", how="inner")
    assert len(pc1_long) == len(ids) * N_PARCELS, "unexpected row loss in PC1 long export"
    pc1_long.to_csv(PC1_LONG_OUT, index=False)
    print(f"saved PC1 long table -> {PC1_LONG_OUT}")

    # ---- sex-specific parcelwise PC1 for cortical brain maps ----
    # Use the mean aligned PC1 loading as the main map value. Median, SD, SE,
    # and n are included so the plotting/statistical choice can be changed later.
    pc1_sex_brain = (
        pc1_long.dropna(subset=["sex"])
        .groupby(["dataset", "sex", "parcel"], as_index=False, sort=False)["pc1"]
        .agg(pc1_mean="mean", pc1_median="median", pc1_sd="std", n="count")
    )
    pc1_sex_brain = pc1_sex_brain.merge(parcel_info, on="parcel", how="left")
    pc1_sex_brain["pc1_se"] = pc1_sex_brain["pc1_sd"] / np.sqrt(pc1_sex_brain["n"])
    pc1_sex_brain.to_csv(PC1_SEX_BRAIN_OUT, index=False)
    print(f"saved sex-specific parcel PC1 maps -> {PC1_SEX_BRAIN_OUT}")

    wn = {}
    centroids = []
    for n in net_names:
        sub_pc1 = aligned_pc1[:, net_idx[n]]
        cent = np.median(sub_pc1, axis=1)
        wn[f"wn_{n}"] = np.sum((sub_pc1 - cent[:, None]) ** 2, axis=1)
        centroids.append(cent)
    wn_df = pd.DataFrame(wn)
    wn_df.insert(0, "sub", ids)

    centroids = np.column_stack(centroids)
    bn = {}
    for i in range(len(net_names)):
        for j in range(i + 1, len(net_names)):
            bn[f"bn_{net_names[i]}__{net_names[j]}"] = np.abs(centroids[:, i] - centroids[:, j])
    bn_df = pd.DataFrame(bn)
    bn_df.insert(0, "sub", ids)

    net_all = wn_df.merge(bn_df, on="sub").merge(
        df[["sub", "dataset", "sex", "age", "meanFD_avgSes", "mean_fc"]], on="sub", how="inner")
    assert len(net_all) == len(ids), "merge lost subjects"

    outcome_cols = [c for c in net_all.columns if c.startswith(("wn_", "bn_"))]

    net_out = net_all[["sub", "dataset", "sex", "age", "meanFD_avgSes", "mean_fc"]].copy()
    for c in outcome_cols:
        net_out[c] = net_all[c]

    net_out.to_csv(NET_DISP_OUT, index=False)
    print(f"\nsaved {net_out.shape[0]} rows x {net_out.shape[1]} cols -> {NET_DISP_OUT}")
    print(f"  {sum(c.startswith('wn_') for c in outcome_cols)} within-network")
    print(f"  {sum(c.startswith('bn_') for c in outcome_cols)} between-network")


if __name__ == "__main__":
    main()
