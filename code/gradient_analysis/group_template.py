#!/usr/bin/env python3
"""
group_template.py -- build group FC template, fit PCA, sign-check PC1
against the S-A axis, save the sign-corrected reference for
procrustes_alignment.py.

Usage: python group_template.py <DATASET>
       DATASET in {PNC, NKI, HCPD, HBN, all_datasets}

OUTPUT:
- TEMPLATE_CACHE (group_template_<DATASET>.npy): Cached group-average FC matrix used to fit the group gradient.
- MATS_CACHE (all_matrices_<DATASET>.npy): Cached subject FC matrices used throughout the gradient pipeline.
- IDS_CACHE (subject_ids_<DATASET>.csv): Subject IDs corresponding, in order, to the cached FC matrices.
- GROUP_PC1_OUT (group_pc1_<DATASET>.npy): Sign-checked group PC1 vector across the 200 parcels.
- GROUP_GRADIENTS_OUT (group_gradients_full_<DATASET>.npy): Full set of group gradients/components, with PC1 sign-corrected; used as the reference for Procrustes alignment.
"""

import json
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import rdata
from brainspace.gradient import GradientMaps

warnings.filterwarnings("ignore", category=UserWarning, module="rdata")

DATASET = sys.argv[1]
N_PARCELS = 200
KEY = "schaefer217_conn"

CONFIG = dict(approach="pca", kernel="normalized_angle",
             n_components=10, random_state=0)
SPARSITY = 0.0

SA_AXIS_CSV = "/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv"
SA_PARCEL_LABELS = "/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv"

ROOT = "/cbica/projects/network_replication/covariate_analyses/sex_diff"
CACHE_DIR = f"{ROOT}/cache"
OUTPUT_DIR = f"{ROOT}/output/{DATASET}/gradient_dispersion"
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

TEMPLATE_CACHE = f"{CACHE_DIR}/group_template_{DATASET}.npy"
MATS_CACHE = f"{CACHE_DIR}/all_matrices_{DATASET}.npy"
IDS_CACHE = f"{CACHE_DIR}/subject_ids_{DATASET}.csv"
GROUP_PC1_OUT = f"{OUTPUT_DIR}/group_pc1_{DATASET}.npy"
GROUP_GRADIENTS_OUT = f"{OUTPUT_DIR}/group_gradients_full_{DATASET}.npy"

# ------------------------------------------------------------------
# Dataset -> FC matrix source.
# HCPD/HBN use covbat-harmonized matrices (multi-site). PNC/NKI are
# single-site, use raw matrices. all_datasets pools all four from the
# single combined harmonized folder (NKI IDs there lack "sub-" prefix).
# ------------------------------------------------------------------
RAW_ROOT = "/cbica/projects/network_replication/manuscript/input"
HARMONIZED_ROOT = f"{ROOT}/output"

DATASET_SOURCE = {
    "PNC":  dict(kind="raw", conn_dir=f"{RAW_ROOT}/PNC/connMatricesData/connectivity_matrices"),
    "NKI":  dict(kind="raw", conn_dir=f"{RAW_ROOT}/NKI/connMatricesData/connectivity_matrices"),
    "HCPD": dict(kind="harmonized", mats_dir=f"{HARMONIZED_ROOT}/HCPD/fc_matrices/matrices_covbat"),
    "HBN":  dict(kind="harmonized", mats_dir=f"{HARMONIZED_ROOT}/HBN/fc_matrices/matrices_covbat"),
    "all_datasets": dict(kind="harmonized", mats_dir=f"{HARMONIZED_ROOT}/all_datasets/fc_matrices/matrices_covbat"),
}

SAMPLE_FILE = {
    "PNC":  f"{ROOT}/input/PNC/sample_info/final_sex_diff_sample.txt",
    "NKI":  f"{ROOT}/input/NKI/sample_info/final_sex_diff_sample.txt",
    "HCPD": f"{ROOT}/input/HCPD/sample_info/final_sex_diff_sample.txt",
    "HBN":  f"{ROOT}/input/HBN/sample_info/final_sex_diff_sample.txt",
    "all_datasets": None,   # built by concatenating the four below
}


def load_sample_ids(dataset):
    if dataset == "all_datasets":
        ids = []
        for ds in ("PNC", "NKI", "HCPD", "HBN"):
            with open(SAMPLE_FILE[ds]) as f:
                ds_ids = [ln.strip().strip('"') for ln in f if ln.strip()]
            ds_ids = [i for i in ds_ids if re.search(r"\d", i)]
            ids.extend(ds_ids)
        return ids
    with open(SAMPLE_FILE[dataset]) as f:
        ids = [ln.strip().strip('"') for ln in f if ln.strip()]
    return [i for i in ids if re.search(r"\d", i)]


def matrix_path(dataset, sid):
    src = DATASET_SOURCE[dataset]
    if src["kind"] == "raw":
        # NKI filenames don't carry the "sub-" prefix that the sample list uses
        file_sid = sid[4:] if (dataset == "NKI" and sid.startswith("sub-")) else sid
        return os.path.join(src["conn_dir"], f"{file_sid}_ConnMatrices.RData")
    # harmonized: covbat_fc_matrices.R saves these as "{sid}_fc_covbat.RData"
    return os.path.join(src["mats_dir"], f"{sid}_fc_covbat.RData")


def load_matrix(dataset, sid):
    src = DATASET_SOURCE[dataset]
    if src["kind"] == "raw":
        raw = rdata.read_rda(matrix_path(dataset, sid))
        if dataset == "NKI":
            # NKI matrices are keyed by session, e.g. "BAS1_schaefer217_conn"
            # -- grep for whichever session key is present for this subject
            candidates = [k for k in raw.keys() if KEY in k]
            if len(candidates) != 1:
                raise ValueError(f"expected exactly 1 matching key for {sid} "
                                 f"in {list(raw.keys())}, got {candidates}")
            return np.asarray(raw[candidates[0]], dtype=np.float64)
        return np.asarray(raw[KEY], dtype=np.float64)
    # harmonized files are saved as a bare matrix (edge_vector_to_matrix
    # output), not a named list -- but fall back gracefully if that changes
    m = rdata.read_rda(matrix_path(dataset, sid))
    if isinstance(m, dict):
        m = m[list(m.keys())[0]]
    return np.asarray(m, dtype=np.float64)


def load_all_matrices(dataset):
    if os.path.exists(MATS_CACHE) and os.path.exists(IDS_CACHE):
        mats = np.load(MATS_CACHE)
        ids = pd.read_csv(IDS_CACHE)["sub"].tolist()
        print(f"loaded cached matrices: {mats.shape}")
        return mats, ids

    ids_in = load_sample_ids(dataset)
    print(f"{dataset}: {len(ids_in)} subjects in sample list(s)")
    print(f"  first ID: {ids_in[0]}")
    print(f"  first path checked: {matrix_path(dataset, ids_in[0])}")
    print(f"  exists: {os.path.exists(matrix_path(dataset, ids_in[0]))}")

    mats, used, missing = [], [], []
    for n, sid in enumerate(ids_in):
        p = matrix_path(dataset, sid)
        if not os.path.exists(p):
            missing.append(sid)
            continue
        m = load_matrix(dataset, sid)
        if m.shape != (N_PARCELS, N_PARCELS):
            print(f"  !! {sid}: unexpected shape {m.shape}, skipping")
            missing.append(sid)
            continue
        np.fill_diagonal(m, 0.0)
        mats.append(m)
        used.append(sid)
        if (n + 1) % 200 == 0:
            print(f"  {n+1}/{len(ids_in)}")

    print(f"loaded {len(used)}, missing {len(missing)}")
    if len(mats) == 0:
        raise RuntimeError(
            f"No matrices loaded for {dataset} -- all {len(ids_in)} subjects "
            f"were missing. Check matrix_path() output above against the "
            f"actual directory contents (ls the conn_dir/mats_dir)."
        )
    mats = np.stack(mats)
    np.save(MATS_CACHE, mats)
    pd.DataFrame({"sub": used}).to_csv(IDS_CACHE, index=False)
    return mats, used


def load_sa_axis():
    sa = pd.read_csv(SA_AXIS_CSV, index_col=0)
    labels = pd.read_csv(SA_PARCEL_LABELS, index_col=0)
    sa["label_norm"] = sa["label"].str.replace(r"^17Networks_", "Networks_", regex=True)
    labels_norm = pd.Series(labels.index).str.replace(r"^17Networks_", "Networks_", regex=True)
    missing = set(labels_norm) - set(sa["label_norm"])
    if missing:
        raise ValueError(f"{len(missing)} labels have no match, e.g. {list(missing)[:5]}")
    sa_reordered = sa.set_index("label_norm").loc[labels_norm].reset_index()
    assert (sa_reordered["label_norm"].to_numpy() == labels_norm.to_numpy()).all()
    return sa_reordered["SA.axis_rank"].to_numpy().astype(np.float64)


def main():
    print(f"{'=' * 60}\nDATASET: {DATASET}\n{'=' * 60}\n")

    mats, ids = load_all_matrices(DATASET)

    if os.path.exists(TEMPLATE_CACHE):
        group_fc = np.load(TEMPLATE_CACHE)
        print("loaded cached group template")
    else:
        group_fc = mats.mean(axis=0)
        np.save(TEMPLATE_CACHE, group_fc)
        print("built and cached group template")

    print(f"group_fc shape {group_fc.shape}, "
          f"symmetric {np.allclose(group_fc, group_fc.T, atol=1e-8)}")

    gm_group = GradientMaps(**CONFIG).fit(group_fc, sparsity=SPARSITY)
    group_pc1 = gm_group.gradients_[:, 0]
    print(f"group PC1 share of top-{CONFIG['n_components']} lambdas: "
          f"{gm_group.lambdas_[0] / gm_group.lambdas_.sum():.3f}")

    sa_axis = load_sa_axis()
    r = np.corrcoef(group_pc1, sa_axis)[0, 1]
    print(f"\ncorr(group PC1, S-A axis) = {r:.4f}")

    if r < 0:
        print("negative -- flipping sign")
        group_pc1 = -group_pc1
        gm_group.gradients_[:, 0] = group_pc1
        r = -r
    else:
        print("positive -- no flip needed")

    print(f"final corr(group PC1, S-A axis) = {r:.4f}")
    assert r > 0.3, "correlation too weak even after sign flip -- check parcel ordering"

    gm_group90 = GradientMaps(**CONFIG).fit(group_fc, sparsity=0.9)
    r90 = abs(np.corrcoef(group_pc1, gm_group90.gradients_[:, 0])[0, 1])
    print(f"\n(robustness) |r| group PC1 at sparsity 0 vs 0.9: {r90:.4f}")

    np.save(GROUP_PC1_OUT, group_pc1)
    np.save(GROUP_GRADIENTS_OUT, gm_group.gradients_)
    print(f"\nsaved sign-checked group PC1        -> {GROUP_PC1_OUT}")
    print(f"saved full sign-corrected reference -> {GROUP_GRADIENTS_OUT}")


if __name__ == "__main__":
    main()