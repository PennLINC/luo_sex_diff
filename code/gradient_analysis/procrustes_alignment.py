#!/usr/bin/env python3
"""
procrustes_alignment.py -- align every subject's gradient to the
sign-checked group template, QC the alignment, check n_components
sensitivity.

Usage: python procrustes_alignment.py <DATASET>

OUTPUT:
- ALIGNED_PC1_OUT (aligned_pc1_<DATASET>.npy): Aligned PC1 values for every subject; shape is subjects*200 parcels.
- ALIGNED_PC1_OUT (aligned_pc1_<DATASET>.csv): CSV version of the same aligned subject-by-parcel PC1 matrix.
- QC_OUT (alignment_qc_<DATASET>.csv(): Subject-level QC containing correlation with the group PC1, correlation with the S-A axis, and whether the subject was flagged for low alignment.
"""

import os
import sys
import numpy as np
import pandas as pd
from brainspace.gradient import GradientMaps

DATASET = sys.argv[1]
CONFIG = dict(approach="pca", kernel="normalized_angle",
             n_components=10, random_state=0)
SPARSITY = 0.0

ROOT = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/covariate_analyses/sex_diff"
CACHE_DIR = f"{ROOT}/cache"
OUTPUT_DIR = f"{ROOT}/output/{DATASET}/gradient_dispersion"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MATS_CACHE = f"{CACHE_DIR}/all_matrices_{DATASET}.npy"
IDS_CACHE = f"{CACHE_DIR}/subject_ids_{DATASET}.csv"
GROUP_GRADIENTS_IN = f"{OUTPUT_DIR}/group_gradients_full_{DATASET}.npy"

ALIGNED_PC1_OUT = f"{OUTPUT_DIR}/aligned_pc1_{DATASET}.npy"
QC_OUT = f"{OUTPUT_DIR}/alignment_qc_{DATASET}.csv"

SA_AXIS_CSV = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/SAaxis/schaefer200x17_SAaxis.csv"
SA_PARCEL_LABELS = "/ceph/projects/sattertt/pennlinc-parcc/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv"


def load_sa_axis():
    sa = pd.read_csv(SA_AXIS_CSV, index_col=0)
    labels = pd.read_csv(SA_PARCEL_LABELS, index_col=0)
    sa["label_norm"] = sa["label"].str.replace(r"^17Networks_", "Networks_", regex=True)
    labels_norm = pd.Series(labels.index).str.replace(r"^17Networks_", "Networks_", regex=True)
    sa_reordered = sa.set_index("label_norm").loc[labels_norm].reset_index()
    return sa_reordered["SA.axis_rank"].to_numpy().astype(np.float64)


def main():
    print(f"{'=' * 60}\nDATASET: {DATASET}\n{'=' * 60}\n")

    mats = np.load(MATS_CACHE, mmap_mode="r")
    ids = pd.read_csv(IDS_CACHE)["sub"].tolist()
    reference = np.load(GROUP_GRADIENTS_IN)
    group_pc1 = reference[:, 0]
    sa_axis = load_sa_axis()
    print(f"{len(ids)} subjects, reference {reference.shape}")

    # n_iter=1 is critical: default (10) runs generalized Procrustes,
    # which replaces the reference with the running mean of aligned
    # subjects after the first pass -- drifting off the sign-checked
    # template from group_template.py.
    print(f"\nfitting + aligning {len(mats)} subjects "
          f"(n_components={CONFIG['n_components']}, n_iter=1)...")
    gm_sub = GradientMaps(alignment="procrustes", **CONFIG)
    gm_sub.fit(list(mats), sparsity=SPARSITY, reference=reference, n_iter=1)

    aligned = np.stack(gm_sub.aligned_)
    aligned_pc1 = aligned[:, :, 0]
    print(f"aligned_pc1 shape: {aligned_pc1.shape}")

    r_group = np.array([np.corrcoef(aligned_pc1[i], group_pc1)[0, 1]
                       for i in range(len(ids))])
    r_sa = np.array([np.corrcoef(aligned_pc1[i], sa_axis)[0, 1]
                     for i in range(len(ids))])

    print(f"\nQC: corr(aligned PC1, group PC1)")
    print(f"  mean {r_group.mean():.4f}  median {np.median(r_group):.4f}  "
          f"range {r_group.min():.4f}-{r_group.max():.4f}")
    print(f"QC: corr(aligned PC1, S-A axis)")
    print(f"  mean {r_sa.mean():.4f}  median {np.median(r_sa):.4f}  "
          f"range {r_sa.min():.4f}-{r_sa.max():.4f}")

    flagged = r_group < 0.3
    print(f"\nflagged (r_group < 0.3): {flagged.sum()} subjects "
          f"({flagged.mean()*100:.1f}%)")
    if flagged.sum() > 0:
        print(f"  flagged IDs (first 10): {[ids[i] for i in np.where(flagged)[0][:10]]}")

    print(f"\n{'-' * 60}\nn_components sensitivity\n{'-' * 60}")
    for k in (3, 5):
        cfg_k = {**CONFIG, "n_components": k}
        ref_k = reference[:, :k]
        gm_sub_k = GradientMaps(alignment="procrustes", **cfg_k)
        gm_sub_k.fit(list(mats), sparsity=SPARSITY, reference=ref_k, n_iter=1)
        aligned_k = np.stack(gm_sub_k.aligned_)[:, :, 0]
        r_k = np.array([np.corrcoef(aligned_pc1[i], aligned_k[i])[0, 1]
                        for i in range(len(ids))])
        print(f"  k={k:<3} vs k=10: mean |r| = {np.abs(r_k).mean():.4f}  "
              f"range {np.abs(r_k).min():.4f}-{np.abs(r_k).max():.4f}")

    np.save(ALIGNED_PC1_OUT, aligned_pc1)
    np.savetxt(ALIGNED_PC1_OUT.replace(".npy", ".csv"), aligned_pc1, delimiter=",")
    qc_df = pd.DataFrame({"sub": ids, "r_group_pc1": r_group,
                         "r_sa_axis": r_sa, "flagged": flagged})
    qc_df.to_csv(QC_OUT, index=False)
    print(f"\nsaved aligned PC1 -> {ALIGNED_PC1_OUT} (+ .csv for spin_test_dispersion.R)")
    print(f"saved QC table    -> {QC_OUT}")


if __name__ == "__main__":
    main()
