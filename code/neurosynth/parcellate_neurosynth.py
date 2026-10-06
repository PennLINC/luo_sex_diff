#!/usr/bin/env python3

"""
Run Neurosynth meta-analyses and parcellate association z-maps into
Schaefer200x17 parcels.

Based on Joëlle Bagautdinova's code:
https://github.com/PennLINC/tractmaps/blob/main/code/data_prep/tm_parcellate.py

and Justine Hansen's code:
https://github.com/netneurolab/ipn-summer-school/tree/main/lectures/2021-07-02/13-15
"""

import os
import gc
import json
import requests
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats as sstats
from nimare.extract import fetch_neurosynth
from nimare.meta.cbma.mkda import MKDAChi2
from neuromaps import images
from neuromaps.parcellate import Parcellater
from nilearn.datasets import fetch_atlas_schaefer_2018


##########################
###### Define Inputs #####
##########################

root = Path(
    "/ceph/projects/sattertt/pennlinc-parcc/network_replication/"
    "covariate_analyses/sex_diff/input"
)

ns_dir = root / "neurosynth/raw"
ns_dir.mkdir(parents=True, exist_ok=True)

parc_dir = root / "neurosynth_annotations/schaefer200"
parc_dir.mkdir(parents=True, exist_ok=True)

# Neurosynth association z-maps
ma_images = ["z_desc-association"]


###############
# Functions
###############

def get_cogatlas_concepts(url=None):
    """Fetch list of concepts from the Cognitive Atlas API."""
    if url is None:
        url = "https://cognitiveatlas.org/api/v-alpha/concept"

    req = requests.get(url)
    req.raise_for_status()

    concepts = set([f.get("name") for f in json.loads(req.content)])
    return concepts


def run_meta_analyses(database, ma_images, use_features=None, outdir=None):
    """Run MKDA Chi-square meta-analyses for each term."""
    if outdir is None:
        outdir = ns_dir / "terms"

    outdir.mkdir(parents=True, exist_ok=True)

    all_features = set(database.get_labels())
    features = [f.replace("terms_abstract_tfidf__", "") for f in all_features]

    if use_features is not None:
        features = set(features) & set(use_features)

    generated = []
    total_terms = len(sorted(features))

    print(f"Running meta-analyses for {total_terms} terms...")

    for i, word in enumerate(sorted(features), 1):
        path = outdir / word.replace(" ", "_")
        path.mkdir(parents=True, exist_ok=True)

        if not all(os.path.exists(path / f"{f}.nii.gz") for f in ma_images):
            ids = database.get_studies_by_label(
                f"terms_abstract_tfidf__{word}",
                0.001
            )

            dset_sel = database.slice(ids)
            dset_unsel = database.slice(
                list(set(database.ids) - set(ids))
            )

            mkda = MKDAChi2(kernel__r=10)
            results = mkda.fit(dset_sel, dset_unsel)
            results.save_maps(path, names=ma_images)

            # MKDA objects can be large. Explicitly release term-specific
            # objects before starting the next meta-analysis to limit peak
            # memory use across the full set of terms.
            del results, mkda, dset_sel, dset_unsel
            gc.collect()

        generated.append(path)

        print(f"[{i}/{total_terms}] {word}")

    print(f"Completed meta-analyses for {total_terms} terms!")

    return generated


def schaeferize(cortical_maps, zscore=True, n_parcels=200, n_networks=17):
    """
    Parcellate cortical maps using the Schaefer 2018 atlas (from Nilearn).
    """

    # Fetch Schaefer atlas
    schaefer = fetch_atlas_schaefer_2018(
        n_rois=n_parcels,
        yeo_networks=n_networks,
        resolution_mm=2
    )

    # Get Schaefer parcel labels
    labels = [
        lbl.decode("utf-8") if isinstance(lbl, bytes) else str(lbl)
        for lbl in schaefer["labels"]
    ]

    # Nilearn may include a background label; remove it
    if len(labels) == n_parcels + 1:
        labels = labels[1:]

    assert len(labels) == n_parcels, (
        f"Expected {n_parcels} parcel labels, got {len(labels)}"
    )

    # Create Parcellater object
    schaefer_parc = Parcellater(
        parcellation=schaefer["maps"],
        space="MNI152",
        resampling_target="parcellation"
    )

    schaefer_maps = {}

    for map_name, value in cortical_maps.items():
        print(f"Parcellating: {map_name}")

        # Current neuromaps versions may return shape (1, n_parcels).
        # Flatten first so z-scoring is explicitly performed across parcels.
        parcellated_map = schaefer_parc.fit_transform(
            data=value["map"],
            space=value["annotation"][2]
        ).ravel()

        if parcellated_map.size != n_parcels:
            raise ValueError(
                f"{map_name}: expected {n_parcels} parcels, "
                f"got {parcellated_map.size}"
            )

        if not np.isfinite(parcellated_map).any():
            raise ValueError(
                f"{map_name}: all parcellated values are non-finite "
                "before z-scoring"
            )

        if zscore:
            parcellated_map = sstats.zscore(
                parcellated_map,
                nan_policy="omit"
            )

        if not np.isfinite(parcellated_map).any():
            raise ValueError(
                f"{map_name}: all parcellated values are non-finite "
                "after z-scoring"
            )

        schaefer_maps[map_name] = parcellated_map

        data = images.load_data(value["map"])
        print(
            f"Original: {data.shape}, "
            f"Parcellated: {parcellated_map.shape}, "
            f"finite: {np.isfinite(parcellated_map).sum()}/{n_parcels}"
        )

    return schaefer_maps, labels


#############################################
# Fetch Neurosynth Data
#############################################

print("Fetching Neurosynth database...")

files = fetch_neurosynth(
    data_dir=str(ns_dir),
    version="7",
    overwrite=False,
    source="abstract",
    vocab="terms",
    return_type="dataset",
)

neurosynth_dset = files[0]

print("Neurosynth Dataset loaded!")


##############################################
# Run meta-analysis
##############################################

generated = run_meta_analyses(
    database=neurosynth_dset,
    ma_images=ma_images,
    use_features=get_cogatlas_concepts(),
)

print("Done with meta-analyses!")


##############################################
# Parcellate Neurosynth (Schaefer 200x17)
##############################################

print("Building term map dictionary...")

root_dir = ns_dir / "terms"

cogatl_maps = {
    term: {
        "annotation": (
            "cognitive atlas",
            term.replace("_", " "),
            "MNI152"
        ),
        "map": str(root_dir / term / "z_desc-association.nii.gz"),
    }
    for term in os.listdir(root_dir)
    if os.path.isdir(root_dir / term)
    and (root_dir / term / "z_desc-association.nii.gz").exists()
}

print(f"Loaded {len(cogatl_maps)} cognitive terms.")

if len(cogatl_maps) == 0:
    raise RuntimeError("No Neurosynth term maps were found.")

print("Parcellating Neurosynth maps to Schaefer200x17...")

schaefer_maps, labels = schaeferize(
    cortical_maps=cogatl_maps,
    zscore=True,
    n_parcels=200,
    n_networks=17,
)


##############################################
# Save parcellated annotations
##############################################

# Save one row per Schaefer parcel and one column per Neurosynth term.
# Keeping the atlas labels in the first column makes the output easy to
# merge with parcel-level analyses in R.
schaefer_df = pd.DataFrame(schaefer_maps)

# Fail rather than silently writing an unusable all-NA term map.
all_missing = schaefer_df.columns[schaefer_df.isna().all()].tolist()

if all_missing:
    raise ValueError(
        f"{len(all_missing)} term(s) contain only NaN values: "
        f"{all_missing[:10]}"
    )

schaefer_df.insert(0, "label", labels)
schaefer_df.insert(0, "regionID", range(1, len(labels) + 1))

outfile = parc_dir / "schaefer200x17_neurosynth_125terms.csv"
schaefer_df.to_csv(outfile, index=False)

print(f"Saved parcellated Neurosynth annotations -> {outfile}")
print(
    f"Output shape: {schaefer_df.shape[0]} parcels x "
    f"{len(schaefer_maps)} terms"
)
print(
    f"Missing annotation values: "
    f"{schaefer_df.iloc[:, 2:].isna().sum().sum()}"
)
print("Done!")
