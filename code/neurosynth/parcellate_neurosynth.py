#!/usr/bin/env python3
"""
Run Neurosynth meta-analyses and parcellate association z-maps into Schaefer200x17 parcels.
Based on Joëlle 's https://github.com/PennLINC/tractmaps/blob/main/code/data_prep/tm_parcellate.py
and Justine Hansen's https://github.com/netneurolab/ipn-summer-school/tree/main/lectures/2021-07-02/13-15
"""

import os
import gzip
import json
import pickle
import requests
import pandas as pd
from pathlib import Path
from scipy import stats as sstats
import nibabel as nib
from nimare.extract import fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.meta.cbma.mkda import MKDAChi2
from neuromaps.datasets import fetch_atlas
from neuromaps import images
from neuromaps.parcellate import Parcellater
from nilearn.datasets import fetch_atlas_schaefer_2018



########################## 
###### Define Inputs# ####
##########################

root = Path("/Users/audluo/PennLINC/sex_diff_local") # some local path since this doesn't work on penn VPN

ns_dir = root / "input/neurosynth/raw"
ns_dir.mkdir(parents=True, exist_ok=True)
parc_dir = root / "input/neurosynth_annotations/schaefer200"
parc_dir.mkdir(parents=True, exist_ok=True)

ma_images = ["z_desc-association"]  # Neurosynth association (z) maps
email_address = "audrey.luo@pennmedicine.upenn.edu"  # optional if fetching abstracts



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
            ids = database.get_studies_by_label(f"terms_abstract_tfidf__{word}", 0.001)
            dset_sel = database.slice(ids)
            dset_unsel = database.slice(list(set(database.ids) - set(ids)))

            mkda = MKDAChi2(kernel__r=10)
            results = mkda.fit(dset_sel, dset_unsel)
            results.save_maps(path, names=ma_images)

        generated.append(path)

        print(f"[{i}/{total_terms}] {word}")

    print(f"✓ Completed meta-analyses for {total_terms} terms!")
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

    labels = [lbl.decode("utf-8") for lbl in schaefer["labels"]]

    # Create Parcellater object
    schaefer_parc = Parcellater(
        parcellation=schaefer["maps"],
        space="MNI152",
        resampling_target="parcellation"
    )

    schaefer_maps = {}
    for map_name, value in cortical_maps.items():
        print(f"Parcellating: {map_name}")

        parcellated_map = schaefer_parc.fit_transform(
            data=value["map"],
            space=value["annotation"][2]
        )

        if zscore:
            parcellated_map = sstats.zscore(parcellated_map, nan_policy="omit")
        schaefer_maps[map_name] = parcellated_map.ravel()

        data = images.load_data(value["map"])
        print(f"Original: {data.shape}, Parcellated: {parcellated_map.shape}")

    return schaefer_maps, labels





#############################################
# Fetch and Convert Neurosynth Data  
#############################################

print("Fetching Neurosynth database...")
files = fetch_neurosynth(
    data_dir=str(ns_dir),
    version="7",
    overwrite=False,
    source="abstract",
    vocab="terms",
)
neurosynth_db = files[0]

output_file = ns_dir / "neurosynth/neurosynth_dataset.pkl.gz"
output_file.parent.mkdir(parents=True, exist_ok=True)

print("Converting Neurosynth database to NiMARE Dataset...")
if not output_file.exists():
    neurosynth_dset = convert_neurosynth_to_dataset(
        coordinates_file=neurosynth_db["coordinates"],
        metadata_file=neurosynth_db["metadata"],
        annotations_files=neurosynth_db["features"],
    )
    neurosynth_dset.save(output_file)
else:
    print("Dataset already exists; skipping conversion.")

print("Loading NiMARE dataset...")
with gzip.open(output_file, "rb") as f:
    neurosynth_dset = pickle.load(f)

##############################################
# Run meta analysis
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
        "annotation": ("cognitive atlas", term, "MNI152"),
        "map": str(root_dir / term / f"z_desc-association.nii.gz"),
    }
    for term in os.listdir(root_dir)
    if os.path.isdir(root_dir / term)
}

print(f"Loaded {len(cogatl_maps)} cognitive terms.")