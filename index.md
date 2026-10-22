---
layout: page
title: Reproducibility Guide
---

# I. Project Information

**Abstract**

Sex differences in human functional connectivity remain debated, and it is unclear whether reported differences follow a common organizing principle across the cortex. Here, we tested whether these differences align with the sensorimotor-association (S-A) axis, which spans unimodal sensorimotor to transmodal association cortices. We examined four large-scale, independent datasets of youth ages 5-23 years (total *N* = 3,355). Across datasets, males showed greater connectivity across the cortex, particularly between functional networks, whereas females showed greater connectivity within networks, most prominently within default mode subnetwork A. More broadly, the sensorimotor and association ends of the S-A axis were more segregated in females and more integrated in males. This difference extended to macroscale cortical organization: the cortical hierarchy was more expanded in females than in males along the principal gradient of functional connectivity. These replicable and generalizable results establish that sex differences in functional connectivity are organized along the cortical hierarchy.

**Team**

| **Role** | **Name** |
| --- | --- |
| **Project Lead** | Audrey C. Luo |
| **Faculty Lead** | Theodore D. Satterthwaite |
| **Analytic Replicator** | Meike Hettwer |
| **Collaborators** | Audrey Luo, Meike Hettwer, Aaron Alexander-Bloch, Joëlle Bagautdinova, Deanna Barch, Andrew Chen, Matthew Cieslak, Elizabeth Flook, Alexandre Franco, Raquel Gur, Ruben Gur, Arielle Keller, Gregory Kiar, Michael Milham, Laura Pritschet, Armin Raznahan, Taylor Salo, Giovanni Salum, Golia Shafiei, Russell Shinohara, Leah Somerville, Valerie Sydnor, Simon Vandekar, Xinyu Zhang, Sheila Shanmugan*, Theodore Satterthwaite* |

**Project Timeline**

| **Project Start Date** | October 2025 |
| --- | --- |
| **Current Project Status** | In prep |

**Code and Communication**

| **Github Repository** | [https://github.com/PennLINC/luo_sex_diff/tree/main](https://github.com/PennLINC/luo_sex_diff/tree/main) |
| --- | --- |
| **Slack Channel** | #network_replication |

**Datasets**

- RBC PNC (health exclude) as discovery dataset, NKI (replication) HCP-D (replication), and HBN (replication)

**Conference Presentations**

- Poster presented at The Society of Biological Psychiatry Annual Meeting, April 2026
- Poster presented at The Organization for Human Brain Mapping Annual Meeting, June 2026

# II. Github Repo Structure


| Folder | Description |
| --- | --- |
| `code/` | where code lives |
| `software/` | various softwares and code for installing them and setting up the project folder |
|  |  |
 

# III. Code Documentation

***Overview of Workflow***

| **Step** | **Task** | **Notes** |
| --- | --- | --- |
| 0 | [Install required software and packages](#0-install-required-software-and-packages) | |
| 1 | [Setup project directory](#1-setup-project-directory) | Make required folders and copy over sample construction files to project directory |
| 2 | [Sample construction](#2-sample-construction) | Sample same as in [Luo et al., 2024](https://www.nature.com/articles/s41467-024-47748-w) |
| 3 | [Format edge and network-pair labels](#3-format-edge-and-network-pair-labels) | |
| 4 | [Compute connectivity measures](#4-compute-connectivity-measures) | Compute edge-level, region-level (functional connectivity strength, within-network and between-network connectivity), and network-level connectivity |
| 5 | [Harmonize multi-site data using CovBat](#5-harmonize-multi-site-data-using-covbat) | Harmonize connectivity metrics computed above and connectivity matrices (for gradient analysis). HCP-D and HBN need to be harmonized. Furthermore, the pooled dataset (combined PNC, HCP-D, NKI, and HBN) also needs to be harmonized. |
| 6 | [Fit GAMs](#6-fit-gams) | |
| 7 | [Compute effect size with RESI](#7-compute-effect-size) | |
| 8 | [Parcellate Neurosynth](#8-parcellate-neurosynth) | |
| 9 | [Significance testing: spin tests](#9-significance-testing-spin-tests) | Spin testing for Neurosynth terms with WNC and BNC maps and for edge-level analysis |
| 10 | [Gradient analysis](#10-gradient-analysis) | |
| 11 | [Visualize results!](#11-visualize-results) | |

### 0. Install required software and packages

- `software/installation/` contains code for:
    - pulling the Docker container for R packages used for statistical analyses (`r-packages-for-parcc_0.0.1.sif` )
    - creating the Python environment used in this project (`sex_diff_py311`)
    - on local computer, installing packages used in Rmd files for figure-making  
    

### 1. Setup project directory
1. Make input and output directories within the replication directory
    
    ```bash
    	cd ~/code/setup
    	bash make_folders.sh
    ```
    
2. Copy over demographic files to replication directory and reformat using `prep_demographic_info.R` 
    - note that this code also sets male as the reference, so after fitting GAMs and using RESI: `positive beta/effect size = Female > Male` and 
    `negative beta/effect size = Male > Female`
    
    

### 2. Sample construction 

- The sample used in this paper is identical to what was used in [Luo et al., 2024](https://www.nature.com/articles/s41467-024-47748-w).
    - The step above copied these files into our project directory

### 3. Format edge and network-pair labels 

- `~/code/labels/make_networkpair_and_edge_labels.R` makes labels for pairs of networks and pairs of regions (i.e. edge labels)

### 4. Compute connectivity measures 

- Using connectivity matrices generated in [Luo et al., 2024](https://www.nature.com/articles/s41467-024-47748-w) (concatenated resting-state and task fMRI in PNC, HCP-D, and HBN; resting-only in NKI due to availability of data), compute edge-level, region-level (functional connectivity strength, within-network and between-network connectivity), and network-level connectivity
    
    ```bash
    cd ~/code/connectivity_measures
    bash submit_compute_connectivity_metrics.sh
    ```
    
- **Output**: subject-by-parcel  `csv's` for each region-level measure, subject-by-edge RData file, and subject-by-networkpair files
    - once jobs are finished, you can inspect the output of an example `.csv`  
    
    ```bash
    cd ~/output/PNC/
    cd FC_strength # or choose your favorite measure!
    head FC_strength_subxparcel_matrix_schaefer200x17.csv
    ```
    

### 5. Harmonize multi-site data using CovBat  
- A few things need to be harmonized: functional connectivity metrics (FC strength, within- and between- network connectivity, network-pair connectivity, and edges) and the connectivity matrices themselves (for the gradient analysis later on).
    - The datasets that need to be harmonized are HCP-D and HBN, as well as the pooled dataset (`all_datasets`) that combines PNC, NKI, HCP-D, and HBN
    
    ```bash
    cd ~/code/covbat_harmonization
    
    # submit covbat for functional connectivity metrics in 
    # HBN and HCP-D - the multi-site datasets
    bash submit_covbat_connectivity_metrics.sh  
    
    # submit covbat harmonization for functional connectivity metrics 
    # from the pooled dataset (PNC, NKI, HCPD, HBN combined)
    bash submit_covbat_all_datasets.sh  
    
    # submit covbat for connectivity matrices in HBN, HCP-D, and 
    # pooled dataset
    bash submit_covbat_fc_matrices.sh
    ```
    
- **Output:**
    - `submit_covbat_connectivity_metrics.sh` and `covbat_connectivity_metrics.R`:
        - harmonized subject-by-parcel WNC and BNC, subject-by-edge, and subject-by-networkpair for `HCP-D` and `HBN`
        
        ```bash
        ROOT="~/parcc/network_replication"
        dataset="HCPD" #for example
        metric="edge" #for example.
         
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/$metric"
        cd $OUTDIR
        ls 
        # should see subxedge_schaefer200x17_covbat.RData
        # and original file: subxedge_schaefer200x17.RData
        
        ```
        
    - `submit_covbat_all_datasets.sh` and `covbat_all_datasets.R`:
        - harmonized subject-by-parcel WNC and BNC, subject-by-edge, and subject-by-networkpair for `all_datasets`
        
        ```bash
        ROOT="~/parcc/network_replication"
        dataset="all_datasets" 
        metric="edge" #for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/$metric"
        cd $OUTDIR
        ls 
        # should see subxedge_schaefer200x17_covbat.RData
        ```
        
    - `submit_covbat_fc_matrices.sh` and `covbat_fc_matrices.R`:
        - harmonized connectivity matrices of subjects in `all_datasets` (aka all subjects across PNC, HCP-D, NKI, and HBN)
        
        ```bash
        ROOT="~/parcc/network_replication"
        dataset="all_datasets" 
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/fc_matrices/matrices_covbat"
        cd $OUTDIR
        ls 
        # should see a bunch of RData files that are the 
        # covbat harmonized FC matrices
        # of each subject
        
        # there should be 3355 RData files (for 3355 total subjects)
        find . -maxdepth 1 -type f -name '*.RData' | wc -l
        ```
        

### 6. Fit GAMs  

- Fit GAMs on PNC, NKI, HCP-D (covbat harmonized), HBN (covbat harmonized) and the pooled dataset (covbat harmonized)
- The GAM scripts do the following:
    1. Examine the main effect of sex: 
    `FC_measure ~ s(age, k=3, fx=T) + sex + head_motion`
    2. Examine the sex-by-age interaction: 
    `FC_measure ~ s(age, k=3, fx=T) + s(age, by=sex, k=3, fx=T) + sex + head_motion`
    3. Examine developmental trajectory for female and males (same equation as #2)
    
    ```bash
    cd ~/code/fit_GAMs
    
    # submit GAMs for edge-, region-, and network-level connectivity for 
    # PNC, NKI, HCP-D, and HBN
    bash submit_fit_GAMs_sex_diff.sh
    
    # submit GAMs for edge-, region-, and network-level connectivity for 
    # pooled dataset
    bash submit_fit_GAMs_all_datasets.sh
    ```
    
- **Output:**
    - GAM outputs for PNC, NKI, HCP-D, HBN, and pooled dataset (`all_datasets` )
        
        ```bash
        ROOT="~/parcc/network_replication"
        dataset="PNC" # for example. Can also check out "all_datasets"
        metric="FC_strength" # for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/GAM/$metric"
        cd $OUTDIR
        ls 
        # should see: ${dataset}_GAM_agebysex_fitted_schaefer200x17.csv  
        # ${dataset}_GAM_agebysex_interaction_schaefer200x17.csv  
        # ${dataset}_GAM_sex_maineffects_schaefer200x17.csv
        
        head ${dataset}_GAM_sex_maineffects_schaefer200x17.csv
        # should show regions + GAM output (i.e. t-values, p-values)
        ```
        

### 7. Compute effect size  

- Compute effect size using [Robust Effect Size Index (RESI)](https://statimagcoll.github.io/RESI/), a standardized measure of association that is comparable across model types and independent of sample size.
    - Positive values indicate female > male connectivity and negative values indicate male > female connectivity.
- This script fits developmental edge-, region-, and network-level GAMs on functional connectivity data, then it computes an effect size using RESI (with the `resiPEse_Th1` and `t2S` methods)
    - `t2S()` gives you the RESI point estimate from the GAM's t-statistic; used for edges since this approach is very fast.
    - `resiPEse_Th1()` estimates the RESI and also estimates its uncertainty (SE and confidence interval); used for region-level and network-level connectivity.
    - Effect size estimates from the two are highly correlated

```bash
cd ~/code/compute_effect_size

# submit script to compute effect size
bash submit_resi_sex_diff.sh
```

- **Output:**
    - RESI outputs for PNC, NKI, HCP-D, HBN, and pooled dataset (`all_datasets` )
        
        ```bash
        ROOT="~/parcc/network_replication"
        dataset="PNC" # for example
        metric="FC_strength" # for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/resi/$metric"
        cd $OUTDIR
        ls 
        # should see: 
        # ${dataset}_resiPEse_Th1_sex_maineffects_schaefer200x17.csv  
        # ${dataset}_t2S_sex_maineffects_schaefer200x17.csv
        
        head ${dataset}_resiPEse_Th1_sex_maineffects_schaefer200x17.csv
        # should show regions + RESIoutput (i.e. Estimate, RESI_Xinyu, etc.)
        ```
        

### 8. Parcellate Neurosynth 

- Run Neurosynth meta-analyses and parcellate association z-maps into Schaefer200x17 parcels.
    - Based on Joëlle Bagautdinova’s [https://github.com/PennLINC/tractmaps/blob/main/code/data_prep/tm_parcellate.py](https://github.com/PennLINC/tractmaps/blob/main/code/data_prep/tm_parcellate.py) and Justine Hansen's [https://github.com/netneurolab/ipn-summer-school/tree/main/lectures/2021-07-02/13-15](https://github.com/netneurolab/ipn-summer-school/tree/main/lectures/2021-07-02/13-15) code

```bash
cd ~/code/neurosynth
bash submit_parcellate_neurosynth.sh
```

- **Output:**
    
    ```bash
    cd ~/input/neurosynth_annotations/schaefer200
    head schaefer200x17_neurosynth_125terms.csv
    ```
    

### 9. Significance testing: spin tests 

1. **Generate spherical rotations for Schaefer200x17 for spin testing (~10 min)**

```bash
cd ~/code/spin_test
bash submit_rotate_schaefer200x17.sh
```

- **Output:**
    
    ```bash
    cd ~/software/rotate_parcellation
    ls # should see :
    # schaefer200x17.coords_sphericalrotations_N10k_seed10.rds  
    # schaefer200x17_sphericalcoords.csv
    ```
    
1. **Spin testing for Neurosynth:** correlate each Neurosynth term map with the regional sex-effect map (for within- and between-network connectivity sex effect maps) and use spin testing to test for significance for each term. Note: FDR correction happens later.  

```bash
cd ~/code/spin_test
bash submit_spin_test_neurosynth.sh
```

- **Output:**
    
    ```bash
    cd ~/output/all_datasets/spin_test/neurosynth
    ls # should see:
    # BNC_neurosynth_spin_test_results_overlap.csv
    # WNC_neurosynth_spin_test_results_overlap.csv
    ```
    
1. **Spin testing for edge-level analysis**
- Spin test to see how much of the variance in edge-wise sex effects is captured by smooth hierarchical organization -- i.e. how strongly the S-A axis explains the spatial distribution of sex differences in edge-level connectivity
- Region-level S-A labels were spatially rotated across the cortical surface and reassigned to the endpoints of each edge. A GAM models the edgewise sex effect as a smooth function of the S-A positions of its two endpoints: 
`y ~ te(p2_SA, p1_SA)` , i.e. S-A rank of parcel 2 and 1 respectively
    - Model addresses: how much of the spatial variation in edge-level sex effects is captured by where the two endpoints fall along the S-A axis?
- For each of 10,000 rotations, we refit the GAM relating edgewise sex effects to the S-A positions of the two connected parcels, generating a null distribution of explained variance. Observed R-sq was then compared with this spatially constrained null distribution.

```bash
cd ~/code/spin_test
bash submit_spin_test_SA_edges.sh
```

- **Output:**
    
    ```bash
    ROOT="~/parcc/network_replication"
    dataset="PNC" # for example
    metric="edge"  
    OUTDIR="$ROOT/covariate_analyses/sex_diff/output/$dataset/spin_test/$metric"
    cd $OUTDIR
    ls 
    # should see: 
    # SA_edges_schaefer200x17_spin_test_results.RData
    ```
    

### 10. Gradient analysis

1. **Generate group template**: build group FC template, fit PCA, sign-check PC1 against the S-A axis, and save the sign-corrected reference  
    
    ```bash
    cd ~/code/gradient_analysis
    bash submit_group_template.sh
    ```
    
    - **Output:**
        
        ```bash
        # Check group-template output
        ROOT="~/parcc/network_replication"
        dataset="PNC" # for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/${dataset}/gradient_dispersion"
      
        # check that output files were written
        cd "$OUTDIR"
        ls
        
        # Should see:
        # group_gradients_full_${dataset}.npy
        # group_pc1_${dataset}.npy
        ```
        
2. **Procrustes alignment:** align every subject's gradient to the sign-checked group template, QC the alignment  

    
    ```bash
    cd ~/code/gradient_analysis
    bash submit_procrustes_alignment.sh
    ```
    
    - **Output:**
        
        ```bash
        # check Procrustes alignment output
        ROOT="~/parcc/network_replication"
        dataset="PNC" # for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/${dataset}/gradient_dispersion"
        
        # check that output files were written
        cd "$OUTDIR"
        ls
        
        # should see additional files:
        # aligned_pc1_${dataset}.npy
        # aligned_pc1_${dataset}.csv
        # alignment_qc_${dataset}.csv
        ```
        
3. **Compute dispersion measures:** compute global, within-network, and between-network dispersion from aligned PC1 for each subject  
    
    ```bash
    cd ~/code/gradient_analysis
    bash submit_compute_dispersion.sh
    ```
    
    - **Output:**
        
        ```bash
        # Check dispersion output
        ROOT="~/parcc/network_replication"
        dataset="PNC" # for example
        OUTDIR="$ROOT/covariate_analyses/sex_diff/output/${dataset}/gradient_dispersion"
        
        # Check that output files were written
        cd "$OUTDIR"
        ls
        
        # Should see additional files:
        # dispersion_global_${dataset}.csv
        # dispersion_network_${dataset}.csv
        # pc1_loadings_long_${dataset}.csv
        # pc1_subject_summary_${dataset}.csv
        # pc1_parcel_by_sex_${dataset}.csv
        ```
        
4. **Compute RESI effect size:** compute RESI sex effect size for global, within- and between-network dispersion  
    
    ```bash
    cd ~/code/gradient_analysis
    bash submit_resi_dispersion.sh
    ```
    
- **Output:**
    
    ```bash
    # Check RESI dispersion output
    ROOT="~/parcc/network_replication"
    dataset="PNC" # for example
    OUTDIR="$ROOT/covariate_analyses/sex_diff/output/${dataset}/gradient_dispersion/resi"
    cd "$OUTDIR"
    ls
    # should see: 
    # ${dataset}_resi_dispersion_meanfc_comparison.csv
    # ${dataset}_resi_sex_dispersion_between_network.csv
    # ${dataset}_resi_sex_dispersion_global.csv
    # ${dataset}_resi_sex_dispersion_within_network.csv
    
    ```
    

### 11. Visualize results!
- Main figures: `~/code/results/main_figures.Rmd`
    - This Rmd uses functions from `~/code/results/main_figures_functions.R`
- Supp figures: `~/code/results/supp_figures.Rmd`
    - This Rmd uses functions from `~/code/results/supp_figures_functions.R`
 
