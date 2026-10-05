#!/usr/bin/env bash

# Submit edgewise S-A spin tests
ROOT="/ceph/projects/sattertt/pennlinc-parcc/network_replication"
PROJDIR="$ROOT/covariate_analyses/sex_diff"

APPTAINER="/vast/parcc/spack/sw/apps/linux-sapphirerapids/apptainer-1.4.1-qpr4lterya7ontg7pkqrx7b3jkab3lcw/bin/apptainer"
SIF="$PROJDIR/software/r_packages/r-packages-for-parcc_0.0.1.sif"
r_script="$PROJDIR/code/spin_test/spin_test_SA_edges.R"

datasets=("PNC" "HCPD" "NKI" "HBN" "all_datasets")
atlases=("schaefer200x17")

cpus=1
time_limit="12:00:00"

for dataset in "${datasets[@]}"; do
    for atlas in "${atlases[@]}"; do

        # Where to save output and error logs
        logs_dir="$PROJDIR/code/logs/spin_test/${dataset}"
        mkdir -p "$logs_dir"

        # Make spin-test output directory
        outputs_root="$PROJDIR/output/${dataset}/spin_test"
        mkdir -p "$outputs_root"

        job_name="${dataset}_${atlas}_SA_edge_spin_test"

        if sbatch \
            --job-name="$job_name" \
            --partition=genoa-std-mem \
            --nodes=1 \
            --ntasks=1 \
            --cpus-per-task="$cpus" \
            --time="$time_limit" \
            --output="${logs_dir}/${job_name}_%j.out" \
            --error="${logs_dir}/${job_name}_%j.err" \
            --wrap="$APPTAINER exec --cleanenv --bind \"$ROOT:$ROOT\" \"$SIF\" Rscript \"$r_script\" \"$dataset\" \"$atlas\""
        then
            echo "Submitted ${dataset} ${atlas} (${cpus} CPU(s), ${time_limit})"
        else
            echo "FAILED to submit ${dataset} ${atlas}" >&2
        fi

    done
done
 