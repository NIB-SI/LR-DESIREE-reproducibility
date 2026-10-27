#!/usr/bin/env bash
# Group the haplotype-tagged genes into syntelogs (Synt_id + synteny category)
# with syntelogfinder: https://github.com/NIB-SI/syntelogfinder
#
# Output used downstream (polyase var_obs_file):
#   output_desiree_liftoff/genespace/De_v1.unitato_liftoff_haplotap_genespace_categories.tsv
#
# Run from a clone of the syntelogfinder repository; params files are in ./params of this folder.
# Env: nextflow, conda; MCScanX installed locally (--mcscanx_path)
set -euo pipefail

PARAMS=$(cd "$(dirname "$0")" && pwd)/params
MCSCANX=MCScanX

# Désirée
nextflow run main.nf -resume \
    -params-file "$PARAMS/desiree_liftoff.json" \
    -c conf/nextflow.config \
    -profile conda \
    --run_blast \
    --mcscanx_path "$MCSCANX"

# Atlantic (Désirée vs. Atlantic comparison)
nextflow run main.nf -resume \
    -params-file "$PARAMS/atlantic_liftoff.json" \
    -c conf/nextflow.config \
    -profile conda \
    --run_blast \
    --mcscanx_path "$MCSCANX"
