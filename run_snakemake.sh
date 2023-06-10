#!/bin/sh
source activate types-pipeline
snakemake  --resources mem_mb=100000 --until all --cores 30
