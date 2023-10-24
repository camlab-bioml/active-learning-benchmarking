#!/bin/sh
snakemake --profile slurm -j 350 --rerun-triggers mtime --rerun-incomplete --use-singularity --singularity-args "-B /home/campbell/mgeuenic/2021-whatsthatcell-analysis-Michael:/home/campbell/mgeuenic/2021-whatsthatcell-analysis-Michael"
