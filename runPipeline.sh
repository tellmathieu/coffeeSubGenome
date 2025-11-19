#!/bin/sh

module load conda
mamba activate braker3

snakemake --snakefile 2snake.smk -j ${SLURM_CPUS_PET_TASK} --rerun-incomplete --use-conda