#!/bin/bash
#SBATCH --partition=BIG
#SBATCH --output=/SCRATCH/hosseini/INVERSION/slurm_output/inversion_%j.out
#SBATCH --workdir=/SCRATCH/hosseini/INVERSION/?????
#SBATCH --job-name="INV-?????"
#SBATCH --ntasks=264
#SBATCH --nodes=11
#SBATCH --mail-type=all
#SBATCH --time=719:59:59
# run compute job

source /home/hosseini/.bashrc
mpirun.openmpi -np 264 mpisolvetomo < in.?????
