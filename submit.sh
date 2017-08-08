#!/bin/bash
#SBATCH -J redoanalysis
#SBATCH -o /n/regal/debivort_lab/claire/bees/analysis.out
#SBATCH -e /n/regal/debivort_lab/claire/bees/analysis.err
#SBATCH -N 1 		
#SBATCH -c 32 		
#SBATCH -t 4-00:00 	
#SBATCH -p general 	
#SBATCH --mem=32000	
#SBATCH --mail-type=END
#SBATCH --mail-user=guerin.claire01@gmail.com

mkdir -p /scratch/$USER/$SLURM_JOB_ID

srun -c 32 matlab-default -nosplash -nodesktop -r "trackchicks;exit"

rm -rf /scratch/$USER/$SLURM_JOB_ID
