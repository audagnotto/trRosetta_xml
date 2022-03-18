#!/bin/bash
#SBATCH --job-name=MIN
#SBATCH --mail-type=ALL                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=martina.audagnotto@astrazeneca.com   # Where to send mail
#SBATCH --time=24:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
module load Rosetta

ls ../../a.initial_peptide/*_initial_peptide_phipsi_*.pdb > list_min
minimize.mpi.linuxgccrelease @minimizer_flags
