#!/bin/bash
#SBATCH --job-name=peptide
#SBATCH --mail-type=ALL                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=martina.audagnotto@astrazeneca.com   # Where to send mail
#SBATCH --time=24:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
module load Rosetta

target=6y75.fasta

BuildPeptide.mpi.linuxgccrelease -in:file:fasta $target -out:file:o $target'_initial_peptide.pdb'

./run_set_phi_psi.sh $target'_initial_peptide.pdb'

