#!/bin/bash
#SBATCH --job-name=3Dpred
#SBATCH --mail-type=ALL                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=martina.audagnotto@astrazeneca.com   # Where to send mail
#SBATCH --time=24:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-10                                     #Array range
module load Rosetta 

#source bashrc with conda.sh, then conda can be used
echo ". /opt/scp/software/Anaconda2/5.2.0/etc/profile.d/conda.sh" >> ~/.bashrc
source ~/.bashrc
#make sure that conda base is activated
conda activate audagnotto_kblr332v2

#rename the files
#python rename_minfiles.py

#create the resfile
python resfile_generation_alatogly.py min_1.pdb
python resfile_generation_glytoala.py min_1.pdb

#run the procedure
./run_mutateGLYtoALA_min_relax.sh 'min_'$SLURM_ARRAY_TASK_ID'.pdb'

