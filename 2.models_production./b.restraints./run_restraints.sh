#!/bin/bash
#SBATCH --job-name=restraints
#SBATCH --mail-type=ALL                                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=martina.audagnotto@astrazeneca.com   # Where to send mail
#SBATCH --time=24:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

a3m=6y75.a3m
npz=6y75.npz

~/.conda/envs/tensorflow_env/bin/python distance_dihedral_angle_restraints_prob0.05.py $a3m $npz
