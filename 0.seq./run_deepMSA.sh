#!/bin/bash

hhsuit2=/projects/cc/audagnotto_kblr332/sync/MSA/deepMSA/hhsuite2/scripts
db=/projects/cc/audagnotto_kblr332/sync/MSA/deepMSA/databases/uniclust30_2017_10_hhsuite/uniclust30_2017_10/uniclust30_2017_10

python $hhsuit2/build_MSA.py 6y75.fasta -hhblitsdb=$db
