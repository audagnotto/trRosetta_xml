#!/bin/bash

#source bashrc with conda.sh, then conda can be used
echo ". /opt/scp/software/Anaconda2/5.2.0/etc/profile.d/conda.sh" >> ~/.bashrc
source ~/.bashrc
#make sure that conda base is activated
source activate tensorflow_env

cp ../0.seq/6y75.hhba3m 6y75.a3m

network=/projects/cc/audagnotto_kblr332/sync/trRosetta/network
model=/projects/cc/audagnotto_kblr332/sync/trRosetta/model2019_07

a3m=6y75.a3m
npz=6y75.npz

~/.conda/envs/tensorflow_env/bin/python $network/predict.py -m $model $a3m $npz
