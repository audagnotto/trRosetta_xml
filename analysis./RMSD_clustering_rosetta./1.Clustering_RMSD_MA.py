#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 16:43:02 2019
Giuseppina La Sala
@author: kfsp227
"""

import os
import scipy
import mdtraj as md
import numpy as np
import pylab as plt
from scipy.cluster.hierarchy import fcluster,dendrogram,linkage
from scipy.spatial.distance import pdist,squareform

##### import trajectory ####

traj = md.load_pdb('allpdbs_lt603.pdb')
#trajProt = md.load_pdb('allpdbs_lt930.pdb', atom_indices=traj.topology.select('name CA'))
trajProt = md.load_pdb('allpdbs_lt603.pdb', atom_indices=traj.topology.select('resid 1 to 153'))
trajProt[0].save_pdb("reference_CA.pdb")


#### if trajectory is in xtc format

#dir='/projects/cc/lasala_kfsp227/W#orkDir/plainMD/LFA-1/PRODUCTION_rep1/'
#traj = md.load(dir+'production.part0002.PBC.xtc', top=dir+'protMG.pdb', stride=10000)
#trajProt = md.load(dir+'production.part0002.PBC.xtc', top=dir+'protMG.pdb', stride=100, atom_indices=traj.topology.select('name CA'))
#trajProt[0].save_pdb("reference_CA.pdb")


#compute RMSD matrix
distances = np.empty((trajProt.n_frames, trajProt.n_frames))
for i in range(trajProt.n_frames):
    distances[i] = md.rmsd(trajProt, trajProt, i)
#print('Max pairwise rmsd: %f nm' % np.max(distances))
reduced_distances = squareform(distances, checks=False)

#compute normal RMSD
Rmsd=md.rmsd(trajProt,trajProt,0)

reduced_distances = squareform(distances, checks=False)

# Plot RMSD matrix
plt.matshow(squareform(reduced_distances))
plt.colorbar()
plt.title('distance matrix')
plt.xlabel('frames')
plt.ylabel('frames')
plt.savefig('RMSDmatrix.png', dpi=300)

# Do the cluster
Z=scipy.cluster.hierarchy.linkage(reduced_distances, method='average')
plt.figure()
K=dendrogram(Z)
sortIndex=K['leaves']
plt.savefig('dendrogramm.png', dpi=300)

# Elbow Method
plt.figure()
plt.title('elbow method')
plt.grid()
plt.xlim(0,20)
plt.plot(np.arange(len(Z[:,2][::-1]))+1,Z[:,2][::-1]) # Plot last column of the Z linkage-matrix (contains the clusters)
plt.savefig('elbow.png', dpi=300)

# Plot the Clusters agaist the RMSD plot
ClusterNum=2  #edit number of clusters

out=fcluster(Z,int(ClusterNum),criterion='maxclust')

# associate frame and cluster
mapping={clusterId:np.where(out==clusterId)[0] for clusterId in set(out)}
centroids=[np.mean(Rmsd[mapping[clusterId]]) for clusterId in mapping.keys()]

#centroid sorting MA 3/16/2021
beta = 1
index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()

file = open('clusters.txt','w')
file.write("CLUSTER ID\tPOPULATION\tCENTROID\tMEMBERS\n")
for clusterId in set(out):
    file.write(str(clusterId) +'\t\t'+ str(len(mapping[clusterId])) + '\t\t' + '\t\t' +  str(mapping[clusterId]) + '\n')
file.close()

currPath=os.getcwd()
for c in np.arange(ClusterNum)+1:

    os.mkdir(os.getcwd()+"/CLUSTER"+str(c))
    os.chdir(os.getcwd()+"/CLUSTER"+str(c))
    for member in mapping[c]:
        trajProt[member].save_pdb("cluster"+str(c)+"_"+str(member)+".pdb")
    os.chdir(currPath)

