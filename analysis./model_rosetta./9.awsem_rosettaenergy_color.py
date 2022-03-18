import os, sys
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import numpy as np
import statistics as st

awsem=open("rmsd_energy_awsem.txt", "r")
rosetta=open("rmsd_energy_rosetta.txt", "r")
exp=open("../../experimental/rescore_awsem/awsem-energy.txt", "r")

awsem_lines=awsem.readlines()
rosetta_lines=rosetta.readlines()
exp_lines=exp.readlines()

rmsd=[]
energy_awsem=[]
for awsem_line in awsem_lines:
	rmsd.append(float(awsem_line.strip().split()[0]))
	energy_awsem.append(float(awsem_line.strip().split()[1]))

energy_rosetta=[]
for rosetta_line in rosetta_lines:
	energy_rosetta.append(float(rosetta_line.strip().split()[1]))

exp_value=[]
for exp_line in exp_lines:
	exp_value.append(float(exp_line.strip().split()[0]))

diff=[float(abs(exp_value[0]-exp_value[1]))]

x=np.array(rmsd)
y=np.array(energy_awsem)
z=np.array(energy_rosetta)

#mean and SD
mean=st.mean(energy_awsem)
st=st.stdev(energy_awsem)
cutoff=-(abs(mean)+st)
print(mean)
print(st)

#plot
plt.hlines(max(exp_value),min(x),max(x), linestyles='dashed', color='black')
plt.hlines(min(exp_value),min(x),max(x), linestyles='dashed', color='black')
print(max(exp_value))
fig=plt.figure(constrained_layout=True)
plt.axhline(y=mean, color='gray', linestyle='--')
plt.axhline(y=max(exp_value), color='black', linestyle='--')
plt.axhline(y=min(exp_value), color='black', linestyle='--')

cm = plt.cm.get_cmap('RdYlBu')
sc=plt.scatter(x, y, c=z, cmap=cm)
#plt.title("Rosetta vs Awsem")
plt.ylabel('Energy [kcal/mol]', fontsize=16, fontweight='normal')
fig.gca().set_xlabel(r'RMSD [$\AA$]', fontsize=16, fontweight='normal')
cbar=plt.colorbar(sc)
cbar.set_label('Rosetta Energy (REU)',fontsize=16, fontweight='normal')
plt.savefig('rmsd_energy_awesm_rosetta.png', dpi=300)

