import os, sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import seaborn as sb

#data
dataMD=open("MD_PC1.xvg", "r")
lines=dataMD.readlines()

pcMD=[]
for line in lines:
	cols=line.split()
	if len(cols) == 2:
		pcMD.append(float(line.strip().split()[1]))
densityMD = gaussian_kde(pcMD)

datacl3=open("CLUSTER2_PC1MD.xvg", "r")
lines=datacl3.readlines()

pccl3=[]
for line in lines:
        cols=line.split()
        if len(cols) == 2:
                pccl3.append(float(line.strip().split()[1]))
densitycl3 = gaussian_kde(pccl3)

#dataMD2=open("../comparison_2eb8_2jho/2eb8_2jhoPC1.xvg", "r")
#lines=dataMD2.readlines()

#pcMD2=[]
#for line in lines:
#        cols=line.split()
#        if len(cols) == 2:
#                pcMD2.append(float(line.strip().split()[1]))
#densityMD2 = gaussian_kde(pcMD2)

#plot

#plt.plot(xs,densityMD(xs), c="grey")
#plt.plot(xs,densitycl1(xs), c='darkorange')
#plt.plot(xs,densitycl3(xs), c='forestgreen')
fig= plt.figure(figsize=(6,6))
sb.kdeplot(pccl3 , bw = 0.15 , color='darkorange', fill=True)
sb.kdeplot(pcMD , bw = 0.15 , color='slategrey', fill=True)
#sb.kdeplot(pcMD2, bw = 0.15 , color='black', fill=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('PC1', fontsize=16, fontweight='normal')
plt.ylabel('Density', fontsize=16, fontweight='normal')
plt.savefig("density_PC1.png", dpi=800)
