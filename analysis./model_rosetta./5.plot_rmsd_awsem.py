import os, sys
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import numpy as np
import statistics as st

rmsd_file=open("rmsd_averaged_T1034.txt","r")
energy_file=open("rescore_AWESM-MD/awsem-energy.txt","r")


lines_rmsd=rmsd_file.readlines()
lines_energy=energy_file.readlines()

outfile=open("rmsd_energy_awsem.txt","w")
order_rmsd=[]
order_energy=[]
order_names=[]
for line_rmsd in lines_rmsd:
	models_rmsd=(line_rmsd.strip().split(" ")[1])
	for line_energy in lines_energy:
		models_en=(line_energy.strip().split(" ")[1])
		if models_en == models_rmsd:
			order_rmsd.append(float(line_rmsd.strip().split(" ")[0]))
			order_energy.append(float(line_energy.strip().split(" ")[0]))
			order_names.append(models_en)

#mean and SD
mean=st.mean(order_energy)
st=st.stdev(order_energy)
cutoff=-(abs(mean)+st)
print(mean)
print(st)

#write out the list of pdbs that have en < -930
output_en=open("list_pdbs_lt%d.txt"%(round(cutoff)), "w")
for en,name in zip(order_energy, order_names):
        if en < cutoff:
                output_en.write("%s\n" %name)

#regression curve
x=np.array(order_rmsd)
y=np.array(order_energy)
res = stats.linregress(x,y)
#print(f"R-squared: {res.rvalue**2:.6f}")
R2=(res.rvalue**2)

#plot
plt.scatter(x, y, label='original data')
plt.plot(x, res.intercept + res.slope *x, 'r', label='fitted line')
plt.axhline(y=cutoff, color='gray', linestyle='--')
plt.annotate("R-sqrt:%.2f"%(round(R2,2)),xy=(25,-600))
plt.title("RMSD vs Energy")
plt.xlabel("RMSD-T1034 (A)")
plt.ylabel("Energy (kcal/mol)")
plt.legend()
plt.savefig('rmsd_energy_awesm.png', dpi=300)

#save data
for el in range(len(order_energy)):
	outfile.write("%s %s %s \n" %(order_rmsd[el],order_energy[el],order_names[el]))
