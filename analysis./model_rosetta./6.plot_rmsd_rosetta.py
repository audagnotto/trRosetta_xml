import os, sys
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from scipy import stats

rmsd_file=open("rmsd_averaged_T1034.txt","r")
energy_file=open("../../2.models_production/c.3dmodel-cst/b.glytoala/score_names_list.txt","r")


lines_rmsd=rmsd_file.readlines()
lines_energy=energy_file.readlines()

outfile=open("rmsd_energy_rosetta.txt","w")
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
#regression curve
x=np.array(order_rmsd)
y=np.array(order_energy)
res = stats.linregress(x,y)
R2=(res.rvalue**2)

#plot
plt.scatter(x, y, label='original data')
plt.plot(x, res.intercept + res.slope *x, 'r', label='fitted line')
plt.annotate("R-sqrt:%.2f"%(round(R2,2)),xy=(25,1800))
plt.title("RMSD vs Energy")
plt.xlabel("RMSD(A)")
plt.ylabel("Energy (REU)")
plt.legend()
plt.savefig('rmsd_energy_rosetta.png', dpi=300)

#save data
for el in range(len(order_energy)):
	outfile.write("%s %s %s \n" %(order_rmsd[el],order_energy[el],order_names[el]))
