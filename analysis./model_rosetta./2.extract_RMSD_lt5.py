import os, sys
import numpy as np
import pylab as plt
import shutil

#read the input file
input_file=open ("RMSD_pairframes.txt", "r")
lines =input_file.readlines()

output=open("list_pdbs_lt5_6y75.txt", "w")
output_rmsd=open("rmsd_lt5.txt", "w")
data=[]
rmsd_data=[]
for line in lines:
	rows=(int((line.strip().split(" ")[0])))
	columns=(int((line.strip().split(" ")[1])))
	rmsd=(float((line.strip().split(" ")[2])))
	name=((line.strip().split(" ")[3]))
	if rmsd < 5:
		data.append(name)	
		rmsd_data.append(rmsd)

data=sorted(set(data))

for rmsds in rmsd_data:
	output_rmsd.write("%f \n" %rmsds)

for el in data:
	output.write("%s \n" %el)
