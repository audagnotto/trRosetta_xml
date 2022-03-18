import os, sys
import numpy as np 
import re

file1=open("clusters.txt", "r")
file2=open("energy_namelt859.txt", "r")

line1=file1.readlines()[1:]
line2=file2.readlines()

indexes=[]
for line in line1:
	#print(re.sub("\s+", ",", line.strip()))
	index=(line.strip().split("\t")[6].strip("\[").strip("\]"))
	print(index)
	indexes.append((re.sub("\s+", ",", index.strip())))

output=open("centroids.txt", "w")
for i in range(len(indexes)):
	new_indexes=(indexes[i].split(","))
	energy=[]
	for el in new_indexes:
		el=int(el)
		energy.append(float(line2[el].split(" ")[0]))
	
	min_en=(min(energy))
	for line in line2:
		value=float(line.strip().split(" ")[0])
		if value == min_en:
			output.write("%f %s \n"%(value, line.strip().split(" ")[1]))

