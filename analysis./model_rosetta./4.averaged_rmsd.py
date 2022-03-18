import os, sys
#average function
def Average(lst):
	return sum(lst)/len(lst)

data=open("RMSD_pairframes.txt", "r")
lines=data.readlines()

name=[]
for line in lines:
	name.append(line.strip().split(" ")[3])

newnames=sorted(set(name))
average=[]
output_average=open("rmsd_averaged_T1034.txt", "w")
for newname in newnames:
	rmsd_frames=[]
	for line in lines:
		name=(line.strip().split(" ")[3])
		rmsd=float(line.strip().split(" ")[2])
		if name == newname:
			rmsd_frames.append(rmsd)
	average.append(Average(rmsd_frames))
		
for el in range(len(newnames)):
	output_average.write("%f %s \n" %(average[el], newnames[el].split("/")[-1][0:-4]))
