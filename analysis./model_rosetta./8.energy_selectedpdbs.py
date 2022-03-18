import os, sys

def lines_that_contain(string,fp):
	return [line for line in fp if string in line]

file1=open("namemodel-frame.txt", "r")

list1=file1.readlines()

energy_line=[]
for line in list1:
	name=(line.strip().split(" ")[1])
	with open("rescore_AWESM-MD/awsem-energy.txt", "r") as fp:
		for line in lines_that_contain(name, fp):
			energy_line.append(line)

output=open("energy_namelt859.txt", "w")
for el in energy_line:
	output.write("%s" %el)
	
