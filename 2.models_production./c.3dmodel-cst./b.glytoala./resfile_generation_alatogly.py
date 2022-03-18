import os, sys
from Bio.PDB import *


p = PDBParser()
structure=p.get_structure('model1','%s' %(sys.argv[1]))

model=structure[0]
chain = model['A']


resgly=[]

for r in structure.get_residues():
    if (r.get_resname())=='GLY':
        resgly.append(r.id[1])

output=open('AtoG.resfile',"w")
output.write("NATRO \n")
output.write("start \n")
for re in resgly:
    #resid chain PIKAA mutation
    output.write("%d A PIKAA G \n" %re)

output.close()
