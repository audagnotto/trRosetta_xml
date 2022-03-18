############################################################################################
# python script to:
# parse the score files and get the models and the coarse score
#source activate audagnotto_kblr332v2
############################################################################################
import os
import sys
import glob
import re
import numpy as np
import shutil
import natsort

#take all the coarse score function
onlyfiles=glob.glob("./step_0_score_generation/*.sc")

score_all=[]
for fname in onlyfiles:
#for fname in name_score:
    print("working on %s" % fname)
    with open(fname) as infile:
        lines_skipfirst=infile.readlines()[2:]
        score_all.append(lines_skipfirst)

#open the file to write the score models list
rmsd_score=("score_names_list.txt")
list_files=("list_files.txt")
#sort the model in ascending order
myfiles=[]
dirFiles = os.listdir('./step_0_model_generation/')
myfiles.append(natsort.natsorted(dirFiles))
#myfiles.append(sorted(dirFiles))
#print(myfiles)
with open(rmsd_score,'w') as rmsdscore_out, open(list_files, "w") as list_out:
    for i in range(len(score_all)):
	#j is the number of models per run
        for j in range(100):
            #'name of the model'
            name_model=(((" ".join((score_all[i][j].split()))).split(" "))[-1])
            name_tofind=name_model+"_0001"
            #coarse_total_score
            score_coarse=float((((" ".join((score_all[i][j].split()))).split(" "))[1]))
            for el in myfiles:
		#k is the total number of files
                for k in range(1000):
                    name=el[k][0:-4]
                    #print(name)
                    if name==name_model:

                        rmsdscore_out.write("%f %s \n" % (score_coarse, name))
                        list_out.write("./step_0_model_generation/%s.pdb \n" % (name))
