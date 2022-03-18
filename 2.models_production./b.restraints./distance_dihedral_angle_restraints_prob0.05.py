##############################################################################################
#python script to convert the distance calculated from the distogram in distance restraints
#based on the utils_ros.py in trRosetta
###############################################################################################
import numpy as np
import sys, os
import random


def read_fasta(file):
    fasta=""
    with open(file, "r") as f:
        for line in f:
            if(line[0] == ">"):
                continue
            else:
                line=line.rstrip()
                fasta = fasta + line;
    return fasta

#define the constants
TDIR='./distance_restraints'
ASTEP=15.0
DSTEP=0.5
DCUT=19.5
ALPHA=1.57
MEFF=0.0001
EBASE=-0.5
EREP=[10.0,3.0,0.5]
DREP=[0.0,2.0,3.5]
PCUT=0.05
seq=read_fasta('%s'%(sys.argv[1]))

#open the pickle file
npz=np.load('%s'%(sys.argv[2]))
dist=npz['dist']
omega=npz['omega']
theta=npz['theta']
phi=npz['phi']

# dictionary to store Rosetta restraints
rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : []}

# define the name of the directory to be created
path = "./restraints"

# define the access rights
access_rights = 0o755

try:
    os.mkdir(path, access_rights)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s" % path)



#data stored in the npz file
########################################################
# dist: 0..20A                                         #
########################################################
nres = dist.shape[0]
#the distance conside are from 4.25 to 19.75 every 0.5 for a total of 32 bins
bins = np.array([4.25+DSTEP*i for i in range(32)])
#only 32/37 bins for the distance are considered based on the bins defintion before
# sum of all the probabilities across the column (axis=1). If axis is negative it counts from the last to the first axis
prob = np.sum(dist[:,:,5:], axis=-1)
#bkgr is the (di/dn)^alpha part of the score: di: dist at the ith bin and dn:dist at the
#Nth bin== DCUT (19.5))
bkgr = np.array((bins/DCUT)**ALPHA)
#MA 4/20/2020 score function definition: translation from probability value to score
#-np.log((dist[:,:,5:]+MEFF)-> pi [MEFF=0.0001 (why it is needed?)]
#dist[:,:,-1][:,:,None]=pN-> probability distance for the last bin @19.5
#dist[:,:,-1][:,:,None]*bkgr[None,None,:]) -> pN*(di/dN)^alpha [EBASE=-0.5 (?)]
attr = -np.log((dist[:,:,5:]+MEFF)/(dist[:,:,-1][:,:,None]*bkgr[None,None,:]))+EBASE
#np.maximum:Compare two arrays and returns a new array containing the element-wise maxima
#1st array: the score for the distance only the first 32 bins array (?)
#2nd array is zero array of 80 elements
repul = np.maximum(attr[:,:,0],np.zeros((nres,nres)))[:,:,None]+np.array(EREP)[None,None,:]
dist = np.concatenate([repul,attr], axis=-1)
#add the missing bins (from 0 to 3.5) and concatenate along the column
bins = np.concatenate([DREP,bins])
#return in an array all the positions i and j for which the probabilities > 0.05
i,j = np.where(prob>PCUT)
prob = prob[i,j]
nbins = 35
step = 0.5
for a,b,p in zip(i,j,prob):
    if b>a:
        name=path+"/%d.%d_distance.txt"%(a+1,b+1)
        with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(dist[a,b])+'\n')
                f.close()
        rst_line = 'AtomPair %s %d %s %d SPLINE TAG %s 1.0 %.3f %.5f'%('CB',a+1,'CB',b+1,name,1.0,step)
        rst['dist'].append([a,b,p,rst_line])

#total number of restraits
print("dist restraints:  %d"%(len(rst['dist'])))

########################################################
# omega: -pi..pi
########################################################
nbins = omega.shape[2]-1+4
bins = np.linspace(-np.pi-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
prob = np.sum(omega[:,:,1:], axis=-1)
i,j = np.where(prob>PCUT)
prob = prob[i,j]
omega = -np.log((omega+MEFF)/(omega[:,:,-1]+MEFF)[:,:,None])
omega = np.concatenate([omega[:,:,-2:],omega[:,:,1:],omega[:,:,1:3]],axis=-1)
for a,b,p in zip(i,j,prob):
    if b>a:
        name=path+"/%d.%d_omega.txt"%(a+1,b+1)
        with open(name, "w") as f:
                f.write('x_axis'+'\t%.5f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.5f'*nbins%tuple(omega[a,b])+'\n')
                f.close()
        rst_line = 'Dihedral CA %d CB %d CB %d CA %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,b+1,name,1.0,ASTEP)

        rst['omega'].append([a,b,p,rst_line])
print("omega restraints: %d"%(len(rst['omega'])))

########################################################
# theta: -pi..pi
########################################################
prob = np.sum(theta[:,:,1:], axis=-1)
i,j = np.where(prob>PCUT)
prob = prob[i,j]
theta = -np.log((theta+MEFF)/(theta[:,:,-1]+MEFF)[:,:,None])
theta = np.concatenate([theta[:,:,-2:],theta[:,:,1:],theta[:,:,1:3]],axis=-1)
for a,b,p in zip(i,j,prob):
    if b!=a:
        name=path+"/%d.%d_theta.txt"%(a+1,b+1)
        with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(theta[a,b])+'\n')
                f.close()
        rst_line = 'Dihedral N %d CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,a+1,b+1,name,1.0,ASTEP)
        rst['theta'].append([a,b,p,rst_line])

print("theta restraints: %d"%(len(rst['theta'])))

########################################################
# phi: 0..pi
########################################################
nbins = phi.shape[2]-1+4
bins = np.linspace(-1.5*ASTEP, np.pi+1.5*ASTEP, nbins)
prob = np.sum(phi[:,:,1:], axis=-1)
i,j = np.where(prob>PCUT)
prob = prob[i,j]
phi = -np.log((phi+MEFF)/(phi[:,:,-1]+MEFF)[:,:,None])
phi = np.concatenate([np.flip(phi[:,:,1:3],axis=-1),phi[:,:,1:],np.flip(phi[:,:,-2:],axis=-1)], axis=-1)
for a,b,p in zip(i,j,prob):
    if b!=a:
        name=path+"/%d.%d_phi.txt"%(a+1,b+1)
        with open(name, "w") as f:
                f.write('x_axis'+'\t%.3f'*nbins%tuple(bins)+'\n')
                f.write('y_axis'+'\t%.3f'*nbins%tuple(phi[a,b])+'\n')
                f.close()
        rst_line = 'Angle CA %d CB %d CB %d SPLINE TAG %s 1.0 %.3f %.5f'%(a+1,a+1,b+1,name,1.0,ASTEP)
        rst['phi'].append([a,b,p,rst_line])

print("phi restraints:   %d"%(len(rst['phi'])))

#RESTRAINTS BASED ON SEQUENCE SEPARATION!!
############################################
#SHORT RESTRAINTS: 1 12 SEQUENCE SEPARATION
############################################
array=[]
sep1=1
sep2=12
pcut=PCUT
array+= [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6]

random.shuffle(array)
minimize_constraints112 = path+'/minimize_1_12.cst'
with open(minimize_constraints112,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()

############################################
#MEDIUM RESTRAINTS: 12 24 SEQUENCE SEPARATION
############################################
array=[]
sep1=12
sep2=24
pcut=PCUT
array+= [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6]

random.shuffle(array)
minimize_constraints112 = path+'/minimize_12_24.cst'
with open(minimize_constraints112,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()

############################################
#LONG RESTRAINTS: 24 to length of the SEQUENCE
############################################
array=[]
sep1=24
sep2=len(seq)
pcut=PCUT
array+= [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6]

random.shuffle(array)
minimize_constraints112 = path+'/minimize_24.cst'
with open(minimize_constraints112,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()

#################################################
#SHORT+MEDIUM RESTRAINTS: 3 to 24 length of the sequence
#################################################
array=[]
sep1=3
sep2=24
pcut=PCUT
array+= [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6]

random.shuffle(array)
minimize_constraints112 = path+'/minimize_3_24.cst'
with open(minimize_constraints112,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()

#################################################
#SHORT+MEDIUM+LONG RESTRAINTS: 1 to length of the sequence
#################################################
array=[]
sep1=1
sep2=len(seq)
pcut=PCUT
array+= [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut]
array += [line for a,b,p,line in rst['omega'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['theta'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.5]
array += [line for a,b,p,line in rst['phi'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and p>=pcut+0.6]

random.shuffle(array)
minimize_constraints112 = path+'/minimize_1_allseq.cst'
with open(minimize_constraints112,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()
#######################################################
#FAST RELAX RESTRAINTS all the seq and probability 0.15
#######################################################
array=[]
sep1=1
sep2=len(seq)
pcut=0.15
array += [line for a,b,p,line in rst['dist'] if abs(a-b)>=sep1 and abs(a-b)<sep2 and seq[a]!='G' and seq[b]!='G' and p>=pcut]
random.shuffle(array)
minimize_constraints_relax = path+'/minimize_relax.cst'
with open(minimize_constraints_relax,'w') as f:
    for line in array:
        f.write(line+'\n')
    f.close()
