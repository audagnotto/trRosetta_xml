import numpy as np
import matplotlib.pyplot as plt


x, y = [], []
with open("./6y75/6Y75_distance_angle.txt") as f:
	for line in f:
		cols=line.split()
		if len(cols) == 3:
			x.append(float(cols[1]))
			y.append(float(cols[0]))


x2, y2 = [], []
with open("./CLUSTER2/cluster2_distance_angle.txt") as m:
        for line in m:
                cols=line.split()
                if len(cols) == 3:
                        x2.append(float(cols[1]))
                        y2.append(float(cols[0]))

xexp1, yexp1 = [], []
with open("6y75_exp/6Y75_exp_distance_angle.txt") as q:
        for line in q:
                cols=line.split()
                if len(cols) == 3:
                        xexp1.append(float(cols[1]))
                        yexp1.append(float(cols[0]))
xexp2, yexp2 = [], []
with open("6tmm_exp/6TMM_exp_distance_angle.txt") as q:
        for line in q:
                cols=line.split()
                if len(cols) == 3:
                        xexp2.append(float(cols[1]))
                        yexp2.append(float(cols[0]))

xaf, yaf = [], []
with open("AF2/AF_distance_angle.txt") as r:
        for line in r:
                cols=line.split()
                if len(cols) == 3:
                        xaf.append(float(cols[1]))
                        yaf.append(float(cols[0]))


xcol, ycol = [], []
with open("T1034_colalab_AF2/colab_AF2_distance_angle.txt") as s:
        for line in s:
                cols=line.split()
                if len(cols) == 3:
                        xcol.append(float(cols[1]))
                        ycol.append(float(cols[0]))


fig= plt.figure(figsize=(6,6))
#plt.scatter(xm,ym,c='lightslategrey', alpha=0.2)
plt.scatter(x,y,c='dimgray', alpha=0.2)
plt.scatter(x2,y2, c='darkorange', alpha=0.5)
plt.scatter(xcol,ycol,c='firebrick',s=130, alpha=0.7)
plt.scatter(xexp1, yexp1, c='black', alpha=0.7)
plt.scatter(xexp2, yexp2, c='lightslategrey', alpha=0.7)
plt.scatter(xaf,yaf,c='firebrick',s=130, marker=(5, 1), edgecolors= "black", alpha=0.7)
plt.xlabel('angle [degree]', fontsize=16, fontweight='normal')
fig.gca().set_ylabel(r'distance [$\AA$]', fontsize=16, fontweight='normal')
plt.savefig('T1034_distance_angle_MD-models.png', dpi=800)
			
