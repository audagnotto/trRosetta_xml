import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


x, y = [], []
x=list(range(2, 153))
with open("MD_RMSF.xvg") as f:
	for line in f:
		cols=line.split()

		if len(cols) == 2:
			y.append(float(cols[1]))

#xm2, ym2 = [], []
#xm2=list(range(1, 154))
#with open("../2jho/MD_RMSF.xvg") as q:
#        for line in q:
#                cols=line.split()
#
#                if len(cols) == 2:
#                        ym2.append(float(cols[1]))


x2, y2 = [], []
x2=list(range(2, 153))
with open("CLUSTER2_RMSF.xvg") as n:
        for line in n:
                cols=line.split()

                if len(cols) == 2:
                        y2.append(float(cols[1]))


#fig = plt.figure()
#ax1=fig.add_subplot(111)
#ax1.plot(x,y, 'bo')plt.hold(True)
#plt.hold(True)
fig= plt.figure(figsize=(6,2))
plt.gca().add_patch(Rectangle((40,0),15,0.5,fill=True, color='teal', alpha=0.2 ,zorder=100, figure=fig))
plt.gca().add_patch(Rectangle((120,0),15,0.5,fill=True, color='purple', alpha=0.2 ,zorder=100, figure=fig))
plt.plot(x2,y2, 'darkorange')
#plt.plot(x2,y2, 'forestgreen')
plt.plot(x,y, 'dimgray')
#plt.plot(xm2,ym2, 'black')
plt.xlabel('residue number')
plt.ylabel('RMSF [nm]')
plt.savefig('RMSF_PC1_gmx.png', dpi=800)
