import os, sys

data=open("score_docking_unfolded_weights.sc", "r")

lines=data.readlines()[1:]

for line in lines:
	
	fa_atr=line.strip()[15]

	print(fa_atr)
