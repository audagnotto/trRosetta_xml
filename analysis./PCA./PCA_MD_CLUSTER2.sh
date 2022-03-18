module load GROMACS/2020.3-fosscuda-2019a

#MD
echo q | gmx make_ndx -f renumbered_2_153.pdb -o index.ndx

echo 3 3 | gmx covar -s renumbered_2_153.pdb -f renumbered_2_153.trr -n index.ndx -xvg none -v MD_eigenvector.trr -o MD_eigenvalue.xvg -last 147

echo 3 3 | gmx anaeig -s renumbered_2_153.pdb -f renumbered_2_153.trr -2d MD_PC1_PC2.xvg -first 1 -last 2 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -n index.ndx 

#compute PC1
echo 3 3 | gmx anaeig -s renumbered_2_153.pdb -f renumbered_2_153.trr -first 1 -last 1 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -n index.ndx -proj MD_PC1.xvg

#RMSF from PC
echo 3 3 | gmx anaeig -s renumbered_2_153.pdb -f renumbered_2_153.trr -first 1 -last 1 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -n index.ndx -rmsf MD_RMSF_PC1.xvg

echo 3 3 | gmx anaeig -s renumbered_2_153.pdb -f renumbered_2_153.trr -first 2 -last 2 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -n index.ndx -rmsf MD_RMSF_PC2.xvg

#compute the RMSF
echo 3 | gmx rmsf -s renumbered_2_153.pdb -f renumbered_2_153.trr -n index.ndx -o MD_RMSF.xvg -xvg none

#models
echo q | gmx make_ndx -f ./CLUSTER2/initial.pdb -o ./CLUSTER2/index.ndx

echo 3 3 | gmx covar -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -xvg none -v CLUSTER2_eigenvector.trr -o CLUSTER2_eigenvalue.xvg -last 147

echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -2d cluster1_PC1_PC2.xvg -first 1 -last 2 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg

#compute PC1
echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -first 1 -last 1 -xvg none -v CLUSTER2_eigenvector.trr -eig CLUSTER2_eigenvalue.xvg -proj CLUSTER2_PC1.xvg

echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -first 1 -last 1 -xvg none -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -proj CLUSTER2_PC1MD.xvg


#RMSF from PC
echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -first 1 -last 1 -xvg none -v CLUSTER2_eigenvector.trr -eig CLUSTER2_eigenvalue.xvg -rmsf CLUSTER2_RMSF_PC1.xvg

echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -first 2 -last 2 -xvg none -v CLUSTER2_eigenvector.trr -eig CLUSTER2_eigenvalue.xvg -rmsf CLUSTER2_RMSF_PC2.xvg

#compute the RMSF
echo 3 | gmx rmsf -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -xvg none -o CLUSTER2_RMSF.xvg

#save the PC1 to smoothy visualize
echo 3 3 | gmx anaeig -s renumbered_2_153.pdb -f renumbered_2_153.trr -first 1 -last 1 -v MD_eigenvector.trr -eig MD_eigenvalue.xvg -n index.ndx -extr MDextractionPC1.pdb -nframes 60


echo 3 3 | gmx anaeig -s ./CLUSTER2/initial.pdb -f ./CLUSTER2/cluster2.trr -n ./CLUSTER2/index.ndx -first 1 -last 1 -xvg none -v CLUSTER2_eigenvector.trr -eig CLUSTER2_eigenvalue.xvg -extr CLUSTER2sextractionPC1.pdb -nframes 60

#compute the PCA internal product
gmx anaeig -v MD_eigenvector.trr -v2 CLUSTER2_eigenvector.trr -inpr MD_CLUSTER2_RMSP.xpm -first 1 -last 6

python convertXPM2DAT.py MD_CLUSTER2_RMSP.xpm > MD_CLUSTER2_RMSP.data
