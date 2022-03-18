#!/bin/bash

/opt/scp/software/Rosetta/3.10-foss-2019a/bin/rosetta_scripts.mpi.linuxgccrelease \
 -database /opt/scp/software/Rosetta/3.10-foss-2019a/database/ \
 -in::file::s $1 \
 -parser::protocol mutateGLYtoALA_min_step0_relax.xml\
 -parser:script_vars mutateGA_resfile_relpath=GtoA.resfile mutateAG_resfile_relpath=AtoG.resfile \
 -run:repeat 3 \
 -ignore_unrecognized_res \
 -out::suffix $1_mincst_relax \
 -nstruct 100 \
 -optimization::default_max_cycles 200 \
 -overwrite
