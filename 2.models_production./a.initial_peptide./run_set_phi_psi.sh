#!/bin/bash

/opt/scp/software/Rosetta/3.10-foss-2019a/bin/rosetta_scripts.mpi.linuxgccrelease \
 -database /opt/scp/software/Rosetta/3.10-foss-2019a/database/ \
 -in::file::s $1 \
 -parser::protocol set_phi_psi.xml \
 -ignore_unrecognized_res \
 -nstruct 10 \
 -out::suffix _phipsi\
 -overwrite
