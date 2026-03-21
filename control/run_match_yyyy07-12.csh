#!/bin/csh
#$ -l arch=lx-amd64
#$ -q scf.q
source $HOME/.cshrc_ceres
/SPG_ops/utils/x86_64/mpich2-1.4/bin/mpirun -np 8
/homedir/dfillmor/mpi/mpitask /homedir/dfillmor/MATCH/runs/$yyyy/runlist_$yyyy07-12
