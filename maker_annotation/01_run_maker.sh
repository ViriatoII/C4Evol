#!/bin/bash

#PBS -l select=1:ncpus=12:mem=40G
#PBS -l walltime=59:00:00
#PBS -A "C4Evol"

module load Augustus/3.3.2 BLAST+/2.2.29 Maker/2.31.8 gcc/8.1.0 
module load /software/modules/compiler/intelmpi/5.0.2.044
module unload Perl 

alias maker2="/gpfs/project/projects/qggp/src/maker2/bin/maker"
source ~/.bashrc 

cd $PBS_O_WORKDIR         # ENTER SCRIPT DIRECTORY

# 2 possible rounds: 	-ESTs for initial annotation using transcriptomes and prot_db_v1
#			-db for final annotation using only prot_db_v2  (made from prot_db_v1 and final ESTs annotations)
MPI="MPI_ESTs" # "_db"


# Run maker  # NOTE: Needs .ctl files
mpirun -n 12  maker2 -base $MPI -fix_nucleotides  ; wait 


cd $MPI.maker.output ; wait

# Generate outputs 
/gpfs/project/projects/qggp/src/maker2/bin/gff3_merge  -d ${MPI}_master_datastore_index.log ; wait 
/gpfs/project/projects/qggp/src/maker2/bin/fasta_merge -d ${MPI}_master_datastore_index.log ; wait

# fold -b70 ../*.fasta > multiliner.fasta

perl /gpfs/project/projects/qggp/src/maker2/bin/AED_cdf_generator.pl -b 0.025 $MPI.all.gff > AED.txt 
