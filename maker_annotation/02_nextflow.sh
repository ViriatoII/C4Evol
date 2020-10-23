#PBS -l select=1:ncpus=12:mem=16G
#PBS -l walltime=34:00:00
#PBS -A "C4Evol"


 module load Augustus/3.3.2 BLAST+/2.2.29 Maker/2.31.8 gcc/8.1.0 Miniconda/3 
 module load /software/modules/compiler/intelmpi/5.0.2.044


# prepare environment
sed 's/#module load /module load/' ~/.bashrc > tmp ; mv tmp ~/.bashrc   # uncomment module load miniconda3 from .bashrc
sed 's/#source/ source/' ~/.profile > tmp ; mv tmp ~/.profile 

# Prefix for ESTs run or db (v2) run 
MPI="MPI_db" # _ESTs" 

if [[ $MPI =~ "ESTs" ]]; then out="nextflow_ESTs" 
elif [[ $MPI =~ "db" ]]; then out="nextflow_db"
else out="nextflow_abInitio" ; fi 


source ~/.bashrc
conda activate nextflow-env
 
cd $PBS_O_WORKDIR         # ENTER SCRIPT DIRECTORY


# Get species name from fasta
#sp=`ls *fasta.index  | sed 's/.fasta.index//'`
sp=`ls *fasta  | sed 's/.fasta//'`


## Separate GFFs of MPI run 
cd $MPI.maker.output

        /gpfs/project/projects/qggp/src/maker2/bin/gff3_merge  -d ${MPI}_master_datastore_index.log
        /gpfs/project/projects/qggp/src/maker2/bin/fasta_merge -d ${MPI}_master_datastore_index.log

        perl /gpfs/project/projects/qggp/src/maker2/bin/AED_cdf_generator.pl -b 0.025 *.all.gff > AED.txt

        # transcript alignments
#        awk '{ if ($2 ~ "est2genome") print $0 }'      $MPI.all.gff > MPI.est2genome.gff  
        # protein alignments
        awk '{ if ($2 ~ "protein2genome") print $0 }'   $MPI.all.gff > MPI.protein2genome.gff
        # repeat alignments
        awk '{ if ($2 ~ "repeat") print $0 }'           $MPI.all.gff > MPI.repeats.gff

cd .. ; wait



# filtered gff uses less RAM 
grep maker $MPI.maker.output/$MPI.all.gff > filtered.gff 


# fold assembly.fasta into multiliner with 50 chars ### -> otherwise nextflow error
 fold -80  $sp.fasta > $MPI.maker.output/multiliner.fasta ; wait

#################### Abinitio training SNAP and AUGUSTUS ################################

NXF_SCRIPT="/home/guerrer/projects/src/nextflow/AbinitioTraining/AbinitioTraining.nf"


nextflow run -c ../pbs.config -profile conda,pbs $NXF_SCRIPT \
    --genome $MPI.maker.output/multiliner.fasta \
    --maker_evidence_gff 'filtered.gff'  \
    --outdir "$out" --species_label $sp

wait


# configure augustus to take nexflow model as input
rm -r ~/projects/src/augustus/config/species/$sp
cp -r  $out/Augustus_training/$sp ~/projects/src/augustus/config/species/.   ; rm filtered.gff


export AUGUSTUS_CONFIG_PATH="/gpfs/project/projects/qggp/src/augustus/config/"  # uses our custom config folders


################### Replace FASTAS by aligned GFFs. BUG??? ############################

# replace parameters in maker_opts: comment out fastas, feed aligned gffs
# Give created models to  maker_opts.ctl
sed    "s/augustus_spe.*/augustus_species=$sp/
	s/^snaphmm.*/snaphmm=$out\/Snap_training\/$sp.hmm/ 
	s/genome=1/genome=0/
	s/always_complete=0/always_complete=1/
	s/protein=\//protein= #\// 
        s/^est=\//est= #\// 
        s/rmlib=/rmlib= #\//
        s/rm_gff=.*/rm_gff=$MPI.maker.output\/MPI.repeats.gff/ 
        s/protein_gff=.*/protein_gff=$MPI.maker.output\/MPI.protein2genome.gff/" maker_opts.ctl > tmp ; mv tmp maker_opts.ctl
	# s/^est_gff=.*/est_gff=$MPI.maker.output\/MPI.est2genome.gff/   # INCLUDE THIS IF MPI="MPI_ESTs"

#######################################################################################

if [[ $MPI =~ "ESTs" ]]; then out="nextflow_ESTs" 
elif [[ $MPI =~ "db" ]]; then out="nextflow_db" 
else out="nextflow" ; fi

# Run 2nd round of maker with nextflow models
mpirun -n 12 maker2 -base $out -fix_nucleotides ; wait


cd $out.maker.output  

	module load Python/3.6.5

	/gpfs/project/projects/qggp/src/maker2/bin/gff3_merge -d *index.log 
	/gpfs/project/projects/qggp/src/maker2/bin/fasta_merge -d *index.log
	
	perl /gpfs/project/projects/qggp/src/maker2/bin/AED_cdf_generator.pl -b 0.025 *.all.gff > AED.txt

	# FILTER BY AEDs and CHANGE GENE NAME 
	sh ../../03_filter_AEDs_change_gen_name.sh

sed 's/.*odule load /#module load /' ~/.bashrc > tmp ; mv tmp ~/.bashrc    # comment module load out again
sed 's/ source/#source/' ~/.profile > tmp ; mv tmp ~/.profile  

