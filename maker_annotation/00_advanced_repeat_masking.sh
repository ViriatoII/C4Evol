#!/bin/bash

#PBS -l select=1:ncpus=12:mem=16G
#PBS -l walltime=59:00:00
#PBS -A "C4Evol"


DIR="/gpfs/project/projects/qggp/src/genometools-1.5.9/bin"
DIR2="/gpfs/project/projects/qggp/src/maker/CRL_Scripts1.0"

module load Augustus/3.3.2 BLAST+/2.2.29 Maker/2.31.8 Muscle gcc/8.1.0 

# Enter script directory. Should be a subfolder with the species name.
cd $PBS_O_WORKDIR          

GENOME="/DIRECTORY/ASSEMBLY.fasta"     # INPUT GENOME
SPECIES=`pwd | sed 's/.*\///'`         # GETS NAME FROM FOLDER



##### Remove aditional information from fasta header:   # NOTE: FASTA file must not have "_" or ".".	
#	sed 's/\(>.*\) edge.*/\1/' $GENOME > $SPECIES\.fasta  # for a 10x assembly
	sed 's/,.*//' $GENOME > $SPECIES\.fasta               # For ARCs scaffolds
	

mkdir repeat-masking ; cd repeat-masking

##################### 1) MITES ###################################################
# Takes a long time

mkdir mites ; cd mites # Detect MITEs repeats
 
/gpfs/project/projects/qggp/src/mite_hunter/MITE_Hunter_manager.pl -i ../../$SPECIES.fasta -g $SPECIES -n 12 -c 12 -S 12345678 ; wait 

cat  *Step8_*.fa > MITE.lib # the tutorial is deceiving on this line, singlet file is already considered by the *, no need to specify.

cd ..

#################### 2) LTR (long terminal repeat) retrotransposons) ####################

mkdir ltr ; cd ltr

ln -s ../../$SPECIES.fasta    # input for next section

$DIR/gt suffixerator -db $SPECIES.fasta -indexname $SPECIES\_index -tis -suf -lcp -des -ssp -dna

mkdir 99 85; cd 99   ; ln -s ../../../$SPECIES.fasta


# 2.1.1. harvest ltrs at 99 threshold

 $DIR/gt ltrharvest -index ../$SPECIES\_index -out $SPECIES.out99 \
 -outinner $SPECIES.outinner99 -gff3 $SPECIES.gff99 -minlenltr 100 \
 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 \
 -maxtsd 5 -motif tgca -similar 99 -vic 10  > $SPECIES.result99 ; wait

#  2.1.2. Find elements with PPT (poly purine tract) or PBS (primer binding site) 

 $DIR/gt gff3 -sort $SPECIES.gff99 > $SPECIES.gff99.sort ; wait

						  # eucaryotic RNAs downloaded from : http://lowelab.ucsc.edu/GtRNAdb/download.html
 $DIR/gt -j 20 ltrdigest -trnas ~/C34_PS/data/databases/tRNA/eukaryotic-tRNAs.fa $SPECIES.gff99.sort ../$SPECIES\_index > $SPECIES.gff99.dgt ; wait


#  2.1.3. Further filtering of the candidate elements 

 perl $DIR2/CRL_Step1.pl  --gff  $SPECIES.gff99.dgt ; wait

 perl $DIR2/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt  --repeatfile $SPECIES.out99 --resultfile $SPECIES.result99 \
 --sequencefile $SPECIES.fasta  --removed_repeats   CRL_Step2_Passed_Elements.fasta ; wait

mkdir fasta_files

mv Repeat_*.fasta fasta_files/. ; mv  CRL_Step2_Passed_Elements.fasta  fasta_files/. 

cd  fasta_files

perl $DIR2/CRL_Step3.pl --directory $PWD --step2 CRL_Step2_Passed_Elements.fasta  --pidentity 60 \
--seq_c 25 ; wait

mv  CRL_Step3_Passed_Elements.fasta  .. ; cd  ..

# 2.1.4. Identify elements with nested insertions 

perl $DIR2/ltr_library.pl --resultfile $SPECIES.result99  --step3 CRL_Step3_Passed_Elements.fasta --sequencefile $SPECIES.fasta ; wait


cat lLTR_Only.lib  ../../mites/MITE.lib  > repeats_to_mask_LTR99.fasta


/gpfs/project/projects/qggp/src/repeatmasker-4-0-9/RepeatMasker -pa 20 -lib repeats_to_mask_LTR99.fasta  \
  -nolow  -dir . $SPECIES.outinner99   ; wait   # E. COLI MESSAGE???

perl $DIR2/cleanRM.pl  $SPECIES.outinner99.out  $SPECIES.outinner99.masked  >  $SPECIES.outinner99.unmasked  ; wait


perl $DIR2/rmshortinner.pl  $SPECIES.outinner99.unmasked  50 > $SPECIES.outinner99.clean  ; wait


makeblastdb -in ~/C34_PS/data/databases/transposons/tpases020812DNA.fa  -dbtype prot
					          #“tpases020812DNA” represents the DNA transposase protein sequence file

blastx  -query  $SPECIES.outinner99.clean -db ~/C34_PS/data/databases/transposons/tpases020812DNA.fa \
 -evalue 1e-10 -num_descriptions 10 -num_threads 20 -out $SPECIES.outinner99.clean_blastx.out.txt


perl  $DIR2/outinner_blastx_parse.pl --blastx  $SPECIES.outinner99.clean_blastx.out.txt  --outinner  $SPECIES.outinner99

# 2.1.5  Building examplars (choosing representative repeats)

perl $DIR2/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile $SPECIES.result99 \
--innerfile passed_outinner_sequence.fasta  --sequencefile $SPECIES.fasta ; wait

makeblastdb -in lLTRs_Seq_For_BLAST.fasta  -dbtype nucl

blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000  \
-out lLTRs_Seq_For_BLAST.fasta.out

makeblastdb -in  Inner_Seq_For_BLAST.fasta  -dbtype nucl

blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta  -evalue 1e-10 -num_descriptions 1000 \
-out  Inner_Seq_For_BLAST.fasta.out

perl $DIR2/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out  \
--step3 CRL_Step3_Passed_Elements.fasta --final LTR99.lib --pcoverage 90 --pidentity 80

cd ../85/  ; ln -s ../../../$SPECIES.fasta

##### THIS PART IS EXPERIMENTAL ######################

 # it's a repetition of previous comands but with 85% threshold. Maybe colapse this with a loop?

# 2.2. Collection of relatively old LTR retrotransposons

# 2.2.1. Harvest ltrs at 85% threshold

  ln -s ../../$SPECIES.fasta    # input for next section

  $DIR/gt ltrharvest -index ../$SPECIES\_index -out $SPECIES.out85 \
  -outinner $SPECIES.outinner85 -gff3 $SPECIES.gff85 -minlenltr 100 \
  -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10  > $SPECIES.result85 ; wait


#  2.2.2. Find elements with PPT (poly purine tract) or PBS (primer binding site) 

  $DIR/gt gff3 -sort $SPECIES.gff85 > $SPECIES.gff85.sort ; wait


   $DIR/gt -j 20 ltrdigest -trnas ~/C34_PS/data/databases/tRNA/eukaryotic-tRNAs.fa $SPECIES.gff85.sort ../$SPECIES\_index > $SPECIES.gff85.dgt  ; wait
       
  
  perl $DIR2/CRL_Step1.pl  --gff  $SPECIES.gff85.dgt ; wait

#  2.2.3. Further filtering of the candidate elements 

  perl $DIR2/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt  --repeatfile $SPECIES.out85 --resultfile $SPECIES.result85 \
 --sequencefile $SPECIES.fasta  --removed_repeats   CRL_Step2_Passed_Elements.fasta ; wait


  mkdir fasta_files2

  mv Repeat_*.fasta fasta_files2/. ; mv  CRL_Step2_Passed_Elements.fasta  fasta_files2/.

  cd  fasta_files2

  perl $DIR2/CRL_Step3.pl --directory .  --step2 CRL_Step2_Passed_Elements.fasta  --pidentity 60 \
  --seq_c 25 ; wait

  mv  CRL_Step3_Passed_Elements.fasta  .. ; cd  ..

# 2.1.4. Identify elements with nested insertions 

  perl $DIR2/ltr_library.pl --resultfile $SPECIES.result85  --step3 CRL_Step3_Passed_Elements.fasta --sequencefile $SPECIES.fasta ; wait


  cat lLTR_Only.lib  ../../mites/MITE.lib  > repeats_to_mask_LTR85.fasta ; wait


  /gpfs/project/projects/qggp/src/repeatmasker-4-0-9/RepeatMasker -pa 10 -lib repeats_to_mask_LTR85.fasta  \
  -nolow -dir . $SPECIES.outinner85  ; wait    # E. COLI MESSAGE???


  perl $DIR2/cleanRM.pl  $SPECIES.outinner85.out  $SPECIES.outinner85.masked  >  $SPECIES.outinner85.unmasked


  perl $DIR2/rmshortinner.pl  $SPECIES.outinner85.unmasked  50 > $SPECIES.outinner85.clean


  makeblastdb -in ~/C34_PS/data/databases/transposons/tpases020812DNA.fa  -dbtype prot

  blastx  -query  $SPECIES.outinner85.clean -db ~/C34_PS/data/databases/transposons/tpases020812DNA.fa \
  -evalue 1e-10 -num_descriptions 10 -num_threads 12  -out $SPECIES.outinner85.clean_blastx.out.txt

  perl  $DIR2/outinner_blastx_parse.pl --blastx  $SPECIES.outinner85.clean_blastx.out.txt  --outinner  $SPECIES.outinner85

# 2.1.5  Building examplars (choosing representative repeats)

perl $DIR2/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile $SPECIES.result85 \
--innerfile passed_outinner_sequence.fasta  --sequencefile $SPECIES.fasta

makeblastdb -in lLTRs_Seq_For_BLAST.fasta  -dbtype nucl

blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000  \
-out lLTRs_Seq_For_BLAST.fasta.out

makeblastdb -in  Inner_Seq_For_BLAST.fasta  -dbtype nucl

blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta  -evalue 1e-10 -num_descriptions 1000 \
-out  Inner_Seq_For_BLAST.fasta.out

perl $DIR2/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out  \
--step3 CRL_Step3_Passed_Elements.fasta --final LTR85.lib --pcoverage 90 --pidentity 80

cd .. ; ln -s 99/LTR99.lib ; ln -s 85/LTR85.lib

#########################################################################################################################

  # Mask LTR99 out of LTR85, so that it is now repeated when concatenating

  /gpfs/project/projects/qggp/src/repeatmasker-4-0-9/RepeatMasker -pa 12 -lib LTR99.lib -dir . LTR85.lib ; wait

  perl $DIR2/remove_masked_sequence.pl  --masked_elements  LTR85.lib.masked  --outfile  FinalLTR85.lib ; wait

  cat LTR99.lib FinalLTR85.lib > allLTR.lib


######################

# 3. Collecting repetitive sequences by RepeatModeler 

  cd .. ; ln -s ../$SPECIES.fasta 

  cat ltr/allLTR.lib mites/MITE.lib > allMITE_LTR.lib

  /gpfs/project/projects/qggp/src/repeatmasker-4-0-9/RepeatMasker -pa 12 -lib allMITE_LTR.lib -dir . $SPECIES.fasta ; wait


  perl $DIR2/rmaskedpart.pl  $SPECIES.fasta.masked  50  >  umseqfile ; wait

  /gpfs/project/projects/qggp/src/repeatmodeler-1.0.11/BuildDatabase -name umseqfiledb -engine ncbi ./umseqfile ; wait

  nohup /gpfs/project/projects/qggp/src/repeatmodeler-1.0.11/RepeatModeler -pa 12 -database umseqfiledb >& umseqfile.out  ; wait

  mv */consensi.fa* . 

  perl /gpfs/project/projects/qggp/src/maker/CRL_Scripts1.0/repeatmodeler_parse.pl --fastafile consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta  \
  --identities repeatmodeler_identities.fasta  ; wait

	# Search if repeats are known
  makeblastdb -in ~/C34_PS/data/databases/transposons/tpases020812.fa  -dbtype prot 

  blastx -query repeatmodeler_unknowns.fasta -db ~/C34_PS/data/databases/transposons/tpases020812.fa  -evalue 1e-10 -num_descriptions 10 \
  -num_threads 12  -out modelerunknown_blast_results.txt   ; wait

  perl $DIR2/transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknowns.fasta ; wait

  mv  unknown_elements.txt  ModelerUnknown.lib

  cat  identified_elements.txt  repeatmodeler_identities.fasta  > ModelerID.lib


########

# 4. Exclusion of gene fragments 


makeblastdb -in ~/C34_PS/data/databases/protein/alluniRefprexp070416  -dbtype prot

	# some .libs are one-liners and others are multiliners. Causes a bug in last step! Correct
cat  allMITE_LTR.lib ModelerID.lib  > KnownRepeats.lib ; wait

cat KnownRepeats.lib  ModelerUnknown.lib > allRepeats.lib ; wait 

blastx -query allRepeats.lib -db /home/guerrer/C34_PS/data/databases/protein/alluniRefprexp070416  \
 -evalue 1e-10 -num_descriptions 10  -num_threads 12 -out allRepeats_blast_results.txt ; wait

 
 #$DIR2/ProtExcluder1.2/ProtExcluder.pl  -option  ModelerUnknown.lib_blast_results.txt  ModelerUnknown.lib # Doesn't work. Replaced by following:
/gpfs/project/projects/qggp/src/maker/CRL_Scripts1.0/ProtExcluder1.2/run_ProtExcluder_ricardo.sh allRepeats_blast_results.txt allRepeats.lib | sh ; wait


## This is just a curiosity, to calculate how much is masked or not:
  /gpfs/project/projects/qggp/src/repeatmasker-4-0-9/RepeatMasker -x -pa 12 -lib allRepeats.lib.noProtFinal -dir $SPECIES.masked_final  $SPECIES.fasta 


cd .. ; ln -s repeat-masking/allRepeats.lib.noProtFinal
