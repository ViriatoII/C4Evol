
        BBMAP="/gpfs/project/projects/qggp/C34_PS/src/bbmap/bbmap.sh"
        IN_10x="/gpfs/project/projects/qggp/C34_PS/data/sequencing/novogene/$Species"

        myin1=`ls $IN_10x/*R1_001.fastq.gz`
        myin2=`ls $IN_10x/*R2_001.fastq.gz`

        out1=`echo "$myin1" | sed 's/.*\///' `   
        out2=`echo "$myin2" | sed 's/.*\///' `

        mkdir $IN_DIR10x/cleaned $IN_DIR10x/chloroplast $IN_DIR10x/mitochondria

        ##filtering chloroplast##
        PLASTID="data/literature/b_nigra/chloroplast/b_nigra_chloro.fasta"

        $BBMAP  in1=$myin1 \
                in2=$myin2  \
                ref=$PLASTID \
                outu1=$IN_10x/unmapped_chloro_R1.fastq.gz \
                outu2=$IN_10x/unmapped_chloro_R2.fastq.gz \
                outm1=$IN_10x/chloroplast/mapped_chloro_R1.fastq.gz \
                outm2=$IN_10x/chloroplast/mapped_chloro_R2.fastq.gz

	# Filter Mitochondreon
        PLASTID=/gpfs/project/projects/qggp/C34_PS/data/literature/e_sativa/mitochondreon/e_sativa_mito.fasta

        
        $BBMAP in1=$IN_10x/unmapped_chloro_R1.fastq.gz \
                in2=$IN_10x/unmapped_chloro_R2.fastq.gz \
                ref=$PLASTID \
                outu1=$IN_10x/cleaned/$out1     \
                outu2=$IN_10x/cleaned/$out2 \
                outm1=$IN_10x/mitochondria/mapped_mito_R1.fastq.gz      \
                outm2=$IN_10x/mitochondria/mapped_mito_R2.fastq.gz