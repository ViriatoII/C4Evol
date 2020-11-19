#PBS -l select=1:ncpus=4:mem=6G
#PBS -l walltime=48:00:00
#PBS -A "PROJECT"

# This script fills empty positions in snp_matrix.csv created by compare_vcfs.py, since those empty positions might mean either of two things:
# a) There is no variant. Genotype has same sequence as reference in that position.
# b) No read maps to that position. There is probably a large SV in that area, or maybe sequencing is incomplete.


module load SamTools bcftools/1.10.2  #bedtools/2.26.0


cd $PBS_O_WORKDIR         # ENTER SCRIPT DIRECTORY


## Variable INPUTs 
contig="LR618883.1" #v4: NC_024468.2"   # Only working for one specific contig/scaffold. Automate? 
ref="extract_regions/refs/B73v5.fasta"
align_dir="extract_regions/pseudomols/" # Should have a folder with all mapped bams (BREED_vs_REF.bam)
refname=`echo "$ref" | sed 's/.*\/// ; s/.fasta//'`

breeds_ite=`head -1  snp_matrix.csv | sed 's/P.*REF,//' | tr "," "\n" `


rm -r  tmp snp_matrix_corrections.csv  ; mkdir tmp 


# Search only positions where "-" exists
for line in `grep "-" snp_matrix.csv `  # ERASE HEAD 
	do
	pos=`echo "$line" | cut -f1 -d "," `
	len=`echo "$line" | cut -f2 -d "," |  wc | sed 's/.* //'  `	

	# prepare a bed file to extract exact nucleotide in question
        echo "$contig $((pos-1)) $((pos+len-2))" | tr " " "\t" > tmp/tmp.bed ; wait

	# Position and ref for later pasting
	#	echo "$pos" >> tmp/positions.txt
	echo $line  | cut -f1,2 -d,  >> tmp/positions.txt  

	colnum=2 # counter starts at 2 (after pos and ref)

	# Get input bams and for each breed
        for breed in $breeds_ite   
                do
		
		# counter updates per breed, resets per line
		let "colnum+=1"

		# Only appy to empty SNPs.
		if [ `echo "$line" | cut -f$colnum -d, ` == "-" ] ; then
	
			#Extract mappings spanning position
                	samtools view  -h $align_dir/ragoo_${breed}_trimreads_vs_$refname/contigs_against_ref.bam $contig:$pos-$pos  > tmp/reads.sam  ; wait

	                #Generate consensus from mapped reads
        	        samtools mpileup  -uf  $ref  tmp/reads.sam | bcftools call -c  --threads 4 | vcfutils.pl vcf2fq > tmp/consensus.fastq ; wait
                	seqtk seq -aQ64 tmp/consensus.fastq > tmp/consensus.fasta  # transform to fasta. To make low qual reads(<20) --> N    -q20 -n N
		
	                #Extract sequence in exact coordinates
        	        sequence2=`bedtools getfasta -fi  tmp/consensus.fasta  -bed tmp/tmp.bed | tail -1 `

		else    #Use called SNP of snp_matrix line.
			sequence2=`echo "$line" | cut -f$colnum -d, ` ; fi 

                echo "$sequence2" >> tmp/$breed.csv ; wait
       		rm tmp/consensus*

                done
	done 

# Concat header with paste of all files (POS,REF) + BREEDS  (letters capitalized)
head  -1 snp_matrix.csv > snp_matrix_corrections.csv 
paste -d "," tmp/positions.txt  tmp/*csv  | sed 's/c/C/g ; s/a/A/g ; s/g/G/g ; s/t/T/g' > t.csv  
#rm -r tmp 

# Add ignored positions to corrections 
sed 's/,.*//' tmp/positions.txt > t2
grep -vf t2  snp_matrix.csv | grep -v "REF" >>  t.csv 

sort -t "," -Vk1 t.csv >>  snp_matrix_corrections.csv 

