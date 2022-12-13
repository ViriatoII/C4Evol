

in_draft_assembly="phasing/purge_haplotigs/$species/curated.fasta"
in_barcoded_reads="longranger/basic/$species/outs/barcoded.fastq.gz"


        # create additional fastq specific for later step (has barcodes on the fastq header)
        gunzip -c  $in_barcoded_reads | \
                perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print "$_\n" if($flag);}' > interleaved.fastq &

        samtools faidx $in_draft_assembly
        bwa index $in_draft_assembly

        bwa mem -t $cpus -p -C $in_draft_assembly $in_barcoded_reads > curated_genome-inter.sam ; wait

        samtools view -Sb curated_genome-inter.sam | samtools sort -@ $cpus2 -tbX - -o mapped_reads_curated_genome-renamed.sortbx.bam    ; wait

        tigmint-molecule mapped_reads_curated_genome-renamed.sortbx.bam | sort -k1,1 -k2,2n -k3,3n > curated_genome-renamed.molecule.bed ; wait

        tigmint-cut -p$cpus -o corrected_assembly.tigmint.fa $in_draft_assembly curated_genome-renamed.molecule.bed; wait

        cat corrected_assembly.tigmint.fa | perl -ne 'chomp;if(/^>/){$ct++;print ">$ct\n";} else {print "$_\n";}' >corrected_assembly.tigmint-renamed.fa ; wait

        samtools faidx corrected_assembly.tigmint-renamed.fa &

        bwa index corrected_assembly.tigmint-renamed.fa; wait

        bwa mem -t $cpus corrected_assembly.tigmint-renamed.fa -p interleaved.fastq  > barcoded_mapped_tigmint.sam ; wait

        samtools view -Sb barcoded_mapped_tigmint.sam  |  samtools sort -n -@ $cpus2 - -o barcoded_mapped_sortedbyname_tigmint.bam ; wait

        rm barcoded_mapped_tigmint.sam

	# put Absolute path of bam file to be input of arcs next:
        readlink -f barcoded_mapped_sortedbyname_tigmint.bam   > list_bam_file.txt   

        src/arcs/Arcs/arcs -f corrected_assembly.tigmint-renamed.fa -a list_bam_file.txt -c 3 -m 10-10000 ; wait 

        python ../makeTSVfile.py corrected_assembly.tigmint-renamed.fa.scaff_s98_c3_l0_d0_e30000_r0.05_original.gv reads.tigpair_checkpoint.tsv corrected_assembly.tigmint-renamed.fa ; wait

        touch empty.fof

        src/LINKS/bin/LINKS -f corrected_assembly.tigmint-renamed.fa -s empty.fof -k 20 -b reads -l 5 -t 2 -a 0.3