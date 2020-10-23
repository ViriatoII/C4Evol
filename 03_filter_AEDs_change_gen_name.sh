#PBS -l select=1:ncpus=2:mem=2G
#PBS -l walltime=6:00:00
#PBS -A "C4Evol"

module load  gcc/8.1.0 Python/3.6.5  #intel/xe2019

maker_scripts="/gpfs/project/projects/qggp/src/maker2/bin/"

# Get species from directory name 
species=`readlink -f . | sed 's/.*\/maker\/// ; s/\/.*//'`

for sp in $species   # change this loop to run from master directory 
	do	

	################ rename maker genes ##################################################
        
        # Get initials from species name 
        if [  $sp == b_nigra ]          ; then PREF="BI"  ; else
        if [  $sp == d_tenuisiliqua ]   ; then PREF="DS"    ; else
        if [  $sp == c_gynandra ]       ; then PREF="BGcg"  ; else
        if [  $sp == d_acris ]          ; then PREF="CAda"  ; else
        if [  $sp == m_spinosa ]        ; then PREF="MP"  ; else
        PREF=`echo "${sp:0:1}${sp:2:1}"  |  tr "[a-z]" "[A-Z]"` ; fi ; fi ; fi ; fi ;fi

	# Create simple gene names
	$maker_scripts/maker_map_ids --prefix ${PREF}_ --justify 8 *.all.gff > gene_names.map
	# Rename maker output files with simple gene names
        $maker_scripts/map_gff_ids   gene_names.map    *.all.gff
        $maker_scripts/map_fasta_ids gene_names.map    *.all.maker.proteins.fasta
        $maker_scripts/map_fasta_ids gene_names.map    *.all.maker.transcripts.fasta


        ############ Filter small and bad proteins out #############################################

        # Get entries based on AED with maker scripts				 	
        /gpfs/project/projects/qggp/src/maker2/bin/quality_filter.pl  -a 0.5 *.all.gff  > filtered.gff ; wait 

        # Get gene ids     			# Note: it needs "-RA" suffix 
        grep $'\tgene\t' filtered.gff  | sed 's/.*ID=// ; s/;Name.*/-RA/' > tmp ; wait

        ## Filter with my custom script  # Filters out proteins with more than 0.5 AED and less than 50 nucltds
        filter_AEDs.py tmp *.all.maker.proteins.fasta ; wait 

        rm tmp filtered.gff

	cd ../..

done
