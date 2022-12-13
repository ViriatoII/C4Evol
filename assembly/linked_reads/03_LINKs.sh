
assembly="phasing/purge_haplotigs/$sp/curated.fasta"     # purged assembly
reads="longreads/$sp/filtered_reads.fasta"


####### Prepare longreads for scaffolding #########

# if ONT data for scaffolding:
cat  $nanopore | NanoFilt -q 10 -l 1000 --headcrop 50 > longreads/$sp/filtered_reads.fastq


# else if PacBio data for scaffolding:
canu -d $sp -p $sp  minReadLength=2000  rawErrorRate=0.15 correctedErrorRate=0.1   \
        -genomeSize=450m -nanopore-raw $inpath/*.fastq \
        gridOptions="-A C4Evol -l walltime=48:00:00 -r y -V" \
        preExec="module load Java/1.8.0 ; module load Perl/5.18.1 ; module load gcc/8.1.0" \
        gridEngineResourceOption="-l select=1:ncpus=THREADS:mem=MEMORY"

ln -s $sp/$sp.correctedReads.fasta longreads/$sp/filtered_reads.fasta


#### LINKs scaffolding from long reads #########


#####Scaffold ##################

# Prepare reads
cat $reads | perl -ne 'chomp;if(/^>/){$ct++;print ">$ct\n";} else {print "$_\n";}' > reads.fa

# Rename assembly contigs
perl -ne 'chomp;if(/^>/){$ct++;print ">$ct\n";} else {print "$_\n";}' $assembly > assembly-renamed.fa

# Input for LINKs
readlink -f reads.fa   > reads.fof   # put ABS path of bam file to be input of arcs next:

# Bloom filter separately? Not necessary
/gpfs/project/projects/qggp/src/LINKS/releases/links_v1.8.7/tools/writeBloom.pl -f  assembly-renamed.fa -k 18

# Iteratively scaffold with different distances # Play with -d and -t parameters
for dist in 1000 2500 5000 7500 10000 15000 20000  
        do
        if [ $dist == 1000 ];  then  paramt=36 ; ln -s assembly-renamed.fa scaffolded.fa ; fi
        if [ $dist == 2500 ];  then  paramt=30 ; fi
        if [ $dist == 5000 ];  then  paramt=18 ; fi
        if [ $dist == 7500 ];  then  paramt=14 ; fi
        if [ $dist == 10000 ]; then  paramt=10 ; fi
        if [ $dist == 15000 ]; then  paramt=5 ; fi
        if [ $dist == 20000 ]; then  paramt=2 ; fi

        mkdir d${dist} ; cd  d${dist}

        /gpfs/project/projects/qggp/src/LINKS/bin/LINKS -f ../scaffolded.fa -d $dist \
                -s ../reads.fof -k 18 -b ${sp}_scaf -l 5 -t $paramt -a 0.3  #-r $sp.scaf.bloom

        cd - ; unlink scaffolded.fa ; ln -s d${dist}/${sp}_scaf.scaffolds.fa scaffolded.fa
