module load /software/modules/tools/Java/1.8.0 /software/modules/tools/Perl/5.18.1  /software/modules/compiler/gcc/8.1.0


canu -d $species -p $species  minReadLength=2000  rawErrorRate=0.15 correctedErrorRate=0.1   \
        -genomeSize=450m -nanopore-raw $inpath/*.fastq \
        gridOptions="-A C4Evol -l walltime=48:00:00 -r y -V" \
        preExec="module load Java/1.8.0 ; module load Perl/5.18.1 ; module load gcc/8.1.0" \
        gridEngineResourceOption="-l select=1:ncpus=THREADS:mem=MEMORY" #30G"



# Map linked reads to assembly for polishing, purging and further scaffolding 

longranger mkref  $pseudohap_assembly; wait
        
mv refdata-$pseudohap_assembly references/refdata-${species}

longranger align --id=${species} \
        --reference="references/$IN_TYPE/refdata-${species}"  \
        --fastqs="$FASTQs"        \
        --localcores=8 --localmem=20   


# Polishing with PILON ###########


java -jar -Xmx300G /software/pilon/1.22/pilon-1.22.jar \
 --genome $ToPilon     --bam $mapped10x \
 --changes --threads 8 --output ${species}_final --outdir $species










