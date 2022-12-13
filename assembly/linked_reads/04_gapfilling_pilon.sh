
module load gcc/8.1.0 Miniconda/3

conda activate gapfilling


assem="/gpfs/project/projects/qggp/C34_PS/experiments/scaffolding/links/$species/d20000/${species}_scaf.scaffolds.fa"
inreads="/gpfs/project/projects/qggp/C34_PS/experiments/scaffolding/links/$species/reads.fa"


tgsgapcloser --ne  --scaff $assem --reads $inreads  --output ${species}_scaff




longranger mkref  $pseudohap_assembly; wait
        
mv refdata-$pseudohap_assembly references/hybrid/refdata-${species}

longranger align --id=${species} \
        --reference="references/hybrid/refdata-${species}"  \
        --fastqs="$raw_linkedreads.fastq"        \
        --localcores=8 --localmem=20   


# Polishing with PILON ###########

java -jar -Xmx300G /software/pilon/1.22/pilon-1.22.jar \
 --genome $ToPilon     --bam $mapped10x \
 --changes --threads 8 --output ${species}_final --outdir $species










