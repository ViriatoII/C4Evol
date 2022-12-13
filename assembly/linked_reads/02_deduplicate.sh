
refpath="$species.fasta"


ln -s longranger/align/$species/outs/possorted_bam.bam  aligned.bam


# Run purge_haplotigs
purge_haplotigs  hist -t 24  -b aligned.bam  -g $refpath
purge_haplotigs  cov  -i aligned.bam.gencov -l 1 -m 22 -h 185
purge_haplotigs  purge  -g $refpath  -c coverage_stats.csv -t 24 -a 70