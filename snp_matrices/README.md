# C4Evol
Comparative Genomics in Brassicaceae Family to understand C3-C4 intermediate evolution

## Comparing VCF files of related breeds.

This pipeline uses bam files of genotypes aligned to a single reference. Those alignments are used to call SNPs in vcf files.     

The various vcf files are used to build a SNP matrix with all called positions and variants of each genotype. However, uncalled variants can have two reasons:    

a) There is no variant. Genotype has same sequence as reference in that position.    
b) No read maps to that position. There is probably a large SV in that area, or maybe sequencing is incomplete. 

The script fill_ambiguous_positions.sh goes back to the bam files and extracts the sequence from the reads, which should equal REF. If it doesn't, the position is left empty, with final sureness (case b).

