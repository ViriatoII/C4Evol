# FIRST BLASTS
src/OrthoFinder/orthofinder.py -t 64  -f pre_blasts -S blast  -og  -n pre_blasts1

# ADD SPECIES (BLAST)
src/OrthoFinder-2.5.1/orthofinder.py -t 64  -f pre_blasts/ -S blast -og  -b pre_blasts/OrthoFinder/Results_pre_blasts1  -n pre_blasts2


# Add rest of species (blast) and prepare msa for next run 
src/OrthoFinder-2.5.1/orthofinder.py -t 64  -f pre_blasts/ -S blast -b pre_blasts/OrthoFinder/Results_pre_blasts1 -n pre_blasts2_tetra_msa -M msa 


# Apply RAXML for gene and species trees
src/OrthoFinder-2.5.1/orthofinder.py -t 64  -M msa -S blast  -fg pre_blasts/OrthoFinder/Results_pre_blasts2_tetra_msa  -n pre_blasts3_tetra_raxml   -T raxml #-ng