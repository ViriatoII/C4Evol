

ortho_run="orthofinder/pre_blasts/OrthoFinder/Results_pre_blasts3_tetra_raxml/WorkingDirectory/"    #7_tetra_raxml





# Colect gene trees of selected HOGs
cat   $ortho_run/Trees_ids/OG00* > all_ogs.tre



cat  all_ogs.tre | tr ":" "\n"  | head -200  |grep "_" | sed 's/.*(\|.*,//' | sed 's/\(.*\)_\(.*\)_/\1  \1_\2_/'  | sort   > gene_sp_translation


sed 's/_.._[0-9][0-9]*-RA//g ; s/_HI3_[0-9][0-9]*-RA//g'  all_ogs.tre > all_ogs.trees 


java -D"java.library.path=src/ASTRAL/A-pro/ASTRAL-MP/lib" -jar src/ASTRAL/A-pro/ASTRAL-MP/astral.1.1.5.jar \
-a ../gene_name_translations2.tsv -i all_ogs.trees -t1 -o all_ogs.final.tree  2>out.log
