#!/usr/bin/env python

import pandas as pd 
import sys 


if len(sys.argv) < 2 : 
    print( '\n  ./extract_orthogroup_genes.py  N0.tsv fasta_dir HOG1 HOG2 HOG3 ... \n\n - fasta_dir: Directory where your species protein/cds files are located ')
    exit()

df  = pd.read_csv(sys.argv[1], sep="\t")

fastadir = sys.argv[2] # "/gpfs/project/projects/qggp/C34_PS/experiments/annotation/orthofinder/abinitio_ann/"
selected = pd.read_csv(sys.argv[3:])
#selected = ["N0.HOG0003237","N0.HOG0025962","N0.HOG0028581","N0.HOG0034511","N0.HOG0034527","N0.HOG0021276","N0.HOG0026267","N0.HOG0030986","N0.HOG0025956","N0.HOG0021770","N0.HOG0022963","N0.HOG0024157","N0.HOG0024382"]

df = df.loc[df.HOG in selected]


# Go through relevant groups of Orthofinder N0
for row in range(0, len(df)):
    
    # Orthogroup in question
    HOG = df.iloc[row,0]

    # Detect species where HOG exists   
    cols = df.iloc[row,3:].notna()
    idxs =[ i for i, val in enumerate(cols) if val == True ]
    species = list(cols[idxs].index)
   
 
    for sp in species:
        # Get all genes of that species in Ogroup
        genes= df.loc[row,sp].split(",")
        
        #### Extract gene from species fasta #############
        infasta = f"{fastadir}/{sp}.fasta"
        
        for seq_record in SeqIO.parse(infasta, "fasta"):
            if seq_record.id in  genes:
                
                # Extract gene
                SeqIO.write(seq_record, HOG + ".fasta", "fasta")
                   

