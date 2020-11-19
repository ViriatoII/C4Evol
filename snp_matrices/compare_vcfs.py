#!/usr/bin/env python3

import vcf     # pip3 install --user  PyVCF
import sys


if len(sys.argv) < 2: 
	print("""compare_vcfs.py  vcf1 vcf2 vcf3
	Expects the breed names to be prefixes (e.g:Breed1_ref.vcf)
	Vcfs must come from alignment to same reference. """)
	exit()

# Gets files and extracts breed name prefix
files = sys.argv[1:]
breeds = [fil.split("_")[0] for fil in files] 



def vcfs_to_dics(breeds=breeds):
	"""Reads input vcf files and saves REF,ALT nucleotides of every position into a sub-dictionary. 
	Each sub-dic is stored in output (all_dics), and can accessed through all_dics[breed] """
	
	# all dictionaries stored here, where keys are breed names
	all_dics =  {}

	for filenum, breed in enumerate(breeds):
		# Each breed is a sub-dictionary
		all_dics[breed]={}
		
		# Read input vcfs 
		vcf_reader = vcf.Reader(open(files[filenum], 'r'))

		# key of dictionaries is position, value is list of [REF, ALT] alleles.
		for record in vcf_reader:
			all_dics[breed][record.POS] = [record.REF, record.ALT ]
	
	return(all_dics)


def dics_to_matrix(all_dics, breeds, outfile="snp_matrix.csv"):
	""" Reads all_dics and creates an output snp_matrix """

	# Append all positions of each subdict (.values() to get subdicts )
	posis = [pos for dic in all_dics.values() for pos in dic.keys()]

	# Unique the positions
	tot_pos = set(posis)

	with open(outfile, "w") as file:
	        head1 = "".join(["," + breed for breed in breeds])

        	# a header
	        file.write(f"POS,REF{head1}\n")

        	#iterate through all SNPs positions
	        for i in sorted(tot_pos): 

        	        nucl = {}

	                for breed in breeds:
                        	# Get reference from any dict
                	        try:  ref = all_dics[breed][i][0]
        	                except: pass

	                        # Get SNPs or assign to "-"  if they don't exist
                        	try:    nucl[breed] = all_dics[breed][i][1][0]
                	        except: nucl[breed] = "-"     # Is it equal to ref or not existant?? -> 00_fill_snp_positions.sh

        	        alleles = "".join(["," + str(nucl[breed]) for breed in breeds])

	                file.write(f"{i},{ref}{alleles}\n")


# Read vcfs 
all_dics = vcfs_to_dics(breeds)

#Create output : snp_matrix 
dics_to_matrix(outfile="snp_matrix.csv", all_dics = all_dics, breeds = breeds )

