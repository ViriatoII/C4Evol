#!/usr/bin/env python3

## calculate N50 from fasta file     ### 
## N50 = contig length so that half of the contigs are longer and 1/2 of contigs are shorter

#import commands
import sys
import os
from itertools import groupby
import numpy
import zipfile
import gzip
# from Bio import SeqIO

lengths = []
n_number = 0

the_input= sys.argv[1]

if str(the_input).endswith("gz"):   # Open gunziped files 

	with gzip.open(the_input) as fasta:
		## parse each sequence by header: groupby(data, key)
		faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

		for record in faiter:
			## join sequence lines
			seq = "".join(s.strip() for s in next(faiter))
			lengths.append(len(seq))
			n_number = n_number + seq.count("N")



else:
	with open(the_input) as fasta:
    	## parse each sequence by header: groupby(data, key)
		faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))

		for record in faiter:
        		## join sequence lines
			seq = "".join(s.strip() for s in next(faiter))
			lengths.append(len(seq))
			n_number = n_number + seq.count("N")




Total_size = int(sum(lengths))

print(f"Contigs: {len(lengths)}")

## sort contigs longest>shortest
all_len = numpy.array( sorted(lengths, reverse=True) )  # sorted lengths ; numpy.array is faster
csum = numpy.cumsum(all_len)


## Count contigs > threshhold
print(f"Contigs >= 5000 bps: {sum(all_len >=   5000)}")
print(f"Contigs >= 10000 bps: {sum(all_len >=  10000)}")
print(f"Contigs >= 25000 bps: {sum( all_len >= 25000)}")
print(f"Contigs >= 50000 bps: {sum( all_len >= 50000)}")
print(f"Contigs >= 100000 bps: {sum( all_len >= 100000)}")

print(f"Largest contig: {all_len[0]}")
print(f"Total size: {Total_size}")

print(f"Length >= 5000 bps: {sum(all_len[all_len >=    5000])}")     # Boolean indexing
print(f"Length >= 10000 bps: {sum(all_len[all_len >=  10000])}")
print(f"Length >= 25000 bps:{ sum(all_len[all_len >=  25000])}")
print(f"Length >= 50000 bps: {sum(all_len[all_len >=  50000])}")
print(f"Length >= 100000 bps: {sum(all_len[all_len >= 100000])}")

## N50
# Half the genome length
n2=int(sum(lengths)/2)

# get index for cumsum >= N/2
csumn2 = min(csum[csum >= n2])
ind    = numpy.where(csum == csumn2)
n50    = all_len[ int(ind[0]) ]

print(f"N50:{n50}")


## N90
nx90=int(sum(lengths)*0.90)

## index for csumsum >= 0.9*N
csumn90 = min(csum[csum >= nx90])
ind90   = numpy.where(csum == csumn90)
n90     = all_len[int(ind90[0]) ]

print(f"N90: {n90}")

L50 = int(ind[0] + 1  )
print(f"L50: {L50}")

L90 = int(ind90[0] +1  )
print(f"L90: {L90}")

print(f"Ns per 100 kbs: {(n_number*100000)/Total_size}")

## write lengths to file (for histogram)
with open('all_lengths.txt', 'w') as handle:
	handle.write('\n'.join(str(i) for i in lengths))
