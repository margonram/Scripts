#!/usr/bin/env python                                                       #
#                                                                           #
# Script BinGenome.py                                                       #
#                                                                           #
#                                                                           #
# Input : window width, name, genome                                        # 
# Output: binning                                                           #
#                                                                           #
#                                                                           #
# by Mar Gonzalez @ CRG (2018)                                              #
#############################################################################

import sys
import os

# Bin width
W = 2000

# Name
name = "mouse"

# Open bed files
genome = "/usr/local/molbio/indexes/mm10/ChromInfo.txt"
f =open(genome)

def get_fields(f):
    chr_name = []
    chr_length = []
    for line in f:
        name,length = line.split()
        chr_name.append(name)
        chr_length.append(length)
    chr_length = map(int, chr_length)
    f.close()
    return chr_name,chr_length

def binarize(chr_name,chr_length,out):
    # Iterate each line (segment)
    for k in range(len(chr_length)):
        i = 1
        j = W
        # Binning of the segment
        while i < chr_length[k] and j < chr_length[k]:
            out.write("%s\t%d\t%d\t%s\n"%(chr_name[k],i,j,name))
            i += W
            j += W
        if i < chr_length[k] and j > chr_length[k]:
            out.write("%s\t%d\t%d\t%s\n"%(chr_name[k],i,chr_length[k],name))

# Get chr, coordinates and state of the segments
chr_name,chr_length = get_fields(f)

# Save output
out = open(name+'_'+`W`+'_bin.bed','w')

# Bin the segments
binarize(chr_name,chr_length,out)
