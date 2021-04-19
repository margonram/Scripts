#!/usr/bin/env python                                                       #
#                                                                           #
# Script Merge2bams.py                                                      #
#                                                                           #
#                                                                           #
# Input : bam1,bam2,name                                                    # 
# Output: merge bam                                                         #
#                                                                           #
#                                                                           #
# by Mar Gonzalez @ CRG (2019)                                              #
#############################################################################

import sys
import os
from termcolor import cprint

# Control arguments
count = len(sys.argv[1:])
if count < 4:
    cprint("4 parameters needed", 'magenta', attrs=['bold'], file=sys.stderr)
    cprint("Merge2bams.py <genome assembly> <bam1> <bam2> <name>", 'magenta', attrs=['bold'], file=sys.stderr)
    exit()

# Genome
genome = sys.argv[1]

# Open bam files
bam1 = sys.argv[2]

bam2 = sys.argv[3]

# output name
name = sys.argv[4]

# control genome
if genome == "mm9":
    maps = "02_map_files/mm9"
    chrominfo = "$MOUSE9/ChromInfo.txt"

elif genome == "mm10":
    maps  = "02_map_files/mm10"
    chrominfo = "$MOUSE10/ChromInfo.txt"

elif genome == "hg19":
    maps  = "02_map_files/hg19"
    chrominfo = "$HUMAN9/ChromInfo.txt"

elif genome == "hg38":
    maps  = "02_map_files/hg38"
    chrominfo = "$HUMAN38/ChromInfo.txt"

else:
    cprint ("incorrect genome", 'grey', attrs=['bold'], file=sys.stderr)
    exit()

# check if files exist
if os.path.isfile(bam1) == 0:
    cprint ("bam1 does not exist", 'grey', attrs=['bold'], file=sys.stderr)
    exit()
if os.path.isfile(bam2) == 0:
    cprint ("bam2 does not exist", 'grey', attrs=['bold'], file=sys.stderr)
    exit()

# Merge bam files
command = 'samtools view -H '+bam1+' > header.sam'
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'samtools cat -h header.sam -o '+maps+'/'+name+'.bam '+bam1+' '+bam2 
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'rm -f header.sam'
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

# Generate profile
command = 'buildChIPprofile -v '+chrominfo+' '+maps+'/'+name+'.bam '+name
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)
