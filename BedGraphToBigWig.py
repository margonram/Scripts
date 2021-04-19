#!/usr/bin/env python                                                       #
#                                                                           #
# ScriptBedGraphToBigWig.py                                                 #
#                                                                           #
#                                                                           #
# Input : genome bedgraph sortedBedgraph name                               # 
# Output: merge bam                                                         #
#                                                                           #
#                                                                           #
# by Mar Gonzalez @ CRG (2019)                                              #
#############################################################################

import sys
import os
import re
from termcolor import cprint

# Control arguments
count = len(sys.argv[1:])
if count < 4:
    cprint("4 parameters needed", 'magenta', attrs=['bold'], file=sys.stderr)
    cprint("BedGraphToBigWig.py <genome> <BedGraph.gz> <sortedBedgraph.gz> <name>", 'magenta', attrs=['bold'], file=sys.stderr)
    exit()

# Genome
genome = sys.argv[1]

# bedgraphs
bedgraph_gz = sys.argv[2]
sortedbedgraph_gz = sys.argv[3]

# Open bam files
name = sys.argv[4]

# control genome
if genome == "mm9":
    chrominfo = "$MOUSE9/ChromInfo.txt"

elif genome == "mm10":
    chrominfo = "$MOUSE10/ChromInfo.txt"

elif genome == "hg19":
    chrominfo = "$HUMAN9/ChromInfo.txt"

elif genome == "hg38":
    chrominfo = "$HUMAN38/ChromInfo.txt"

else:
    cprint ("incorrect genome", 'grey', attrs=['bold'], file=sys.stderr)
    exit()


# check if files exist
if os.path.isfile(bedgraph_gz) == 0:
    cprint ("BedGraph does not exist", 'grey', attrs=['bold'], file=sys.stderr)
    exit()

# Decompressed names
bedgraph=re.sub(r'.gz','',bedgraph_gz )
sortedbedgraph=re.sub(r'.gz','',sortedbedgraph_gz )

# Obtain Bigwig
command = 'mkdir Bigwigs'
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'gzip -d '+bedgraph_gz
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'grep -v track '+bedgraph+' | sort -k1,1 -k2,2n > '+sortedbedgraph
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'bedGraphToBigWig '+sortedbedgraph+' '+chrominfo+' Bigwigs/'+name+'.bw'
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'gzip '+bedgraph
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)

command = 'gzip '+sortedbedgraph
cprint (command, 'magenta', attrs=['bold'], file=sys.stderr)
os.system(command)
