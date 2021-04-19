#!/usr/bin/env python                                                       #
#                                                                           #
# Script MergePeaks.py                                                      #
#                                                                           #
#                                                                           #
# Input : two bed files                                                     # 
# Output: merged peaks,merged peaks 1, merged peaks 2                       #
#                                                                           #
#                                                                           #
# by Mar Gonzalez @ CRG (2020)                                              #
#############################################################################


import sys
import os
import re


# Open bed files
peaks1 = sys.argv[1]
f1 =open(peaks1)

peaks2 = sys.argv[2]
f2 =open(peaks2)

name1 = sys.argv[3]
name2 = sys.argv[4]

# Define function to get fields from 4 column files

def get3_fields(f):
    chrn = []
    start = []
    end = []
    next(f)
    for line in f:
        chrni,starti,endi = line.split('\t')
        chrn.append(chrni)
        start.append(starti)
        end.append(endi)
    f.close()
    #for i in range(len(chrn)):
    #    chrn[i]=re.sub(r'chr','',chrn[i])
    #chrn = map(int, chrn)
    start = map(int, start)
    end = map(int, end)
    return chrn,start,end

# Define function to unify coordinates of peaks

def unify_coordinates(chrn,start,end):
    chrnu = []
    startu = []
    endu = []
    k = 1
    i = 0
    while i < len(chrn):
        endsave = []
        try:
            chrn[i]
        except IndexError:
            break
        try:
            chrn[i+k]
        except IndexError:
            break
        endsave.append(end[i])
        while chrn[i] == chrn[i+k] and end[i] >= start[i+k]:
            endsave.append(end[i+k])
            end[i] = max(endsave)
            k += 1
            try:
                chrn[i+k]
            except IndexError:
                break
        chrnu.append(chrn[i])
        startu.append(start[i])
        endu.append(max(endsave))
        i += k
        k = 1
    return chrnu,startu,endu

    

# Define function to overlap peaks

def overlap_peaks(chrnu,startu,endu,chrnS,startS,endS):
    chrnuS=[]
    startuS=[]
    enduS=[]
    for i in range(len(chrnu)):
        for j in range(len(chrnS)):
            if chrnu[i] == chrnS[j] and startu[i] <= startS[j] and endu[i] >= endS[j]:
                chrnuS.append(chrnu[i])
                startuS.append(startu[i])
                enduS.append(endu[i])
    return chrnuS,startuS,enduS


# Get fields from peak files

chrn1,start1,end1 = get3_fields(f1)
chrn2,start2,end2 = get3_fields(f2)

chrn = chrn1 + chrn2
start = start1 + start2
end = end1 + end2

# Order peaks

indexS = sorted(range(len(start)), key=lambda k: start[k])

chrnsort = []
startsort = []
endsort = []
for i in indexS:
    chrnsort.append(chrn[i])
    startsort.append(start[i])
    endsort.append(end[i])

indexC = sorted(range(len(chrnsort)), key=lambda k: chrnsort[k])

chrnsortsort = []
startsortsort = []
endsortsort = []
for i in indexC:
    chrnsortsort.append(chrnsort[i])
    startsortsort.append(startsort[i])
    endsortsort.append(endsort[i])

# Unify coordinates

chrnu,startu,endu = unify_coordinates(chrnsortsort,startsortsort,endsortsort)

# Translate original peaks to unified coordinates

chrnu1,startu1,endu1 = overlap_peaks(chrnu,startu,endu,chrn1,start1,end1)
chrnu2,startu2,endu2 = overlap_peaks(chrnu,startu,endu,chrn2,start2,end2)

# Save outputs

folder = name1+'_'+name2+'_outputs'
command = 'mkdir '+folder
print command
os.system(command)

out1_tmp = open(folder+'/'+name1+'_tmp.txt','w')

for i in range(len(chrnu1)):
    out1_tmp.write("%s\t%s\t%s\n"%(chrnu1[i],startu1[i],endu1[i]))

out1_tmp.close()

command = 'sort '+folder+'/'+name1+'_tmp.txt | uniq > '+folder+'/'+name1+'_unified.bed'
print command
os.system(command)

out2_tmp = open(folder+'/'+name2+'_tmp.txt','w')

for i in range(len(chrnu2)):
    out2_tmp.write("%s\t%s\t%s\n"%(chrnu2[i],startu2[i],endu2[i]))

out2_tmp.close()

command = 'sort '+folder+'/'+name2+'_tmp.txt | uniq > '+folder+'/'+name2+'_unified.bed'
print command
os.system(command)

out = open(folder+'/'+name1+'_'+name2+'_unified.bed','w')

for i in range(len(chrnu)):
    out.write("%s\t%s\t%s\n"%(chrnu[i],startu[i],endu[i]))

out2_tmp.close()

command = 'rm '+folder+'/*tmp*'
print command
os.system(command)

