################################################################
#                                         
#Calculate GC content                     
#Usage: python getGC_content.py fasta.file > GC_content.output 
################################################################

import sys

fastaFile = open(sys.argv[1], 'r')

gc = 0
at = 0
unknown = 0

for line in fastaFile:
    if line.startswith('>'):
        if(gc + at) > 0:
            total = gc + at
            total2 = total + unknown
            percentage = round((float(gc / total2))*100, 3 )
            print (seq_id,percentage)
        seq_id = line.strip()[1:]

        #reset counts
        gc = 0
        at = 0
        unknown = 0

    else:
        nuc_str = list(line.strip())
        for n in nuc_str:
            if n == 'G' or n == 'g' or n== 'C' or n == 'c':
                gc += 1
            elif n== 'A' or n =='a' or n == 'T' or n == 't':
                at += 1
            elif n== 'N' or n=='n':
                unknown += 1

#calculates for last sequence in file
total = gc +at
total2 = total + unknown
percent = round((float(gc / total2))*100, 3 ) 
print (seq_id,percent)

