import sys

import montyseqlib

seqs = []
for defline, seq in montyseqlib.readfasta(sys.argv[1]):
    seqs.append(seq)

for i in range(len(seqs[0])):
    counts = {'A':0, 'C':0, 'G':0, 'T':0}
    for seq in seqs:
        counts[seq[i]]+=1
    print(counts)