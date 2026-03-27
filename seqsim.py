import argparse
import math
import random
import sys

import montyseqlib

def aseq(length, variations, padding=' ', symbol='*'):
	seq = [padding] * length
	for pos, nt in variations: seq[pos] = symbol
	return ''.join(seq)

parser = argparse.ArgumentParser()
parser.add_argument('chromlen', type=float, metavar='<float>',
	help='length of the chromosomal region (e.g. 1e5)')
parser.add_argument('readlen', type=float, metavar='<float>',
	help='length of sequencing reads (e.g. 1e3)')
parser.add_argument('depth', type=float, metavar='<float>',
	help='sequencing depth (e.g. 10 means 10x)')
parser.add_argument('--output', '-o', metavar='<file>',
	help='output file [default is stdout]')
parser.add_argument('--snp_rate', '-s', type=float, default=0.001,
	metavar='<float>', help='sequencing error rate [%(default).2f]')
parser.add_argument('--err_rate', '-e', type=float, default=0.1,
	metavar='<float>', help='sequencing error rate [%(default).2f]')
parser.add_argument('--randomseed', '-r', type=int, metavar='<int>',
	help='use a specific random seed [default is random]')
arg = parser.parse_args()
if arg.randomseed: random.seed(arg.randomseed)
arg.chromlen = int(arg.chromlen)
arg.readlen  = int(arg.readlen)

# create parental genotypes
src = ''.join(random.choices('ACGT', k=arg.chromlen)) # reference seq
mom, vmom = montyseqlib.create_variant(src, arg.snp_rate)
dad, vdad = montyseqlib.create_variant(src, arg.snp_rate)

# create sequencing reads
reads = {} # organized by starting offset
readn = int(arg.depth * arg.chromlen / arg.readlen)
for _ in range(readn):
	off = random.randint(0, arg.chromlen - arg.readlen)
	if off not in reads: reads[off] = []
	if random.random() < 0.5:
		reg = mom[off:off+arg.readlen]
		s, v = montyseqlib.create_variant(reg, arg.err_rate, lower=True)
		reads[off].append( ('m', s) )
	else:
		reg = dad[off:off+arg.readlen]
		s, v = montyseqlib.create_variant(reg, arg.err_rate, lower=True)
		reads[off].append( ('d', s) )

# display alignment
print('R', src)
print('M', mom)
print('v', aseq(arg.chromlen, vmom, symbol='m'))
print('D', dad)
print('v', aseq(arg.chromlen, vdad, symbol='d'))
for off, aligns in sorted(reads.items()):
	for parent, seq in aligns:
		head = ' ' * off
		print(parent, ' ', head, seq, sep='')

# decode alignments...
# for each column of the alignment
# count the letters and sort into decreasing order
# look up p(het) in heterozygous.lut
# if a column is het, label the reads as belonging to mom or dad
# compute accuracy of labeling