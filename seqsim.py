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

src = ''.join(random.choices('ACGT', k=arg.chromlen)) # reference seq
mom, vmom = montyseqlib.create_variant(src, arg.snp_rate)
dad, vdad = montyseqlib.create_variant(src, arg.snp_rate)

reads = []
readn = int(arg.depth * arg.chromlen / arg.readlen)
for _ in range(readn):
	i = random.randint(0, arg.chromlen - arg.readlen)
	if random.random() < 0.5:
		reg = mom[i:i+arg.readlen]
		s, v = montyseqlib.create_variant(reg, arg.err_rate, lower=True)
		reads.append( (i, 'm', s ) )
	else:
		reg = dad[i:i+arg.readlen]
		s, v = montyseqlib.create_variant(reg, arg.err_rate, lower=True)
		reads.append( (i, 'd', s) )
reads = sorted(reads, key=lambda x: x[0])

print(' ', src)
print(' ', mom)
print(' ', aseq(arg.chromlen, vmom, symbol='m'))
print(' ', dad)
print(' ', aseq(arg.chromlen, vdad, symbol='d'))

for offset, parent, seq in reads:
	head = ' ' * offset
	print(parent, ' ', head, seq, sep='')

