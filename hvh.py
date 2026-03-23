import argparse
import math
import random

import montyseqlib

parser = argparse.ArgumentParser()
parser.add_argument('iterations', type=int)
parser.add_argument('depth', type=int)
parser.add_argument('--err_rate', type=float, default=0.12, metavar='<float>',
	help='sequencing error rate [%(default).2f]')
parser.add_argument('--seed', type=int)
arg = parser.parse_args()
if arg.seed: random.seed(arg.seed)


hom_count = {}
het_count = {}
for _ in range(arg.iterations):
	# hom: mom is A (differences are sequencing error)
	hom = {'A':0, 'C':0, 'G':0, 'T':0}
	for i in range(arg.depth):
		if random.random() < arg.err_rate:
			nt = random.choice('CGT')
			hom[nt] += 1
		else:
			hom['A'] += 1
	vals = list(hom.values())
	vals.sort(reverse=True)
	sig = f'{vals[0]}.{vals[1]}.{vals[2]}.{vals[3]}'
	if sig not in hom_count: hom_count[sig] = 0
	hom_count[sig] += 1

	# het: mom is A, dad is T
	het = {'A':0, 'C':0, 'G':0, 'T':0}
	for i in range(arg.depth):
		if random.random() < 0.5:
			nt = 'A' # mom
			if random.random() < arg.err_rate: nt = random.choice('CGT')
		else:
			nt = 'T' # dad
			if random.random() < arg.err_rate: nt = random.choice('ACG')
		het[nt] += 1
	vals = list(het.values())
	vals.sort(reverse=True)
	sig = f'{vals[0]}.{vals[1]}.{vals[2]}.{vals[3]}'
	if sig not in het_count: het_count[sig] = 0
	het_count[sig] += 1


for sig in hom_count or het_count:
	hom = hom_count[sig] if sig in hom_count else 0
	het = het_count[sig] if sig in het_count else 0
	print(sig, hom, het, hom / (hom+het), sep='\t')
