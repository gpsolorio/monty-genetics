import argparse
import random

import montyseqlib

parser = argparse.ArgumentParser()
parser.add_argument('length', type=int)
parser.add_argument('depth', type=int)
parser.add_argument('snps', type=int)
parser.add_argument('--seed', type=int)
parser.add_argument('--sub_rate', type=float, default=0.12, metavar='<float>', help='Substitution Rate [%(default).2f]')
arg = parser.parse_args()

if arg.seed: random.seed(arg.seed)

# create mom and dad haplotypes
mom_seq = ''.join(random.choices('ACGT', k=arg.length))

dad_seq = list(mom_seq)
not_mom = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
for i in range(arg.snps):
    dad_seq[i] = not_mom[dad_seq[i]]
dad_seq = ''.join(dad_seq)

# generate reads from mom and dad
for i, read in enumerate(montyseqlib.create_reads(mom_seq, arg.depth, sub_rate=arg.sub_rate)):
    print(f'>m.{i}')
    print(read)

for i, read in enumerate(montyseqlib.create_reads(dad_seq, arg.depth, sub_rate=arg.sub_rate)):
    print(f'>d.{i}')
    print(read)