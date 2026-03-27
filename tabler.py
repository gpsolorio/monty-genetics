import argparse
import glob
import sys

def custom_sort(sig):
	a, b, c, d = sig.split('.')
	a = int(a)
	b = int(b)
	c = int(c)
	d = int(d)
	return a + b + c + d + 0.1 * a + 0.01 * b + 0.001 * c + 0.001 * d

parser = argparse.ArgumentParser()
parser.add_argument('table_dir')
parser.add_argument('--minobs', type=int, metavar='<int>', default=5000,
	help='minimum number of observations [%(default)i]')
arg = parser.parse_args()

data = {}
for file in glob.glob(f'{arg.table_dir}/*'):
	with open(file) as fp:
		header = next(fp)
		for line in fp:
			sig, hom, het, phom = line.split()
			hom = int(hom)
			het = int(het)
			if max(het, hom) < arg.minobs:
				print(f'skiping {sig} {hom} {het} due to low counts', 
					file=sys.stderr)
				continue
			data[sig] = het / (hom + het)

ordered = list(data.keys())
ordered.sort(key=custom_sort)
for sig in ordered: print(sig, data[sig])