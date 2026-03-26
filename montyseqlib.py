import random
import sys
import gzip

def create_variant(parent, var_rate=0.001, lower=False):
	'''Create a duplicate sequence with random variants'''
	if lower:
		notA = 'cgt'
		notC = 'agt'
		notG = 'act'
		notT = 'acg'
	else:
		notA = 'CGT'
		notC = 'AGT'
		notG = 'ACT'
		notT = 'ACG'
	seq = []
	variants = []
	for i, nt in enumerate(parent):
		if random.random() < var_rate:
			if   nt == 'A': nt = random.choice(notA)
			elif nt == 'C': nt = random.choice(notC)
			elif nt == 'G': nt = random.choice(notG)
			elif nt == 'T': nt = random.choice(notT)
			else: sys.exit('Error')
			variants.append( (i, nt) )
		seq.append(nt)
	return ''.join(seq), variants

def create_reads(parent, num_reads, sub_rate=0.1, lower=False):
	'''Generate reads with random sequencing errors'''
	for _ in range(num_reads):
		s, v = create_variant(parent, var_rate=sub_rate, lower=lower)
		yield s

##################
## File Readers ##
##################

def getfp(filename):
	"""Returns a file pointer for reading based on file name"""
	if   filename.endswith('.gz'): return gzip.open(filename, 'rt')
	elif filename == '-':          return sys.stdin
	else:                          return open(filename)

def readfasta(filename):
	"""Simple fasta file iterator: yields defline, seq"""
	name = None
	seqs = []
	fp = getfp(filename)
	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield name, seq
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield name, ''.join(seqs)
	fp.close()
