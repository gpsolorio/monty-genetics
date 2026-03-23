import random
import sys

def create_reads(parent, num_reads, sub_rate=0.1):
    '''Generate reads with random sequencing errors'''
    for _ in range(num_reads):
        read = []
        for nt in parent:
            if random.random() > sub_rate:
                read.append(nt)
            else:
                if nt == 'A': read.append(random.choice('CGT'))
                elif nt == 'C': read.append(random.choice('AGT'))
                elif nt == 'G': read.append(random.choice('ACT'))
                elif nt == 'T': read.append(random.choice('AGC'))
                else: sys.exit('Error')
        yield ''.join(read)

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
