import numpy as np
import glob
import bz2

# Transposes a list of strings
def transpose(seqls):
	return [''.join([seqls[i][p]for i in xrange(len(seqls))]) for p in xrange(len(seqls[0]))]
	
# Calculates the entropy of a MSA column
# Example input: ['L','L', '.', 'I', 'L', 'I', 'I', '.', '.', 'V']
def calcEntropyPosition(ls):
	ls = ''.join(ls)
	ls = np.array([ls.count(i)for i in 'QWERTYIPASDFGHKLCVNM'], dtype=np.float64)
	ls = ls[ls>0]/ls.sum()
	return abs(sum((np.log(ls)*ls)/np.log(20)))

# Parses a HSSP file, extracts the MSA sequences and 
# returns them in a list. 
def parsHssp(string, hsspid):
	# Find the start of the MSA
	msaStart = string.find('\n%s' % (hsspid))
	string = string[msaStart:]
	string = string.split('\n')
	sequences = []
	pos = 0
	for line in string[1:-2]:
		if line == '': # Empty line
			continue
		if line.startswith('#=GC'): 
			pos = 0 
		else:
			try:
				line = line.split(' ')[-1]
				sequences[pos] += line
			except: 
				sequences += [line]
			pos += 1
	return sequences

# Function that reads and decompresses hssp files and returns the MSA sequences
def openHSSP(hsspLocation):
	hsspid = hsspLocation.split('/')[-1]
	hsspid = hsspid.split('.')[0]
	file0  = open(hsspLocation)
	compressed = file0.read()
	file0.close()
	data = bz2.decompress(compressed)
	return parsHssp(data, hsspid)

a= glob.glob('hg-hssp/*')
print len(a)
print a[0]
data = openHSSP(a[500])
data = transpose(data)
print calcEntropyPosition(data[0])

