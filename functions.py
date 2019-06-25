import numpy as np
import glob
import bz2

# Transposes a list of strings
def transpose(seqls):
	return [''.join([seqls[i][p]for i in xrange(len(seqls))]) for p in xrange(len(seqls[0]))]
	
# Calculates the entropy of a MSA column
# Example input: ['L','L', '.', 'I', 'L', 'I', 'I', '.', '.', 'V']
def calcEntropyPosition(ls):
	ls = convert(ls)
	ls = ls[ls>0]
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

# Converts a column of a MSA to a 20 long vector of percentages per amino-acid
# Amino-acid percentages in the following order: ARDNCEQGHILKMFPSTWYV
def convert(ls):
	ls = ''.join(ls)
	ls = ls.upper()
	ls = np.array([ls.count(i)for i in 'ARDNCEQGHILKMFPSTWYV'], dtype=np.float64)
	if ls.sum()>0:
		return ls/ls.sum()
	return ls

# Creating a dataset
def createDataset(hsspFiles, datasetName='dset'):
	dataset = np.zeros(int(10e6),20)
	ind = 0
	for hsspFile in hsspFiles:
		sequences = openHSSP(hsspFile)
		columns   = transpose(sequences)
		for column in comlumns:
			# only columns with more than 25% non-canonical AA  and gaps are added to dataset
			if sum([column.count(i)for i in 'ARDNCEQGHILKMFPSTWYV'])/float(len(column)) < .25:
				continue
			# Amino-acids to frequencies 
			column = convert(column)
			# Add column to dataset
			dataset[ind] = column
			ind += 1
	# Save dataset
	np.save(datasetName, dataset[:ind])
			
# Function that reads and decompresses hssp files and returns the MSA sequences
def openHSSP(hsspLocation):
	hsspid = hsspLocation.split('/')[-1]
	hsspid = hsspid.split('.')[0]
	file0  = open(hsspLocation)
	compressed = file0.read()
	file0.close()
	data = bz2.decompress(compressed)
	return parsHssp(data, hsspid)
