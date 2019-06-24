from matplotlib import pyplot as plt
from torch.autograd import Variable
import torch.nn.functional as F
import torch.optim as optim
import torch.nn as nn
import numpy as np
import argparse
import random
import torch
import glob
import time

# Module that learns one encoding and one decoding layer in a seperate  
# class since layers are trained in a "greedy pretrained" manner
class Auto(nn.Module):
	def __init__(self, inp, hid):
		super(Auto, self).__init__()
		self.fc1 = nn.Linear(inp, hid) # Encoding layer
		self.fc2 = nn.Linear(hid, inp) # Decoding layer
		self.bn = nn.BatchNorm1d(hid, affine=True) 
	def encode(self, x):
		x = self.bn(self.fc1(x))
		return F.sigmoid(x)
	def decode(self, x):
		x = self.fc2(x)
		return F.sigmoid(x)
	def forward(self, x):
		x1 = self.encode(x)
		x2 = self.decode(x1)
		return x1, x2

# The class that takes all the pretrained layers and combines them into a
# single model.  
class Total(nn.Module):
	def __init__(self, models):
		super(Total, self).__init__()
		self.models = models
	def forward(self, x):
		# Encode the input
		x = self.encode(x)
		# Reconstruct the encoding 
		x = self.decode(x)
		return x
	def encode(self, x):
		for i in xrange(len(self.models)):
			x = self.models[i].encode(x)
		return x
	def decode(self, x):
		for i in range(len(self.models))[::-1]:
			x = self.models[i].decode(x)
		return x
	# return the parameters of all layers for optimization
	def parameters(self):
		ls = [i for i in models[0].parameters()]
		for i in xrange(1, len(self.models)):
			ls += [i for i in models[i].parameters()]
		return ls
	def saveModels(self):
		[torch.save(self.models[index], 'aa2vec/model%i.pt' % (index)) \
							for index in xrange(len(self.models))]
	# Puts the total autoencoder in evaluation mode
	def eval(self):
		for i in self.models:
			i.eval()
	# Puts the total autoencoder in training mode
	def train(self):
		for i in self.models:
			i.train()

# Greedy pretraining of the individual layers
def trainLayers(models):
  newset = dset.copy()
  for index in xrange(len(models)):
	for x in xrange(max(0, index-1), index):
		newset = Variable(models[x].encode(newset).data)
	model = models[index]
	optimizer = optim.Adam(model.parameters())
	for epoch in xrange(epochs):
		model.train()
		rand = torch.randperm(newset.size()[0])
		t0 = time.time()
		for batch in xrange((newset.size()[0]-batchsize)/batchsize):
			ind = rand[batch*batchsize:(batch+1)*batchsize]
			optimizer.zero_grad()
			inp = newset[ind]
			hidden, output = model(inp)
			loss = F.binary_cross_entropy(output, inp) 
			loss.backward()
			optimizer.step()
		print 'Train %i : %.4f time : %.2f' % (epoch+1, loss.data[0], time.time()-t0)
		torch.save(model, 'aa2vec/model%i.pt' % (index))

# Function that finetunes the total autoencoder after pretraining
def trainModel():
	optimizer = optim.Adam(model.parameters())
	for epoch in xrange(epochs):
		rand = torch.randperm(dset.size()[0])
		t0 = time.time()
		model.train()
		for batch in xrange((dset.size()[0]-batchsize)/batchsize):
			ind = rand[batch*batchsize:(batch+1)*batchsize]
			optimizer.zero_grad()
			inp = dset[ind]
			output = torch.squeeze(model(inp))
			loss = F.binary_cross_entropy(output, inp)
			loss.backward()
			optimizer.step()
		print 'Train %i : %.4f time : %.2f' % (epoch+1, loss.data[0], time.time()-t0)
		torch.save(model, 'TotalModel.pt')
		model.saveModels()

# Transposes a list of strings
def transpose(seqls):
	return [''.join([seqls[i][p]for i in xrange(len(seqls))]) for p in xrange(len(seqls[0]))]
	
# Calculates the entropy of a MSA column
def calcEntropyPosition(ls):
	ls = ''.join(ls)
	ls = np.array([ls.count(i)for i in 'QWERTYIPASDFGHKLCVNM'], dtype=np.float64)
	ls = ls[ls>0]/ls.sum()
	return abs(sum((np.log(ls)*ls)/np.log(20)))

if __name__ == "__main__":
	# Load the HSSP dataset
	# Each row/sample represents a column of a MSA
	# There are 20 values per row representing the amino-acids percentages
	dset = np.load('dset.npy')
	# Sort the dataset, so it will learn about the amino-acid distributions
	# If not sorted, the autoencoder will learn amino-acid properties
	dset.sort() 
	# Prepare dataset for the neural network
	dset = Variable(torch.Tensor(dset))
	# Define the number of hidden layers and the number of neurons per layer
	auto = [20, 16, 10, 5, 2]
	# Create the layers in a list
	models = [Auto(auto[i], auto[i+1]) for i in xrange(0, len(auto)-1)]
	# the number of times the layers will run over the complete dataset
	epochs = 2
	# The number of samples the layers will see per trainingstep
	batchsize = 100						  
	
	# Greedy pretraining of the layers
	trainLayers(models)
	# Create the full autoencoder
	model = Total(models)
	# Finetune the complete autoencoder
	trainModel()
	
	
	
	
	
	

