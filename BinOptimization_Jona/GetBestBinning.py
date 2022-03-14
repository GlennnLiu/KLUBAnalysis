import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

class Binning:
	def __init__(self,stra,nbins,binning,limit,width):
		self.stra=stra
		self.nbins=nbins
		self.binning=binning
		self.limit=limit
		self.width=width

#	def __init__(self,filename):
#		if not os.path.exists(filename):
#			print "File not exists! ",filename
#			return -1
#		f=open(filename,"r")
#		self.stra=f.readline().split(": ")[1]
#		self.nbins=float(f.readline().split(": ")[1])
#		self.limit=float(f.readline().split(": ")[1])
#		self.width=float(f.readline().split(": ")[1])

	def Write(self,filename):
		f=open(filename,"w")
		f.truncate()
		f.write("Strategy: "+self.stra)
		f.write("\nNumber of bins: "+str(self.nbins))
		f.write("\nBinning scheme: \n")
		f.write(self.binning)
		f.write("\nExpected limit: "+str(self.limit))
		f.write("\n2 sigma width: "+str(self.width))
		f.close()
	
	def Compare(self,second):
		if self.limit < second.limit:
			return True
		elif self.limit > second.limit:
			return False
		else:
			if self.width < second.width:
				return True
			elif self.width > second.width:
				return False
			else:
				if self.nbins < second.nbins:
					return True
				elif self.nbins > second.nbins:
					return False
				else:
					return True

	def Replace(self,second):
		self.stra=second.stra
                self.nbins=second.nbins
                self.binning=second.binning
                self.limit=second.limit
                self.width=second.width

def parseFile(filename, CL='50.0', exp=True):
    f = open(filename)
    matches = []
    for line in f:
        search = ('Expected %s%%: r <'%CL)
        if not exp: search = 'Observed Limit: r <'

        if not search in line:
            continue
        val = line.replace(search, '')
        val = float(val)
        matches.append(val)

    if len(matches) == 0:
        print "did not find any expected in file: " , filename, 'CL=', CL, 'exp?=', exp
        return 1e10
    else:
        return matches[-1]


channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]

mHH=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']

year="2018"

var="DNNoutSM_kl_1"

bchmk=["300","400","600","850","1250","1250"]

strategies=['ConstSize','FlatS','FlatB','FlatSB']

nbins_test=range(1,55)

k=0

for ch in channels:
        if ch == 'TauTau':
                n=5
        else:
                n=6
        for ca in categories:
                for i in range(n):
			if ca == 'sboostedLLMcut' and i == 0:
				benchmark="320"
			else:
				benchmark=bchmk[i]
			path="{0}_{1}{2}/".format(ch,ca,mHH[i])
			BestBin=Binning(" ",100," ",1e10,1e10)
			for stra in strategies:
				for nbins in nbins_test:
					Input="{0}{1}_{2}/comb.Radion{3}.log".format(path,stra,str(nbins),benchmark)
					if not os.path.exists(Input):
						#print "Warning","{0}{1}_{2}/comb.Radion{3}.log".format(path,stra,str(nbins),benchmark)
						continue
					binning=open("{0}{1}_{2}/Binning.txt".format(path,stra,str(nbins))).readline()
					limit=10.*parseFile(Input)
					width=10*abs(parseFile(Input,CL="97.5")-parseFile(Input,CL=" 2.5"))
					temp=Binning(stra,nbins,binning,limit,width)
					if not BestBin.Compare(temp):
						BestBin.Replace(temp)
			print "{0}_{1}{2}: {3}_{4}".format(ch,ca,mHH[i],BestBin.stra,str(BestBin.nbins))
			BestBin.Write(path+"BestBin.txt")
			Bin=BestBin.binning
			print BestBin.nbins, len(Bin.split(','))-1
