import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array
import argparse

channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]

def ArrayToStr(a):
    return ', '.join([ str(x) for x in a ])


for ch in channels:
	print ch
	for ca in categories:
		if ch == 'TauTau':
			folder='{0}_{1}'.format(ch,ca)
		else:
			folder='{0}_{1}DR'.format(ch,ca)
		f=open('BinOptimizer_mHH_Output/{0}/BestBin.txt'.format(folder),'r')
		binning=np.array([int(float(b)) for b in f.readline().split(', ')])
		bchmks=np.array([b for b in f.readline().split(', ')])
		selections=[]
		criteria=[]
		for i in range(len(bchmks)):
			sel='{0}mHH{1}_{2}'.format(ca,str(binning[i]),str(binning[i+1]))
			selections.append(sel)
			if ch == 'TauTau':
				cri='{0}, (HHKin_mass_raw >= {1} && HHKin_mass_raw < {2})'.format(ca,str(binning[i]),str(binning[i+1]))
			else:
				cri='{0}, DR, (HHKin_mass_raw >= {1} && HHKin_mass_raw < {2})'.format(ca,str(binning[i]),str(binning[i+1]))
			criteria.append(cri)
		if ch != 'TauTau':
			selections.append('{0}InvDR'.format(ca))
			criteria.append('{0}, InvDR'.format(ca))
		print '","'.join(selections)
		#for i in range(len(selections)):
			#print '"{0}",'.format(selections[i])
			#print '{0} = {1}'.format(selections[i],criteria[i])
		#print ArrayToStr(selections)
		#print '\n'
	#print '\n'
