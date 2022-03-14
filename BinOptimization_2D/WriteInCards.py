import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]

names=['1btag','2btag','boost','VBF']

massPoints=[250,260,270,280,300,320,350,400,450,500,550,600,650,700,750,800,850,900,1000,1250,1500,1750,2000,2500,3000]

regions=['SR', 'SStight', 'OSinviso', 'SSinviso']

backgrounds=["DY","TT","QCD","Others"]


for ch in channels:
	print ch
        selections=[]
        for ca in categories:
                if ch == 'TauTau':
                        folder='{0}_{1}'.format(ch,ca)
                else:
                        folder='{0}_{1}DR'.format(ch,ca)
                f=open('BinOptimizer_mHH_Output/{0}/BestBin.txt'.format(folder),'r')
                binning=np.array([int(float(b)) for b in f.readline().split(', ')])
                for i in range(len(binning)-1):
                        sel='{0}mHH{1}_{2}'.format(ca,str(binning[i]),str(binning[i+1]))
                        selections.append(sel)
                if ch != 'TauTau':
                        selections.append('{0}InvDR'.format(ca))
	
	print '"'+'","'.join(selections)+'"'
