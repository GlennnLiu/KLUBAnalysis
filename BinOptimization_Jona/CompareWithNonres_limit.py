import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

class Binning:
#        def __init__(self,stra,nbins,binning,limit,width):
#                self.stra=stra
#                self.nbins=nbins
#                self.binning=binning
#                self.limit=limit
#                self.width=width

	def __init__(self,filename):
		if not os.path.exists(filename):
			print "File not exists! ",filename
			return -1
		f=open(filename,"r")
		self.stra=f.readline().split(": ")[1]
		self.nbins=float(f.readline().split(": ")[1])
		f.readline()
		self.binning=f.readline()
		self.limit=float(f.readline().split(": ")[1])
		self.width=float(f.readline().split(": ")[1])

        def Write(self,filename):
                f=open(filename,"w")
                f.truncate()
                f.write("Strategy: "+self.stra)
                f.write("\nNumber of bins: "+str(self.nbins))
                f.write("\nBinning scheme: \n")
                f.write(binning)
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
names=['1btag','2btag','boost','VBF']

mHH=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']

bchmk=["300","400","600","850","1250","1250"]

for ch in channels:

	c1 = rt.TCanvas("c1", "c1", 650, 500)
	c1.SetFrameLineWidth(3)
	c1.SetBottomMargin (0.6)
	c1.SetRightMargin (0.05)
	c1.SetLeftMargin (0.15)
	c1.SetGridx()
	c1.SetGridy()
	
	pad1=rt.TPad("pad1","pad1",0,0.5,1,1,-1,0,0)
	pad2=rt.TPad("pad2","pad2",0,0.0,1,0.5,-1,0,0)
	pad1.SetBottomMargin(0.0)
	pad1.SetGridx()
	pad1.SetGridy()
	pad2.SetBottomMargin(0.5)
	pad2.SetTopMargin(0)
	pad2.SetGridx()
	pad2.SetGridy()
	pad1.Draw()
	pad2.Draw()

	if ch == 'TauTau':
                n=5
        else:
                n=6
	h_new=rt.TH1D("","",4*n,0,4*n)
        h_new.SetMarkerColor(rt.kBlue+2)
	h_new.SetMarkerStyle(24)
	h_new.SetStats(0)
	
	h_old=rt.TH1D("","",4*n,0,4*n)
        h_old.SetMarkerColor(rt.kAzure+10)
	h_old.SetMarkerStyle(24)
	h_old.SetStats(0)
	
	ratio=rt.TH1D("","",4*n,0,4*n)
	ratio.SetMarkerStyle(34)
	ratio.SetMarkerSize(1.6)
	ratio.SetStats(0)
	ratio.GetYaxis().SetTitleSize(0.11*3/5)
	ratio.GetXaxis().SetTitleSize(0.13*3/5)
	ratio.GetYaxis().SetLabelSize(0.1*3/5)
	ratio.GetXaxis().SetLabelSize(0.1*3/5)
	ratio.GetXaxis().SetLabelOffset(0.012)
	ratio.GetYaxis().SetTitleOffset(0.7)
	ratio.GetXaxis().SetTitleOffset(0.9)
	ratio.GetYaxis().SetTitle("#frac{Bin Optimized}{Previous}")
	ratio.SetMaximum(1.3)
	ratio.SetMinimum(0.8)
	
	hframe=rt.TH1D("","",4*n,0,4*n)
	hframe.GetYaxis().SetTitleSize(0.047*7/5)
	hframe.GetXaxis().SetTitleSize(0.055*7/5)
	hframe.GetYaxis().SetLabelSize(0.045*7/5)
	hframe.GetXaxis().SetLabelSize(0.055*7/5)
	hframe.GetXaxis().SetLabelOffset(0.012)
	hframe.GetYaxis().SetTitleOffset(0.7)
	hframe.GetXaxis().SetTitleOffset(1)
	hframe.GetYaxis().SetTitle("95% CL on #sigma #times #bf{#it{#Beta}}(S#rightarrowHH#rightarrow bb#tau#tau) [fb]")
	hframe.SetStats(0)
	hframe.SetMaximum(10000.)
	hframe.SetMinimum(1.1)
	
	k=1
        for j,ca in enumerate(categories):
		for i in range(n):
			if ca == 'sboostedLLMcut' and i == 0:
                                benchmark="320"
                        else:
                                benchmark=bchmk[i]
			
			fnew="{0}_{1}{2}/BestBin.txt".format(ch,ca,mHH[i])
			limit_new=Binning(fnew).limit
			
			fold="/data_CMS/cms/liugeliang/HHbbtautau_limits/cards_{0}Legacy2018_06Apr2021/{1}{2}DNNoutSM_kl_1/combined_out/comb.Radion{3}.log".format(ch,ca,mHH[i],benchmark)
			limit_old=10.*parseFile(fold)
			
			if limit_new < limit_old*0.1:
				limit_new=limit_new*100
			h_new.SetBinContent(k,limit_new)
			h_old.SetBinContent(k,limit_old)
			ratio.SetBinContent(k,limit_new/limit_old)
			ratio.GetXaxis().SetBinLabel(k,names[j]+"_"+mHH[i])
			k=k+1
	
	pad1.cd()
	hframe.Draw()
	h_new.Draw("Psame")
	h_old.Draw("Psame")
	pad1.SetLogy(1)
	c1.Update()
	
	pad2.cd()
	ratio.LabelsOption("v")
	ratio.Draw("P")
	#pad2.SetLogy(1)
	c1.Update()
	
	c1.SaveAs("Output/Comparison_{0}.pdf".format(ch),"pdf")
