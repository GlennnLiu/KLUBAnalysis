import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

class Binning1:
#        def __init__(self,stra,nbins,binning,limit,width):
#                self.stra=stra
#                self.nbins=nbins
#                self.binning=binning
#                self.limit=limit
#                self.width=width

	def __init__(self,filename):
		if not os.path.exists(filename):
			print "File not exists! ",filename
			return None
		f=open(filename,"r")
		self.stra=f.readline().split(": ")[1]
		self.nbins=int(f.readline().split(": ")[1])
		f.readline()
		self.binning=f.readline()
		self.limit=float(f.readline().split(": ")[1])
		self.width=float(f.readline().split(": ")[1])
		f.close()

class Binning2:
	def __init__(self,filename):
                if not os.path.exists(filename):
                        print "File not exists! ",filename
                        self.nbins=55
			self.binning=""
			self.limit=1e7
			return None
                f=open(filename,"r")
                #self.stra=f.readline().split(": ")[1]
                self.nbins=int(f.readline().split(": ")[1])
                f.readline()
                self.binning=f.readline()
                self.limit=float(f.readline().split(": ")[1])
		f.close()
                #self.width=float(f.readline().split(": ")[1])

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

#categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","VBFloose"]
categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]
#names=['1btag','2btag','VBF']
names=['1btag','2btag','boost','VBF']

mHH1=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']
mHH2=['mHH250_625', 'mHH625_775', 'mHH775_1100','mHH1100_3500','InvDR']

bchmk1=["300","400","600","850","1250","1250"]
bchmk2=["450","700","900","1250","1250"]
bchmk={}
for i in range(len(mHH1)):
        bchmk[mHH1[i]]=bchmk1[i]
for i in range(len(mHH2)):
        bchmk[mHH2[i]]=bchmk2[i]

pt = rt.TPaveText(0.1663218-0.02,0.886316,0.3045977-0.02,0.978947,"brNDC")
pt.SetBorderSize(0)
pt.SetTextAlign(12)
pt.SetTextFont(62)
pt.SetTextSize(0.075)
pt.SetFillColor(0)
pt.SetFillStyle(0)
pt.AddText("CMS #font[52]{Internal}" )

pt2 = rt.TPaveText(0.736111,0.9066667,0.847222,0.954641,"brNDC")
pt2.SetBorderSize(0)
pt2.SetFillColor(0)
pt2.SetTextSize(0.06)
pt2.SetTextFont(42)
pt2.SetFillStyle(0)
pt2.AddText("2018 - 59.7 fb^{-1} (13 TeV)")


for ch in channels:

	if ch == 'TauTau':
		n=19
	else:
		n=23
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
	
        pt4 = rt.TPaveText(0.4819196+0.036-0.05,0.7780357+0.015+0.02,0.9008929+0.036-0.05,0.8675595+0.015,"brNDC")
        pt4.SetTextAlign(12)
        pt4.SetFillColor(rt.kWhite)
        pt4.SetFillStyle(1001)
        pt4.SetTextFont(42)
        pt4.SetTextSize(0.075)
        pt4.SetBorderSize(0)
        pt4.SetTextAlign(32)
	if ch == 'ETau':
		pt4.AddText('bb e#tau_{h}')
	elif ch == 'MuTau':
		pt4.AddText('bb #mu#tau_{h}')
	elif ch == 'TauTau':
		pt4.AddText('bb #tau_{h}#tau_{h}')
		
	h_new=rt.TH1D("","",n,0,n)
        h_new.SetMarkerColor(rt.kBlue+2)
	h_new.SetMarkerStyle(24)
	h_new.SetStats(0)
	
	h_old=rt.TH1D("","",n,0,n)
        h_old.SetMarkerColor(rt.kAzure+10)
	h_old.SetMarkerStyle(24)
	h_old.SetStats(0)
	
	ratio=rt.TH1D("","",n,0,n)
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
	ratio.GetYaxis().SetTitle("#frac{tot_bkg_thr=1}{tot_bkg_thr=0.18}")
	#ratio.SetMaximum(1.2)
	#ratio.SetMinimum(0.8)
	
	hframe=rt.TH1D("","",n,0,n)
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
		if ch == 'TauTau':
                        if ca == 'sboostedLLMcut':
                                mHH=['mHH250_625', 'mHH625_775', 'mHH775_1100','mHH1100_3500']
                        else:
                                mHH=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500']
                else:
                        if ca == 'sboostedLLMcut':
                                mHH=['mHH250_625', 'mHH625_775', 'mHH775_1100','mHH1100_3500','InvDR']
                        else:
                                mHH=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']
		for m in mHH:
			#fold="../BinOptimization_Jona/{0}_{1}{2}/BestBin.txt".format(ch,ca,mHH[i])
			fold1='{0}_{1}{2}/BestBin.txt'.format(ch,ca,m)
			fnew2='/data_CMS/cms/liugeliang/HHbbtautau_BinOpt/Kons/backup8/{0}_{1}{2}/BestBin.txt'.format(ch,ca,m)
			fnew1='backup9/{0}_{1}{2}/BestBin.txt'.format(ch,ca,m)
                        fold2='/data_CMS/cms/liugeliang/HHbbtautau_BinOpt/Kons/backup7/{0}_{1}{2}/BestBin.txt'.format(ch,ca,m)
			#fnew2="backup5/{0}_{1}{2}/BestBin.txt".format(ch,ca,mHH[i])
			#fnew2="backup/{0}_{1}{2}/BestBin.txt".format(ch,ca,mHH[i])
			#limit_Jona=Binning1(fold).limit
			#limit_new=min(Binning2(fnew1).limit,Binning2(fnew2).limit)
			limit_new1=Binning2(fnew1).limit
			limit_new2=Binning2(fnew2).limit
			limit_old1=Binning2(fold1).limit
                        limit_old2=Binning2(fold2).limit
			limit_new=min(limit_new1,limit_new2)
			limit_old=min(limit_old1,limit_old2)

			#limit_old=Binning2(fnew2).limit
			#temp=open('{0}_{1}{2}/Load.txt'.format(ch,ca,mHH[i]),'r')
                        #temp.readline()
                        #same='e' in temp.readline()
                        #same=same and 'e' in temp.readline()
                        #same=same and not 'e' in temp.readline()
                        #print same
                        #if limit_new > limit_old and same:
			#	limit_new=limit_old
			
			#fold="/data_CMS/cms/liugeliang/HHbbtautau_limits/cards_{0}Legacy2018_06Apr2021/{1}{2}DNNoutSM_kl_1/combined_out/comb.Radion{3}.log".format(ch,ca,m,bchmk[m])
			#limit_veryold=10.*parseFile(fold)
			
			#if limit_new < limit_old*0.1:
			#	limit_new=limit_new*100
			h_new.SetBinContent(k,limit_new)
			h_old.SetBinContent(k,limit_old)
			ratio.SetBinContent(k,limit_new/limit_old)
			ratio.GetXaxis().SetBinLabel(k,names[j]+"_"+m)
			k=k+1
	
	legend=rt.TLegend(0.1,0.75,0.45,0.9)
	legend.SetBorderSize(0)
	legend.SetFillStyle(0)
	legend.AddEntry(h_new,'tot_bkg_thr=1','P')
	legend.AddEntry(h_old,'tot_bkg_thr=0.18','P')
	
	pad1.cd()
	hframe.Draw()
	h_new.Draw("Psame")
	h_old.Draw("Psame")
	pt.Draw()
	pt2.Draw()
	pt4.Draw()
	pad1.SetLogy(1)
	legend.Draw('Same')
	c1.Update()
	
	pad2.cd()
	ratio.LabelsOption("v")
	ratio.Draw("P")
	#pad2.SetLogy(1)
	c1.Update()
	
	c1.SaveAs("Output/Comparison_{0}.pdf".format(ch),"pdf")
