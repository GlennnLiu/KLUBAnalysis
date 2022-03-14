import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

class Binning:
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

def StrToFloat(binning):
	new_bins=[]
	for b in binning.split(", "):
		new_bins.append(float(b))
	new_bins=np.array(new_bins)
	return new_bins

def OldBins(f,ca,re):
	oldbins=[]
	h=f.Get("data_obs_"+ca+"_"+re+"_"+var)
	for i in range(h.GetNbinsX()+1):
		oldbins.append(round(h.GetBinLowEdge(i+1),4))
	oldbins=np.array(oldbins)
	return oldbins

def GetHist(f,pr,ca,re):
        fullname=pr+"_"+ca+"_"+re+"_"+var
        if not f.GetListOfKeys().Contains(fullname):
                print "*** WARNING: histo " , fullname , " not available"
                return -1
        else:
                return f.Get(fullname)

def GetTotalBkg(f,ca,re):
        for i,bkg in enumerate(backgrounds):
                fullname=bkg+"_"+ca+"_"+re+"_"+var
                if not f.GetListOfKeys().Contains(fullname):
                        print "*** WARNING: histo " , fullname , " not available"
                        return -1
                if i == 0:
                        h=f.Get(fullname).Clone()
                else:
                        h.Add(f.Get(fullname))
        return h

def EmptyBinProtection(f,ca,re,new_bins):
        #h_DY=GetHist(f,"DY",ca,re).Rebin(len(new_bins)-1,"",new_bins)
        #h_TT=GetHist("TT",re).Rebin(len(new_bins)-1,"",new_bins)
        h_bkg=GetTotalBkg(f,ca,re).Rebin(len(new_bins)-1,"",new_bins)
        for i in range(len(new_bins)-1):
                if h_bkg.GetBinContent(i+1) < 4*h_bkg.GetBinError(i+1):
                        return False
        return True

channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]
names=['1btag','2btag','boost','VBF']

mHH=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']

massPoints=[250,260,270,280,300,320,350,400,450,500,550,600,650,700,750,800,850,900,1000,1250,1500,1750,2000,2500,3000]

regions=['SR', 'SStight', 'OSinviso', 'SSinviso']

global backgrounds
backgrounds=["DY","TT","QCD","Others"]

bchmk=["300","400","600","850","1250","1250"]

global var
var='DNNoutSM_kl_1'

fin={}
fin2={}
fout={}

for ich,ch in enumerate(channels):
	fin[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}_14Apr2021/analyzedOutPlotter.root".format(ch),"read")
	fin2[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_06Apr2021/{0}/analyzedOutPlotter_backup.root".format(ch),"read")
	fout[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_10May2021/{0}/analyzedOutPlotter.root".format(ch),"recreate")

	if ch == 'TauTau':
		n=5
	else:
		n=6
	h={}
	for ica,ca in enumerate(categories):
		for i in range(n):
			if ca == 'sboostedLLMcut' and i == 0:
                                benchmark="320"
                        else:
                                benchmark=bchmk[i]
			#path='{0}_{1}{2}/'.format(ch,ca,mHH[i])
			BestBin=Binning('{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i]))
			newlimit=BestBin.limit
			oldlimit=10.*parseFile("/data_CMS/cms/liugeliang/HHbbtautau_limits/cards_{0}Legacy2018_06Apr2021/{1}{2}DNNoutSM_kl_1/combined_out/comb.Radion{3}.log".format(ch,ca,mHH[i],benchmark))
			
			newbins=StrToFloat(BestBin.binning)
			#print BestBin.nbins
			#print newbins
			oldbins=OldBins(fin2[ch],ca+mHH[i],'SR')
			#print oldbins
			print ch,ca+mHH[i]
			prot=EmptyBinProtection(fin[ch],ca+mHH[i],'SR',oldbins)
			if newlimit > oldlimit and prot:
			#if 0:
				print "Still old binning"
				for re in regions:
                                        h['data_obs'+ca+mHH[i]+re]=GetHist(fin2[ch],'data_obs',ca+mHH[i],re).Clone()
                                        #fout[ch].Write()
                                        for m in massPoints:
                                                h['Radion'+str(m)+ca+mHH[i]+re]=GetHist(fin2[ch],'Radion'+str(m),ca+mHH[i],re).Clone()
                                                #print re,m
                                                #fout[ch].Write()
                                        for bkg in backgrounds:
                                                if bkg == 'QCD' and re is not 'SR':
                                                        continue
                                                h[bkg+ca+mHH[i]+re]=GetHist(fin2[ch],bkg,ca+mHH[i],re).Clone()


			else:
				for re in regions:
                                	h['data_obs'+ca+mHH[i]+re]=GetHist(fin[ch],'data_obs',ca+mHH[i],re).Rebin(len(newbins)-1,"",newbins)
                                	#fout[ch].Write()
                                	for m in massPoints:
                                        	h['Radion'+str(m)+ca+mHH[i]+re]=GetHist(fin[ch],'Radion'+str(m),ca+mHH[i],re).Rebin(len(newbins)-1,"",newbins)
						#print re,m
                                        	#fout[ch].Write()
                                	for bkg in backgrounds:
                                        	if bkg == 'QCD' and re is not 'SR':
                                                	continue
                                        	h[bkg+ca+mHH[i]+re]=GetHist(fin[ch],bkg,ca+mHH[i],re).Rebin(len(newbins)-1,"",newbins)
                                        	#fout[ch].Write()

			#print ch,ca+mHH[i],'newlimit better?',newlimit < oldlimit,'EmptyBinProtection?',prot
	
	fout[ch].Write()	
	fin[ch].Close()
	fin2[ch].Close()
	fout[ch].Close()
