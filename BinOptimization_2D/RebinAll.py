import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

class Binning2:
        def __init__(self,filename):
                if not os.path.exists(filename):
                        print "File not exists! ",filename
                        self.nbins=55
                        self.binning=""
                        self.limit=1e7
                        self.target=-(self.nbins+1e3*self.limit)
                        return None
                f=open(filename,"r")
                #self.stra=f.readline().split(": ")[1]
                self.nbins=int(f.readline().split(": ")[1])
                f.readline()
                self.binning=f.readline()
                self.limit=float(f.readline().split(": ")[1])
                self.target=-(self.nbins+1e3*self.limit)
                f.close()
                #self.width=float(f.readline().split(": ")[1])


def StrToFloat(binning):
        new_bins=[]
        for b in binning.split(", "):
                new_bins.append(float(b))
        new_bins=np.array(new_bins)
        return new_bins

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
        h.SetName('totalBkg_'+ca+"_"+re+"_"+var)
        h.SetTitle('totalBkg_'+ca+"_"+re+"_"+var)
        return h



channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]

names=['1btag','2btag','boost','VBF']

massPoints=[250,260,270,280,300,320,350,400,450,500,550,600,650,700,750,800,850,900,1000,1250,1500,1750,2000,2500,3000]

regions=['SR', 'SStight', 'OSinviso', 'SSinviso']

backgrounds=["DY","TT","QCD","Others"]

var='DNNoutSM_kl_1'

tag='09Jun2021'

os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}".format(tag))

fin={}
fin2={}
fout={}

for ch in channels:
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

	fin[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}_01Jun2021/analyzedOutPlotter.root".format(ch),"read")
	os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}/{1}".format(tag,ch))
        fout[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}/{1}/analyzedOutPlotter.root".format(tag,ch),"recreate")
	
	h={}
	for sel in selections:
		fnew='{0}_{1}/BestBin.txt'.format(ch,sel)
		bnew=Binning2(fnew)
		newbins=StrToFloat(bnew.binning)
		print ch,sel
		for re in regions:
			h['data_obs'+sel+re]=GetHist(fin[ch],'data_obs',sel,re).Rebin(len(newbins)-1,"",newbins)
                        for m in massPoints:
                        	h['Radion'+str(m)+sel+re]=GetHist(fin[ch],'Radion'+str(m),sel,re).Rebin(len(newbins)-1,"",newbins)
                	for bkg in backgrounds:
                                if bkg == 'QCD' and re is not 'SR':
                                	continue
                                h[bkg+sel+re]=GetHist(fin[ch],bkg,sel,re).Rebin(len(newbins)-1,"",newbins)
                h['totalBkg'+sel]=GetTotalBkg(fin[ch],sel,'SR').Rebin(len(newbins)-1,"",newbins)
	
	fout[ch].Write()
        fin[ch].Close()
	fout[ch].Close()

