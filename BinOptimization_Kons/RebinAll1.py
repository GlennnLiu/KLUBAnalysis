import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array

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

def GetBinning(template):
	nPoints=template.GetNbinsX()
	old_bins=[template.GetBinLowEdge(1)]
	for ibin in range(1,nPoints+1):
        	old_bins.append(template.GetBinLowEdge(ibin+1))
	old_bins=np.array(old_bins)
	return old_bins


channels=['ETau','MuTau','TauTau']

categories=["s1b1jresolvedMcut","s2b0jresolvedMcut","sboostedLLMcut","VBFloose"]
names=['1btag','2btag','boost','VBF']

mHH1=['mHH250_335','mHH335_475','mHH475_725','mHH725_1100','mHH1100_3500','InvDR']
mHH2=['mHH250_625', 'mHH625_775', 'mHH775_1100','mHH1100_3500','InvDR']

massPoints=[250,260,270,280,300,320,350,400,450,500,550,600,650,700,750,800,850,900,1000,1250,1500,1750,2000,2500,3000]

regions=['SR', 'SStight', 'OSinviso', 'SSinviso']

global backgrounds
backgrounds=["DY","TT","QCD","Others"]

bchmk1=["300","400","600","850","1250","1250"]
bchmk2=["450","700","900","1250","1250"]
bchmk={}
for i in range(len(mHH1)):
        bchmk[mHH1[i]]=bchmk1[i]
for i in range(len(mHH2)):
        bchmk[mHH2[i]]=bchmk2[i]

global var
var='DNNoutSM_kl_1'


for ch in channels:
	fin=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}_12Oct2021/analyzedOutPlotter.root".format(ch),"read")
	fbin=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_24May2021/{0}/analyzedOutPlotter.root".format(ch),"read")
	os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_12Oct2021")
	os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_12Oct2021/{0}".format(ch))
	fout=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_12Oct2021/{0}/analyzedOutPlotter.root".format(ch),"recreate")
	
	h={}
	
	for ca in categories:
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
		
		for mhh in mHH:
			print ch,ca,mhh
			temp=GetHist(fbin,'DY',ca+mhh,'SR')
			newbins=GetBinning(temp)
			for re in regions:
				h['data_obs'+ca+mhh+re]=GetHist(fin,'data_obs',ca+mhh,re).Rebin(len(newbins)-1,"",newbins)
                                for m in massPoints:
                                	h['Radion'+str(m)+ca+mhh+re]=GetHist(fin,'Radion'+str(m),ca+mhh,re).Rebin(len(newbins)-1,"",newbins)
                                for bkg in backgrounds:
                                        if bkg == 'QCD' and re is not 'SR':
                                                continue
                                        h[bkg+ca+mhh+re]=GetHist(fin,bkg,ca+mhh,re).Rebin(len(newbins)-1,"",newbins)
                        h['totalBkg']=GetTotalBkg(fin,ca+mhh,'SR').Rebin(len(newbins)-1,"",newbins)
	fout.Write()
	fin.Close()
	fbin.Close()
	fout.Close()
