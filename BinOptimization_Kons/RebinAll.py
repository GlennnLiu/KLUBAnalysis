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
	h.SetName('totalBkg_'+ca+"_"+re+"_"+var)
        h.SetTitle('totalBkg_'+ca+"_"+re+"_"+var)
        return h

def EmptyBinProtection(f,ca,re,new_bins):
        h_DY=GetHist(f,"DY",ca,re).Rebin(len(new_bins)-1,"",new_bins)
        h_TT=GetHist(f,"TT",ca,re).Rebin(len(new_bins)-1,"",new_bins)
        h_bkg=GetTotalBkg(f,ca,re).Rebin(len(new_bins)-1,"",new_bins)
        for i in range(len(new_bins)-1):
                if h_bkg.GetBinContent(i+1) < 0.18 or h_DY.GetBinContent(i+1) < 1e-7 or h_TT.GetBinContent(i+1) < 1e-7:#2*h_bkg.GetBinError(i+1):
                        return False
        return True

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

fin={}
fin2={}
fout={}

tag='23May2021'

os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}".format(tag))

for ich,ch in enumerate(channels):
	fin[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}_23May2021/analyzedOutPlotter.root".format(ch),"read")
	#fin2[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_06Apr2021/{0}/analyzedOutPlotter.root".format(ch),"read")
	os.system("mkdir /data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}/{1}".format(tag,ch))
	fout[ch]=rt.TFile("/data_CMS/cms/liugeliang/HHbbtautau_histos/Legacy2018_{0}/{1}/analyzedOutPlotter.root".format(tag,ch),"recreate")

	
	h={}
	for ica,ca in enumerate(categories):
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
                for i in range(len(mHH)):
			#path='{0}_{1}{2}/'.format(ch,ca,mHH[i])
			fold1='{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i])
                        fnew2='/data_CMS/cms/liugeliang/HHbbtautau_BinOpt/Kons/backup8/{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i])
                        fnew1='backup9/{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i])
                        fold2='/data_CMS/cms/liugeliang/HHbbtautau_BinOpt/Kons/backup7/{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i])
			
			new1=Binning2(fnew1)
                        new2=Binning2(fnew2)
                        old1=Binning2(fold1)
                        old2=Binning2(fold2)
			
			#BestBin1=Binning2('{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i]))
			#BestBin2=Binning2('backup5/{0}_{1}{2}/BestBin.txt'.format(ch,ca,mHH[i]))
			#temp=open('{0}_{1}{2}/Load.txt'.format(ch,ca,mHH[i]),'r')
			#temp.readline()
			#same='e' in temp.readline()
			#same=same and 'e' in temp.readline()
			#same=same and not 'e' in temp.readline()
			#print same
			#same=False
			#if BestBin1.limit > BestBin2.limit and same:
                        #        BestBin=BestBin2
		#		print "old"
                #        else:
                #                BestBin=BestBin1
		#	newlimit=BestBin.limit
		#	oldlimit=10.*parseFile("/data_CMS/cms/liugeliang/HHbbtautau_limits/cards_{0}Legacy2018_06Apr2021/{1}{2}DNNoutSM_kl_1/combined_out/comb.Radion{3}.log".format(ch,ca,mHH[i],benchmark))
			if old1.target>old2.target:
				BestBin=old1
			else:
				BestBin=old2
			
			newbins=StrToFloat(BestBin.binning)
			#print BestBin.nbins
			#print newbins
			#oldbins=OldBins(fin2[ch],ca+mHH[i],'SR')
			#print oldbins
			print ch,ca+mHH[i]
			#prot=EmptyBinProtection(fin[ch],ca+mHH[i],'SR',oldbins)
			prot=False
			if prot:
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
				h['totalBkg']=GetTotalBkg(fin[ch],ca+mHH[i],'SR').Rebin(len(newbins)-1,"",newbins)
                                        	#fout[ch].Write()

			#print ch,ca+mHH[i],'newlimit better?',newlimit < oldlimit,'EmptyBinProtection?',prot
	
	fout[ch].Write()	
	fin[ch].Close()
	#fin2[ch].Close()
	fout[ch].Close()
