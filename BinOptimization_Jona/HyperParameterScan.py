import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array
import argparse
#from ConfigReader import *

def GetHist(pr,re):
	fullname=pr+"_"+ca+"_"+re+"_"+var
	if not fin.GetListOfKeys().Contains(fullname):
		print "*** WARNING: histo " , fullname , " not available"
		return -1
	else:
		return fin.Get(fullname)

def GetTotalBkg(re):
	for i,bkg in enumerate(backgrounds):
		fullname=bkg+"_"+ca+"_"+re+"_"+var
		if not fin.GetListOfKeys().Contains(fullname):
			print "*** WARNING: histo " , fullname , " not available"
                	return -1
		if i == 0:
			h=fin.Get(fullname).Clone()
		else:
			h.Add(fin.Get(fullname))
	return h

def ShiftToOldBins(d):
	if d < var_min or d > var_max:
		print "Out of range!!! ",d
		return -991
	close=np.digitize(np.array([d]),old_bins)[0]
	if abs(d-old_bins[close%len(old_bins)]) > abs(d-old_bins[close-1]):
		return old_bins[close-1]
	else:
		return old_bins[close%len(old_bins)]

def NewBin(sig,bkg,stra,nbins):
	if stra == 'ConstSize':
		new_bins=np.array([x for x in np.linspace(num=nbins+1,start=var_min,stop=var_max)])
		#make sure each edge in new binning is in old binning, otherwise modify it to be the closest edge in old binning
		for i in range(len(new_bins)):
			new_bins[i]=ShiftToOldBins(new_bins[i])
		return new_bins
	if stra == 'FlatS':
		template=sig.Clone()
	elif stra == 'FlatB':
		template=bkg.Clone()
	elif stra == 'FlatSB':
		template=sig.Clone()
		template.Scale(bkg.Integral()/sig.Integral())
		template.Add(bkg)
	else:
		print "No strategy "+stra+"!!!"
		return -1
	tot=template.Integral()
	step=tot/float(nbins)
	new_bins=[var_min]
	y=0
	for i in range(template.GetNbinsX()):
		y+=template.GetBinContent(i+1)
		if new_bins[len(new_bins)-1]> 0.25 and old_bins[i]-new_bins[len(new_bins)-1] < 0.009:
			continue
		if y > step:
			if abs(y-step) > abs(y-template.GetBinContent(i+1)-step) and old_bins[i] > new_bins[len(new_bins)-1]+0.0000001 and (old_bins[i]-new_bins[len(new_bins)-1] >= 0.009 or new_bins[len(new_bins)-1] < 0.25):
				new_bins.append(old_bins[i])
				y=template.GetBinContent(i+1)
				#y=0
			else:
				new_bins.append(old_bins[i+1])
				y=0
	if new_bins[len(new_bins)-1]<var_max-0.000001:
		new_bins.append(var_max)
	new_bins=np.array(new_bins)
	return new_bins
	
def QuantileLimitProtection(new_bins,stra,nbins):
	n=int(len(new_bins)-1)
	nbins=int(nbins)
	if n >= nbins:
		return True,new_bins
	elif nbins-n == 1:
		new_bins=np.insert(new_bins,1,ShiftToOldBins((new_bins[0]+new_bins[1])/2))
	elif nbins-n == 2:
		new_bins=np.insert(new_bins,1,ShiftToOldBins(new_bins[0]+2.*(new_bins[1]-new_bins[0])/3))
		new_bins=np.insert(new_bins,1,ShiftToOldBins(new_bins[0]+1.*(new_bins[1]-new_bins[0])/3))
	elif nbins-n == 3:
		new_bins=np.insert(new_bins,1,ShiftToOldBins(new_bins[0]+3.*(new_bins[1]-new_bins[0])/4))
		new_bins=np.insert(new_bins,1,ShiftToOldBins(new_bins[0]+2.*(new_bins[1]-new_bins[0])/4))
		new_bins=np.insert(new_bins,1,ShiftToOldBins(new_bins[0]+1.*(new_bins[1]-new_bins[0])/4))
	else:
		print "Quantile limit protection overload!!!",stra,nbins,ch,ca,nbins-n
		return False,-992
	return True,new_bins
	
def EmptyBinProtection(new_bins,stra,nbins,re):
	h_DY=GetHist("DY",re).Rebin(len(new_bins)-1,"",new_bins)
	h_TT=GetHist("TT",re).Rebin(len(new_bins)-1,"",new_bins)
	h_bkg=GetTotalBkg(re).Rebin(len(new_bins)-1,"",new_bins)
	for i in range(len(new_bins)-1):
		if h_bkg.GetBinContent(i+1) < 1.2*h_bkg.GetBinError(i+1):
		#if h_bkg.GetBinContent(i+1) < 0.1:
		if h_bkg.GetBinContent(i+1) < 0.18 or h_DY.GetBinContent(i+1) < 1e-7 or h_DY.GetBinContent(i+1) < 1e-7:
		#print h_DY.GetBinContent(i+1)*h_DY.GetEntries()/h_DY.Integral()
		#if h_DY.GetBinContent(i+1)*h_DY.GetEntries()/h_DY.Integral() < 1 or h_TT.GetBinContent(i+1)*h_TT.GetEntries()/h_TT.Integral() < 1 or h_bkg.GetBinContent(i+1) < 0.15:
			print "Fail empty bin protection:",stra,nbins,ch,ca,"(",h_bkg.GetBinLowEdge(i+1),",",h_bkg.GetBinLowEdge(i+2),")"
			return False
	return True

def TooScatteredProtection(new_bins,stra,nbins):
	for i in range(len(new_bins)-1):
		if new_bins[i]>0.25 and new_bins[i+1]-new_bins[i]<0.009:
			print "Fail too scattered distribution protection:",stra,nbins,ch,ca
			return False
	return True

parser = argparse.ArgumentParser(description='Command line parser of hyperparameters scan options')
parser.add_argument('--in', dest='fin', help='Input histo files')
parser.add_argument('--ch', dest='ch', help='channels to use (str of channels names separated by space)', default=None)
parser.add_argument('--year', dest='year', help='year to use', default=None)
parser.add_argument('--var', dest='var', help='variable to use', default=None)
parser.add_argument('--ca', dest='ca', help='categories to use (str of selection complete names separated by space)')
#parser.add_argument('--dir', dest='dir', help='directory where to find the distributions with minibins')
parser.add_argument('--analysis', dest='analysis', help='type of analysis we are doing: GGF or VBF')
#parser.add_argument('--cfg', dest='config', help='config files, used to read bkgs')
parser.add_argument('--bchmk', dest='bchmk',help='name of benchmark chosen for optimization')
args = parser.parse_args()

global fin
global ch
global year
global var
global ca
global analysis
global benchmark

finname=args.fin
ch=args.ch
year=args.year
var=args.var
ca=args.ca
analysis=args.analysis
#configname=args.config
#config=ConfigReader(configname)
benchmark=args.bchmk

fin=rt.TFile(finname,"read")

global backgrounds
backgrounds=['DY',"TT","Others","QCD"]

global strategies
strategies=['ConstSize','FlatS','FlatB','FlatSB']
#strategies=['FlatSB']

global regions
regions=['SR','SStight','OSinviso','SSinviso']
re='SR'

global nbins_test
nbins_test=range(1,55)
#nbins_test=[9]

template=GetHist(benchmark,re)
nPoints=template.GetNbinsX()
global old_bins
old_bins=[template.GetBinLowEdge(1)]
for ibin in range(1,nPoints+1):
	old_bins.append(template.GetBinLowEdge(ibin+1))
old_bins=np.array(old_bins)


global var_min
global var_max
var_min=round(old_bins[0],0)
var_max=round(old_bins[nPoints],0)
print var_min,var_max

parent='/home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_Jona/'
if not os.path.exists(parent+ch+"_"+ca):
	os.system("mkdir "+parent+ch+"_"+ca)


h_bkg=GetTotalBkg(re)
h_sig=GetHist(benchmark,re)
for stra in strategies:
	for nbins in nbins_test:
		new_bins=NewBin(h_sig,h_bkg,stra,nbins)
		#prot1,new_bins=QuantileLimitProtection(new_bins,stra,nbins)
		#print prot1,new_bins
		#if prot1:
		#	prot2=EmptyBinProtection(new_bins,stra,nbins,'SR')
			#prot3=TooScatteredProtection(new_bins,stra,nbins)
		#if prot1 and prot2 and prot3:
		#if prot1 and prot2:
		nbins_real=len(new_bins)-1
		if EmptyBinProtection(new_bins,stra,nbins_real,'SR'):
			path=ch+"_"+ca+"/"+stra+"_"+str(nbins_real)
			#print stra,nbins,len(new_bins)-1,new_bins
			if os.path.exists(parent+path):
				continue
			
			os.system("mkdir "+parent+path)
			
			finf=open(parent+path+"/Binning.txt","w")
			finf.truncate()
			for i,bins in enumerate(new_bins):
				finf.write(str(round(bins,4)))
				if i < len(new_bins)-1:
					finf.write(", ")
			finf.close()
			
			fout=rt.TFile(parent+path+"/analyzedOutplotter.root","recreate")
			h={}
			for pr in ['data_obs',benchmark] + backgrounds:
				for re in regions:
					if pr is 'QCD' and re is not 'SR':
						continue
					h[pr]=GetHist(pr,re).Rebin(len(new_bins)-1,"",new_bins)
			h['TotalBkg']=GetTotalBkg('SR').Rebin(len(new_bins)-1,"",new_bins)
			h['TotalBkg'].SetName("TotalBkg_"+ca+"_"+re+"_"+var)
			fout.Write()
			fout.Close()
			
			os.system('python '+parent+'write_res_card.py -f '+parent+path+'/analyzedOutplotter.root -o BinOpt_'+ch+ca+' -c '+ch+' -s '+ca+' -y '+year+' --signal '+benchmark+' -b 0 -u 0 -i /home/llr/cms/liugeliang/HH_bbtautau/CMSSW_11_1_0_pre6/src/KLUBAnalysis/config/mainCfg_'+ch+'_Legacy2018.cfg -r 1 -p '+parent)
			os.system('rm -r {1}{0}/*BinOpt*'.format(path,parent))
			os.system('mv {1}*BinOpt_{2}{3} {1}{0}/'.format(path,parent,ch,ca))
			os.system('combineCards.py -S {2}{0}/*BinOpt*/*/hhres*.{1}.txt > {2}{0}/comb.{1}.txt'.format(path,benchmark,parent))
			#os.system("echo 'SignalScale rateParam * {0} 0.01' > {1}add.txt".format(benchmark,parent))
			#os.system("cat {2}add.txt >> {2}{0}/comb.{1}.txt".format(path,benchmark,parent))
			os.system("echo 'SignalScale rateParam * {1} 0.01' >> {2}{0}/comb.{1}.txt".format(path,benchmark,parent))
			os.system("text2workspace.py {2}{0}/comb.{1}.txt -o {2}{0}/comb.{1}.root".format(path,benchmark,parent))
		
#limit=open(parent+ch+"_"+ca+"/get_limits.sh","w")
#limit.truncate()
for stra in strategies:
	for nbins in nbins_test:
		path=ch+"_"+ca+"/"+stra+"_"+str(nbins)
		if os.path.exists(parent+path):
			os.system("combine -M AsymptoticLimits /home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_Jona/{0}/comb.{1}.root --run blind --noFitAsimov -m {2} -n .{3}.{4} --freezeParameters SignalScale > /home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_Jona/{0}/comb.{1}.log".format(path,benchmark,benchmark.split('n')[1],ch,ca))
#limit.close()
#os.system('chmod +x '+parent+ch+"_"+ca+"/get_limits.sh")
#os.system('t3submit '+parent+ch+"_"+ca+"/get_limits.sh")

