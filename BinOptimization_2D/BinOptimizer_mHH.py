import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array
import argparse

dnn_min=0.
dnn_max=1.
mHH_min=250.
mHH_max=3500.

Tot_bkg_thr=5.
tot_bkg_thr=1.
dy_thr=1e-7
#tt_thr=1e-7
#qcd_thr=1e-7

min_step=0.001
max_nsteps=1000
min_width=0.001

def ArrayToStr(a):
    return ', '.join([ str(x) for x in a ])

def GetHist(pr,re):
        fullname=pr+"_"+ca+"_"+re+"_"+var
        if not fin.GetListOfKeys().Contains(fullname):
                print ("*** WARNING: histo " , fullname , " not available")
                return -1
        else:
                return fin.Get(fullname).Clone()

def GetTotalBkg(re):
        for i,bkg in enumerate(backgrounds):
                fullname=bkg+"_"+ca+"_"+re+"_"+var
                if not fin.GetListOfKeys().Contains(fullname):
                        print ("*** WARNING: histo " , fullname , " not available")
                        return -1
                if i == 0:
                        h=fin.Get(fullname).Clone()
                else:
                        h.Add(fin.Get(fullname))
        h.SetName('totalBkg_'+ca+"_"+re+"_"+var)
        h.SetTitle('totalBkg_'+ca+"_"+re+"_"+var)
        return h

def HistToArray(h):
        return np.array([h.GetBinContent(i+1) for i in range(h.GetNbinsX())])

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
        print ("did not find any expected in file: " , filename, 'CL=', CL, 'exp?=', exp)
        return 1e10
    else:
        return matches[-1]

def TwoSides(h,mass):
	tot=h.Integral()
	N=h.GetNbinsX()
	if mass < 600:
		r_thr=0.1
	else:
		r_thr=0.05
	l_thr=0.2
	temp=0
	for i in range(1,N+1):
		temp+=h.GetBinContent(i)
		if temp>tot*l_thr:
			break
	leftside=i+1
	temp=0
	for i in reversed(range(1,N+1)):
		temp+=h.GetBinContent(i)
                if temp>tot*r_thr:
                        break
	rightside=i
	return leftside,rightside

def NewBin(sig,bkg,stra,nbins):
	if stra == 'ConstSize':
		new_bins=np.array([round(x,3) for x in np.linspace(num=nbins+1,start=dnn_min,stop=dnn_max)])
		#make sure each edge in new binning is in old binning, otherwise modify it to be the closest edge in old binning
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
	new_bins=[dnn_min]
	y=0
	for i in range(template.GetNbinsX()):
		y+=template.GetBinContent(i+1)
		#if new_bins[len(new_bins)-1]> 0.25 and old_bins[i]-new_bins[len(new_bins)-1] < 0.009:
		#	continue
		if y > step:
			if abs(y-step) > abs(y-template.GetBinContent(i+1)-step) and old_bins[i] > new_bins[len(new_bins)-1]+0.0000001:
				new_bins.append(old_bins[i])
				y=template.GetBinContent(i+1)
				#y=0
			else:
				new_bins.append(old_bins[i+1])
				y=0
	if new_bins[len(new_bins)-1]<dnn_max-0.000001:
		new_bins.append(dnn_max)
	new_bins=np.array(new_bins)
	return new_bins

class Binning:
	def __init__(self,edges):
		self.edges=edges
		self.nbins=len(edges)-1
		self.limit=None
		self.target=None
	
	def Protection(self,tot_bkg_y,ref_bkg_y,tt_thr,qcd_thr):
		for i in range(len(self.edges)):
			self.edges[i]=round(self.edges[i],3)
		self.edges=np.sort(self.edges)
		if self.edges[0] != 0:
			self.edges[0]=0
		if self.edges[-1] != 1:
			self.edges[-1]=1
		
		for i in reversed(range(self.nbins)):
			if self.edges[i] > self.edges[i+1]-min_width:
				self.edges[i]=self.edges[i+1]-min_width
			if self.edges[i] < 0:
				print ("Protection cannot be achieved!")
				return False
			while np.sum(tot_bkg_y[int(round(max_nsteps*self.edges[i],0)):int(round(max_nsteps*self.edges[i+1],0))]) < tot_bkg_thr or np.sum(ref_bkg_y['TT'][int(round(max_nsteps*self.edges[i],0)):int(round(max_nsteps*self.edges[i+1],0))]) < tt_thr or np.sum(ref_bkg_y['DY'][int(round(max_nsteps*self.edges[i],0)):int(round(max_nsteps*self.edges[i+1],0))]) < dy_thr or np.sum(ref_bkg_y['QCD'][int(round(max_nsteps*self.edges[i],0)):int(round(max_nsteps*self.edges[i+1],0))]) < qcd_thr:
				self.edges[i]=self.edges[i]-min_step
				if self.edges[i] < 0:
					print ("Protection cannot be achieved!")
                                	return False
			if self.edges[i] < 0:
				print ("Protection cannot be achieved!")
				return False
		return True

	def ComputeLimit(self,out,histos,benchmark):
		if not os.path.exists(out):
			os.system("mkdir "+out)
		fbin=open(out+"/Binning.txt","w")
		fbin.truncate()
		fbin.write("Number of bins: "+str(self.nbins)+"\n")
		fbin.write(ArrayToStr(self.edges))
		fbin.close()
		
		fout=rt.TFile(out+"/analyzedOutplotter.root","recreate")
		htemp={}
		for pr in ['data_obs','totalBkg',benchmark] + backgrounds:
			for re in regions:
				if pr == 'QCD' and re != 'SR':
					continue
				if pr == 'totalBkg' and re != 'SR':
                                        continue
				htemp[pr+re]=histos[pr+re].Rebin(self.nbins,"",self.edges)
				htemp[pr+re].SetName('{0}_{1}_{2}_{3}'.format(pr,ca,re,'DNNoutSM_kl_1'))
				htemp[pr+re].SetTitle('{0}_{1}_{2}_{3}'.format(pr,ca,re,'DNNoutSM_kl_1'))
				htemp[pr+re].Write()
		#fout.Write()
		#if out == "{0}/Bin1_Point1/ConstSize1".format(path):
		#	c=rt.TCanvas()
		#	htemp[bchmk+'SR'].Draw()
		#	c.SaveAs('test.png')
		fout.Close()
		os.system('python '+parent+'write_res_card.py -f '+out+'/analyzedOutplotter.root -o BinOpt_'+ch+ca+' -c '+ch+' -s '+ca+' -y '+year+' --signal '+benchmark+' -b 0 -u 0 -i /home/llr/cms/liugeliang/HH_bbtautau/CMSSW_11_1_0_pre6/src/KLUBAnalysis/config/mainCfg_'+ch+'_Legacy2018.cfg -r 1 -p '+parent)
		os.system('rm -r {0}/*BinOpt*'.format(out))
		os.system('mv {1}*BinOpt_{2}{3} {0}/'.format(out,parent,ch,ca))
		os.system('combineCards.py -S {0}/*BinOpt*/*/hhres*.{1}.txt > {0}/comb.{1}.txt'.format(out,benchmark,parent))
		os.system("echo 'SignalScale rateParam * {1} 0.001' >> {0}/comb.{1}.txt".format(out,benchmark))
		os.system("text2workspace.py {0}/comb.{1}.txt -o {0}/comb.{1}.root".format(out,benchmark))
		os.system("combine -M AsymptoticLimits {0}/comb.{1}.root --run blind --noFitAsimov -m {2} -n .{3}.{4} --freezeParameters SignalScale > {0}/comb.{1}.log".format(out,benchmark,benchmark.split('n')[1],ch,ca))
		self.limit=parseFile("{0}/comb.{1}.log".format(out,benchmark))
		self.target=-(self.limit*1e3+self.nbins)
                os.system("rm -r {0}/*BinOpt*".format(out))
		os.system("rm {0}/comb*root".format(out))
	
	def isEquivalent(self, other):
        	if len(self.edges) != len(other.edges):
            		return False
        	return np.count_nonzero(self.edges == other.edges) == len(self.edges)
    
	def isBetter(self, other):
		if other == None or other.limit is None:
			return True
		return self.target > other.target


class GetBestBinning:
	def __init__(self,work_area,Histos,low_bin,high_bin,bchmk):
		self.work_area=work_area
		self.load=self.work_area+'Load.txt'
		self.bchmk=bchmk
		self.histos={}
		for pr in ['data_obs','totalBkg',bchmk]+backgrounds:
                        for re in regions:
                                if pr == 'QCD' and re != 'SR':
                                        continue
				if pr == 'totalBkg' and re != 'SR':
                                        continue
                                self.histos[pr+re]=Histos[pr+re].ProjectionX("",low_bin,high_bin)
				self.histos[pr+re].SetName('{0}_{1}_{2}_{3}'.format(pr,ca,re,'DNNoutSM_kl_1'))
				self.histos[pr+re].SetTitle('{0}_{1}_{2}_{3}'.format(pr,ca,re,'DNNoutSM_kl_1'))

		
		self.tot_bkg_y=HistToArray(self.histos['totalBkgSR'])
		self.ref_bkg_y={}
		for pr in ['DY','TT','QCD']:
			self.ref_bkg_y[pr]=HistToArray(self.histos[pr+'SR'])
		self.binnings=[]
		self.BestBinning=None
		
		self.histos['noTT']=self.histos['DYSR'].Clone()
		self.histos['noTT'].Add(self.histos['QCDSR'])
		self.histos['noTT'].Add(self.histos['OthersSR'])
		self.histos['noQCD']=self.histos['DYSR'].Clone()
		self.histos['noQCD'].Add(self.histos['TTSR'])
		self.histos['noQCD'].Add(self.histos['OthersSR'])
		tot_e=self.histos['totalBkgSR'].Rebin(1,'',np.array([0.,1.])).GetBinError(1)
		tot_y=self.histos['totalBkgSR'].Rebin(1,'',np.array([0.,1.])).GetBinContent(1)
		nott_e=self.histos['noTT'].Rebin(1,'',np.array([0.,1.])).GetBinError(1)
		nott_y=self.histos['noTT'].Rebin(1,'',np.array([0.,1.])).GetBinContent(1)
		noqcd_e=self.histos['noQCD'].Rebin(1,'',np.array([0.,1.])).GetBinError(1)
		noqcd_y=self.histos['noQCD'].Rebin(1,'',np.array([0.,1.])).GetBinContent(1)
		if nott_e > 0.9*tot_e and nott_y > 0.8*tot_y:
			self.tt_thr=-1
		else:
			self.tt_thr=1e-7
		if noqcd_e > 0.9*tot_e and noqcd_y > 0.8*tot_y:
			self.qcd_thr=-1
		else:
			self.qcd_thr=1e-7
		
		if not os.path.exists(self.work_area):
			os.system("mkdir "+self.work_area)
		f=open(self.load,'w')
		f.truncate()
		f.write('tot_bkg_thr='+str(tot_bkg_thr)+'\n')
		f.write('dy_thr='+str(dy_thr)+'\n')
		f.write('tt_thr='+str(self.tt_thr)+'\n')
		f.write('qcd_thr='+str(self.qcd_thr)+'\n\n')
		f.close()

	def Find(self,binning):
                for ib,b in enumerate(self.binnings):
                        if binning.isEquivalent(b):
                                return True,ib
                return False,-1

	def Scan(self,Bin):
		for stra in ['ConstSize','FlatB','FlatS','FlatSB']:
			for nbins in range(1,max_nbins+1):
				if Bin == 1:
					r=1
				elif Bin == 2:
					r=2
				else:
					r=3
				if 'boost' not in ca:
					r=3
				if nbins%r != 0:
					continue
				new_bins=NewBin(self.histos[self.bchmk+'SR'],self.histos['totalBkgSR'],stra,nbins)
				binning=Binning(new_bins)
				f=open(self.load,'a')
				f.write('{0}_{1}:\n'.format(stra,str(nbins)))
				prot=binning.Protection(self.tot_bkg_y,self.ref_bkg_y,self.tt_thr,self.qcd_thr)
				if not prot:
					f.write('Protection cannot be achieved!\n\n')
					f.close()
					continue
				find,ib=self.Find(binning)
				if find:
					f.write('Already found!\n\n')
                                        f.close()
					continue
				binning.ComputeLimit(self.work_area+stra+str(nbins),self.histos,self.bchmk)
				self.binnings.append(binning)
				f.write('Limit: {0}\nNumber of bins: {1}\n\n'.format(binning.limit,binning.nbins))
				f.close()
				if self.BestBinning == None:
					self.BestBinning=binning
				else:
					if not self.BestBinning.isBetter(binning):
						self.BestBinning=binning


parser = argparse.ArgumentParser(description='Command line parser of hyperparameters scan options')
parser.add_argument('--in', dest='fin', help='Input histo files')
parser.add_argument('--ch', dest='ch', help='channels to use (str of channels names separated by space)', default=None)
parser.add_argument('--year', dest='year', help='year to use', default=None)
parser.add_argument('--ca', dest='ca', help='categories to use (str of selection complete names separated by space)')
parser.add_argument('--analysis', dest='analysis', help='type of analysis we are doing: GGF or VBF')
parser.add_argument('--max_nbins', dest='max_nbins',help='maximum number of bins allowed')
args = parser.parse_args()

finname=args.fin
fin=rt.TFile(finname,"read")
ch=args.ch
year=args.year
ca=args.ca
analysis=args.analysis
max_nbins=int(args.max_nbins)
var='DNNoutSM_kl_1_HHKin_mass_raw'

backgrounds=['DY',"TT","Others","QCD"]
regions=['SR','SStight','OSinviso','SSinviso']

parent='/home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_2D/'
path=parent+ch+"_"+ca
if not os.path.exists(path):
        os.system("mkdir "+path)

massPoints=[250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
Radions=['Radion'+str(mass) for mass in massPoints]

Histos={}
for mass in massPoints:
	for re in regions:
		Histos['Radion{0}'.format(mass)+re]=GetHist('Radion{0}'.format(mass),re)
for pr in ['data_obs']+backgrounds:
	for re in regions:
		if pr == 'QCD' and re is not 'SR':
			continue
		Histos[pr+re]=GetHist(pr,re)
Histos['totalBkg'+'SR']=GetTotalBkg('SR')

old_mHH_bins=[mHH_min]
for ibin in range(1,Histos['totalBkgSR'].GetNbinsY()+1):
	old_mHH_bins.append(Histos['totalBkgSR'].ProjectionY().GetBinLowEdge(ibin+1))
old_mHH_bins=np.array(old_mHH_bins)

old_bins=[dnn_min]
for ibin in range(1,Histos['totalBkgSR'].GetNbinsX()+1):
	old_bins.append(Histos['totalBkgSR'].ProjectionX().GetBinLowEdge(ibin+1))
old_bins=np.array(old_bins)

f=open('{0}/Load.txt'.format(path),'w')
f.truncate()
f.write('old_mHH_bins:\n[{0}]\n\n'.format(ArrayToStr(old_mHH_bins)))
f.close()

ibchmk=19
bchmk=Radions[ibchmk]
bchmks=[bchmk]
high_bin=Histos['totalBkgSR'].GetNbinsY()
new_mHH_bins=[old_mHH_bins[high_bin]]
low_bin,ex=TwoSides(Histos[bchmk+'SR'].ProjectionY(),massPoints[ibchmk])
Targets=[]
check=0
Bin=1
Point=1

while low_bin>0:
	if Histos['totalBkgSR'].ProjectionX("",low_bin,high_bin).Integral()<5:
                low_bin-=1
                continue
	f=open('{0}/Load.txt'.format(path),'a')
	f.write('Bin{0}_Point{1}:\n'.format(str(Bin),str(Point)))
	f.write('Range: [{0},{1}]\n'.format(str(old_mHH_bins[low_bin-1]),str(old_mHH_bins[high_bin])))
	f.write('Benchmark: {0}\n'.format(bchmk))
	opt=GetBestBinning("{0}/Bin{1}_Point{2}/".format(path,str(Bin),str(Point)),Histos,low_bin,high_bin,bchmk)
	opt.Scan(Bin)
	if opt.BestBinning == None:
		f.write('Protection cannot be achieved!\n\n')
		f.close()
		low_bin-=1
		Point+=1
		continue
	else:
		f.write('Best limit: {0}\n'.format(str(opt.BestBinning.limit)))
		f.write('Number of bins: {0}\n'.format(str(opt.BestBinning.nbins)))
	if len(Targets)>0:
		if opt.BestBinning.target<max(Targets):
			check+=1
		else:
			check=0
			best_low_bin=low_bin
	else:
		best_low_bin=low_bin
	f.write('Check point: {0}\n\n'.format(str(check)))
	f.close()
	Targets.append(opt.BestBinning.target)
	if check>=2:
		high_bin=best_low_bin-1
		new_mHH_bins.append(old_mHH_bins[high_bin])
		left,right=TwoSides(Histos[bchmk+'SR'].ProjectionY(),massPoints[ibchmk])
		while right>high_bin:
			ibchmk-=1
			if ibchmk < 0:
				new_mHH_bins[-1]=mHH_min
				break
			bchmk=Radions[ibchmk]
			#bchmks.append(bchmk)
			left,right=TwoSides(Histos[bchmk+'SR'].ProjectionY(),massPoints[ibchmk])
		if ibchmk < 0:
			break
		bchmks.append(bchmk)
		low_bin,ex=TwoSides(Histos[bchmk+'SR'].ProjectionY(),massPoints[ibchmk])
		Targets=[]
		check=0
		Bin+=1
		Point=1
	else:
		low_bin-=1
		Point+=1

if new_mHH_bins[-1]>mHH_min:
	new_mHH_bins.append(mHH_min)

new_mHH_bins=np.array(list(reversed(new_mHH_bins)))
bchmks=np.array(list(reversed(bchmks)))

f=open('{0}/BestBin.txt'.format(path),'w')
f.write('{0}\n'.format(ArrayToStr(new_mHH_bins)))
f.write('{0}\n'.format(ArrayToStr(bchmks)))
f.close()
