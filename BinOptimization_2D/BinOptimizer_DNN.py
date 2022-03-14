import os
import time
import threading
import ROOT as rt
import numpy as np
import pandas as pd
from array import array
import argparse
import bayes_opt

dnn_min=0.
dnn_max=1.

min_step=0.0001
max_nsteps=10000
#max_nbins=25

tot_bkg_thr=1.
dy_thr=1e-7
tt_thr=1e-7
qcd_thr=1e-7
#binwidth_thr=0.2
min_width=0.0001

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

def NewBin(sig,bkg,stra,nbins):
	if stra == 'ConstSize':
		new_bins=np.array([round(x,4) for x in np.linspace(num=nbins+1,start=dnn_min,stop=dnn_max)])
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
	
	@staticmethod
	def fromPoint(point):
		temp=[0.]
		for n in range(len(point)):
			if temp[-1]+point[str(n)]<1:
				temp.append(temp[-1]+point[str(n)])
			else:
				break
		temp.append(1.)
		temp=np.array(temp)
		return Binning(temp)
	
	@staticmethod
	def fromEntry(path):
		if not os.path.exists(path):
                        print ("Folder not exists! ",path)
                        return -1
		f=open(path+"/Binning.txt","r")
		temp=np.array([float(i) for i in f.readline().split(", ")])
		binning=Binning(temp)
		#binning.limit=10*parseFile("{0}/comb.{1}.log".format(path,benchmark))
		return binning
	
	def Protection(self,binwidth_thr,tot_bkg_y,ref_bkg_y):
		for i in range(len(self.edges)):
			self.edges[i]=round(self.edges[i],4)
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
		
		#if smooth is False:
		#	return True
		flat=False
		while flat == False and self.nbins>3:
			flat=True
			for i in range(self.nbins):
				if i == 0:
					temp=(self.edges[i+1]-self.edges[i])>binwidth_thr*(self.edges[i+2]-self.edges[i+1]) or (self.edges[i+2]-self.edges[i+1])>0.5
				elif i == self.nbins-1:
					temp=(self.edges[i+1]-self.edges[i])>binwidth_thr*(self.edges[i]-self.edges[i-1]) or (self.edges[i]-self.edges[i-1])>0.5
				else:
					temp=((self.edges[i+1]-self.edges[i])>binwidth_thr*(self.edges[i+2]-self.edges[i+1]) or (self.edges[i+2]-self.edges[i+1])>0.5) and ((self.edges[i+1]-self.edges[i])>binwidth_thr*(self.edges[i]-self.edges[i-1]) or (self.edges[i]-self.edges[i-1])>0.5)
				if self.edges[i]<0.25:
					temp=True
				flat=flat and temp
				if temp == False:
					if i == 0:
						e=self.edges[i]
					elif i == self.nbins-1:
						e=self.edges[i+1]
					elif self.edges[i+2]-self.edges[i+1] > self.edges[i]-self.edges[i-1]:
						e=self.edges[i+1]
					elif self.edges[i+2]-self.edges[i+1] < self.edges[i]-self.edges[i-1]:
						e=self.edges[i]
					else:
						e=round((self.edges[i+1]+self.edges[i])/2,4)
					self.edges[i+1]=e
					self.edges[i]=e
			delete=[]
			for i in range(self.nbins):
				if self.edges[i]==self.edges[i+1]:
					delete.append(i)
			self.edges=np.delete(self.edges,delete)
			self.nbins=len(self.edges)-1
			
		
		return True
	
	def toPoint(self):
		point={}
		for i in range(min(self.nbins,max_nbins)):
			point[str(i)]=self.edges[i+1]-self.edges[i]
		for i in range(self.nbins,max_nbins):
			point[str(i)]=np.random.rand()*0.991+0.009
		return point

	def ComputeLimit(self,out):#,tot_bkg_y,ref_bkg_y):
		if not os.path.exists(out):
			os.system("mkdir "+out)
		#prot=self.Protection(tot_bkg_y,ref_bkg_y)
		#if not prot:
		#	return -1
		fbin=open(out+"/Binning.txt","w")
		fbin.truncate()
		fbin.write("Number of bins: "+str(self.nbins)+"\n")
		fbin.write(ArrayToStr(self.edges))
		fbin.close()
		
		fout=rt.TFile(out+"/analyzedOutplotter.root","recreate")
		h={}
		for pr in ['data_obs',benchmark] + backgrounds:
			for re in regions:
				if pr == 'QCD' and re != 'SR':
					continue
				h[pr]=GetHist(pr,re).Rebin(self.nbins,"",self.edges)
		h['totalBkg']=GetTotalBkg('SR').Rebin(self.nbins,'',self.edges)
		fout.Write()
		fout.Close()
		os.system('python '+parent+'write_res_card.py -f '+out+'/analyzedOutplotter.root -o BinOpt_'+ch+ca+' -c '+ch+' -s '+ca+' -y '+year+' --signal '+benchmark+' -b 0 -u 0 -i /home/llr/cms/liugeliang/HH_bbtautau/CMSSW_11_1_0_pre6/src/KLUBAnalysis/config/mainCfg_'+ch+'_Legacy2018.cfg -r 1 -p '+parent)
		os.system('rm -r {0}/*BinOpt*'.format(out))
		os.system('mv {1}*BinOpt_{2}{3} {0}/'.format(out,parent,ch,ca))
		os.system('combineCards.py -S {0}/*BinOpt*/*/hhres*.{1}.txt > {0}/comb.{1}.txt'.format(out,benchmark,parent))
		os.system("echo 'SignalScale rateParam * {1} 0.001' >> {0}/comb.{1}.txt".format(out,benchmark))
		os.system("text2workspace.py {0}/comb.{1}.txt -o {0}/comb.{1}.root".format(out,benchmark))
		os.system("combine -M AsymptoticLimits {0}/comb.{1}.root --run blind --noFitAsimov -m {2} -n .{3}.{4} --freezeParameters SignalScale > {0}/comb.{1}.log".format(out,benchmark,benchmark.split('n')[1],ch,ca))
		self.limit=parseFile("{0}/comb.{1}.log".format(out,benchmark))
		os.system("rm -r {0}/*BinOpt*".format(out))
                os.system("rm {0}/comb*root".format(out))

	def isEquivalent(self, other):
        	if len(self.edges) != len(other.edges):
            		return False
        	return np.count_nonzero(self.edges == other.edges) == len(self.edges)
    
	def isBetter(self, other):
		if other == None or other.limit is None:
			return True
		return self.limit*1e3+self.nbins < other.limit*1e3+other.nbins
        	#if self.limit != other.limit:
            	#	return self.limit < other.limit
        	#return len(self.edges) < len(other.edges)

class BayesianOptimization:
	def __init__(self,work_area,acq,kappa,xi):
		self.work_area=work_area
		self.work_point=0
		self.stable=0
		self.fin=fin
		self.tot_bkg_y=HistToArray(GetTotalBkg("SR"))
		self.ref_bkg_y={}
		self.ref_bkg_y['TT']=HistToArray(GetHist("TT","SR"))
		self.ref_bkg_y['DY']=HistToArray(GetHist("DY","SR"))
		self.ref_bkg_y['QCD']=HistToArray(GetHist("QCD","SR"))
		self.binnings=[]
		self.BestBinning=None
		
		self.bounds={}
		for i in range(max_nbins):
			self.bounds[str(i)]=(min_width,1.)
		
		self.optimizer=bayes_opt.BayesianOptimization(f=None, pbounds=self.bounds, random_state=12345, verbose=1)
		self.utilities = [
            		bayes_opt.util.UtilityFunction(kind='ucb', kappa=kappa, xi=xi),# kappa_decay=0.99),
            		bayes_opt.util.UtilityFunction(kind='ei', kappa=kappa, xi=xi),
            		bayes_opt.util.UtilityFunction(kind='poi', kappa=kappa, xi=xi),
        		]
		
		self.bounds_transformer=bayes_opt.SequentialDomainReductionTransformer()
                self.bounds_transformer.initialize(self.optimizer.space)	
		
		f=open(self.work_area+"/Load.txt","w")
		f.truncate()
		f.write("tot_bkg_thr:"+str(tot_bkg_thr))
		f.write("\ndy_thr:"+str(dy_thr))
		f.write("\ntt_thr:"+str(tt_thr))
		f.write("\nqcd_thr:"+str(qcd_thr)+'\n')
		f.close()
		
		self.RegisterOld()
		#self.RegisterRandom(10-len(self.binnings))
		self.stable=0
		self.prot=0
		#self.BestBinning=None
		
		
		for count in range(300):
		#while self.stable < 30 and self.prot < 30:
			f=open(self.work_area+"/Load.txt","a")
			
			nPoints=self.RegisterSuggest()
			if nPoints == -1:
				f.write("Protection not achieved!\n")
			elif nPoints == -2:
				f.write("Already found!\n")
			else:
				f.write("Number of bins: "+str(self.binnings[-1].nbins))
				f.write("	Target: "+str(self.optimizer.space.target[-1])+"\n")
			
			#prot=binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
			#temp_binning=self.binnings[-1]
			#prot=temp_binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
			#if prot:
			#	find,ib=self.Find(temp_binning)
			#	find=False
			#	if not find:
			#		temp_binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
			#		self.work_point+=1
			#		self.Register(temp_binning)
			#		f.write("Number of bins: "+str(self.binnings[-1].nbins))
			#		f.write("       Target: "+str(self.optimizer.space.target[-1])+"\n")
			#	else:
			#		f.write("Already found!"+str(ib)+"\n")
			#else:
			#	f.write("Protection not achieved!\n")
			
			nPoints=self.binnings[-1].nbins-1
			if count%2 == 1:
				self.UpdateBounds(f,nPoints)	
		
		#	nPoints=self.RegisterSmoothedBest()
                #        if nPoints == -1:
                #                f.write("Protection not achieved!\n")
                #        elif nPoints == -2:
                #                f.write("Already found!\n")
                #        else:
                #                f.write("Number of bins: "+str(self.binnings[-1].nbins))
                #                f.write("       Target: "+str(self.optimizer.space.target[-1])+"\n\n")
		
			nPoints=self.RegisterMergedBest()
			if nPoints == -1:
                                f.write("Protection not achieved!\n")
                        elif nPoints == -2:
                                f.write("Already found!\n")
                        else:
				f.write("Number of bins: "+str(self.binnings[-1].nbins))
                        	f.write("       Target: "+str(self.optimizer.space.target[-1])+"\n\n")
			
	#		if count%2 == 0:
	#			self.UpdateBounds(f)
			f.close()
		
	def Find(self,binning):
		for ib,b in enumerate(self.binnings):
			if binning.isEquivalent(b):
				return True,ib
		return False,-1
	
	def RegisterOld(self):
		for stra in ['ConstSize','FlatB','FlatS','FlatSB']:
			for nbins in range(1,max_nbins-20):
				if nbins%5 != 0 and nbins>10:
					continue
				#if os.path.exists("{0}/{1}_{2}".format(old_input,stra,str(nbins))):
				f=open(self.work_area+"/Load.txt","a")
				f.write(stra+"_"+str(nbins)+"\n")
				#binning=Binning.fromEntry("{0}/{1}_{2}".format(old_input,stra,str(nbins)))
				new_bins=NewBin(GetHist(benchmark,'SR'),GetTotalBkg('SR'),stra,nbins)
				binning=Binning(new_bins)
				print binning.limit,binning.nbins
				prot=binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
				if not prot:
					f.write('Protection cannot be achieved!\n\n')
					f.close()
					continue
				find,ib=self.Find(binning)
				if find:
					f.write('Already found!\n\n')
					f.close()
					continue
				binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
				self.work_point+=1
				self.Register(binning)
				f.write("Number of bins: "+str(self.binnings[-1].nbins))
				f.write("       Target: "+str(self.optimizer.space.target[-1])+"\n\n")
				f.close()
	#			self.UpdateBounds(binning.nbins-1)

	def RegisterRandom(self,n):
		for i in range(n):
			temp=np.sort(np.random.rand(self.nbins-1))
			temp=np.concatenate((np.array([0.]),np.array([round(j,4) for j in temp]),np.array([1.])))
			binning=Binning(temp)
			prot=binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
			if not prot:
				continue
			find,ib=self.Find(binning)
			if find:
				continue
			binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
			self.work_point+=1
			self.Register(binning)
		#self.stable=0

	def RegisterSuggest(self):
		point=self.Suggest(0)
		#os.system("echo {0} > {1}/Point_{2}/Point.txt".format(str(point),self.work_area,str(self.work_point)))
		binning=Binning.fromPoint(point)
		#f.write(ArrayToStr(binning.edges)+'\n')
		prot=binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
		#f.write(ArrayToStr(binning.edges)+'\n')
		if not prot:
			self.prot+=1
			#self.RegisterRandom(1)
			return -1,binning
		find,ib=self.Find(binning)
		if find:
			return -2,binning
		else:
			#self.stable+=1
			#return -3
			binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
			self.work_point+=1
			self.Register(binning)
		
		#merged_binning=Binning(np.delete(binning.edges,[1]))
		#find,ib=self.Find(merged_binning)
		#if not find:
		#	merged_binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
                #        self.work_point+=1
                #        self.Register(merged_binning)
		
		return binning.nbins

	def RegisterSmoothedBest(self):
		target_space=self.optimizer.space
                opt_binning=Binning.fromPoint(target_space.array_to_params(target_space.params[np.argmax(target_space.target)]))
		#merged_binning=Binning(np.delete(opt_binning.edges,[1]))
		prot=opt_binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
		if not prot:
			return -1
		find,ib=self.Find(opt_binning)
		if find:
			return -2
		else:
			opt_binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
                        self.work_point+=1
                        self.Register(opt_binning)
		return opt_binning.nbins
	
	def RegisterMergedBest(self):
                target_space=self.optimizer.space
                opt_binning=Binning.fromPoint(target_space.array_to_params(target_space.params[np.argmax(target_space.target)]))
                merged_binning=Binning(np.delete(opt_binning.edges,[1]))
                prot=merged_binning.Protection(0.2,self.tot_bkg_y,self.ref_bkg_y)
                if not prot:
                        return -1
                find,ib=self.Find(merged_binning)
                if find:
                        return -2
                else:
                        merged_binning.ComputeLimit('{0}/Point_{1}'.format(self.work_area,str(self.work_point)))
                        self.work_point+=1
                        self.Register(merged_binning)
                return merged_binning.nbins
	
	def Register(self,binning):
		if binning.limit == None:
			print ("No limit, cannot register!")
			return -2
		loss=-(binning.limit*1e3+binning.nbins)
		point=binning.toPoint()
		self.optimizer.register(params=point,target=loss)
		self.binnings.append(binning)
		if self.BestBinning == None:
			self.BestBinning=binning
		else:
			if self.BestBinning.isBetter(binning):
				self.stable+=1
			else:
				self.stable=0
				self.BestBinning=binning
		#self.UpdateBounds()
	
	def UpdateBounds(self,f,nPoints):
		target_space=self.optimizer.space
		opt_binning=Binning.fromPoint(target_space.array_to_params(target_space.params[np.argmax(target_space.target)]))
		#nPoints=opt_binning.nbins-1
		new_bounds=self.bounds_transformer.transform(self.optimizer.space)
		prev=True
		for i in range(max_nbins):
			if prev:
				prev=i < nPoints
				#prev=new_bounds[str(i)][1]-new_bounds[str(i)][0]<0.1 and i < nPoints
			else:
				new_bounds[str(i)]=np.array([min_width,1.])
		for i in range(max_nbins):
			if new_bounds[str(i)][0] > new_bounds[str(i)][1]:
				new_bounds[str(i)]=np.array([new_bounds[str(i)][1],new_bounds[str(i)][0]])
			#new_bounds[str(i)][0]*=0.5
			#new_bounds[str(i)][1]*=1.5
			if new_bounds[str(i)][0]<min_width:
				new_bounds[str(i)][0]=min_width
			if new_bounds[str(i)][1]>1:
				new_bounds[str(i)][1]=1.
			#if new_bounds[str(i)][1]-new_bounds[str(i)][0]<0.001:
			#	mid=(new_bounds[str(i)][1]+new_bounds[str(i)][0])/2
			#	new_bounds[str(i)][1]=mid+0.0005
			#	new_bounds[str(i)][0]=mid-0.0005
			#if new_bounds[str(i)][1]-new_bounds[str(i)][0]<0.1:
			#	new_bounds[str(i)][1]=new_bounds[str(i)][0]+0.1
		#print new_bounds
		self.optimizer.set_bounds(new_bounds)
		f.write("nbins of best binning: "+str(opt_binning.nbins)+"\n")
		#for i in range(max_nbins):
		#	f.write("{0}: [{1},{2}]; ".format(str(i),str(round(old_bounds[i][0],3)),str(round(old_bounds[i][1],3))))
                #f.write("\n")
		for i in range(max_nbins):
			f.write("{0}: [{1},{2}]; ".format(str(i),str(round(new_bounds[str(i)][0],4)),str(round(new_bounds[str(i)][1],4))))
		f.write("\n\n")
	
	def Suggest(self,util_idx):
		point=self.optimizer.suggest(self.utilities[util_idx])
		return point

	def Write(self):
		f=open(self.work_area+"/BestBin.txt","w")
		f.truncate()
		f.write("Number of bins: "+str(self.BestBinning.nbins))
		f.write("\nBinning scheme:")
		f.write("\n"+ArrayToStr(self.BestBinning.edges))
		f.write("\nExpected limit: "+str(self.BestBinning.limit))


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
parser.add_argument('--max_nbins', dest='max_nbins',help='maximum number of bins allowed')
args = parser.parse_args()

global fin
global ch
global year
global var
global ca
global analysis
global benchmark
global max_nbins

finname=args.fin
ch=args.ch
year=args.year
var=args.var
ca=args.ca
analysis=args.analysis
benchmark=args.bchmk
max_nbins=int(args.max_nbins)

fin=rt.TFile(finname,"read")

global backgrounds
backgrounds=['DY',"TT","Others","QCD"]

global regions
regions=['SR','SStight','OSinviso','SSinviso']

parent='/home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_2D/'
path=parent+ch+"_"+ca
if not os.path.exists(path):
        os.system("mkdir "+path)

#old_input='/home/llr/cms/liugeliang/HH_bbtautau/KLUBAnalysis/BinOptimization_Jona/OnlyBinnings/{0}_{1}'.format(ch,ca)
#old_input='/data_CMS/cms/liugeliang/HHbbtautau_BinOpt/Jona/backup2/{0}_{1}'.format(ch,ca)

template=GetHist(benchmark,'SR')
template.SetName('template')
template.SetTitle('template')
nPoints=template.GetNbinsX()
global old_bins
old_bins=[template.GetBinLowEdge(1)]
for ibin in range(1,nPoints+1):
        old_bins.append(template.GetBinLowEdge(ibin+1))
old_bins=np.array(old_bins)


BestBinning=None

h={}
h['TotalBkg']=GetTotalBkg('SR').Rebin(1,'',np.array([0.,1.]))
for bkg in backgrounds:
	h[bkg]=GetHist(bkg,'SR').Rebin(1,'',np.array([0.,1.]))
h['noTT']=h['DY'].Clone()
h['noTT'].Add(h['QCD'])
h['noTT'].Add(h['Others'])
h['noQCD']=h['DY'].Clone()
h['noQCD'].Add(h['TT'])
h['noQCD'].Add(h['Others'])
tot_bkg_e=h['TotalBkg'].GetBinError(1)
tot_bkg=h['TotalBkg'].GetBinContent(1)
nott_e=h['noTT'].GetBinError(1)
nott=h['noTT'].GetBinContent(1)
noqcd_e=h['noQCD'].GetBinError(1)
noqcd=h['noQCD'].GetBinContent(1)
if nott_e > 0.9*tot_bkg_e and nott > 0.8*tot_bkg:
	tt_thr=-1
if noqcd_e > 0.9*tot_bkg_e and noqcd > 0.8*tot_bkg:
	qcd_thr=-1

os.system("mkdir "+path)
ob=BayesianOptimization(path,'ucb',2.57,0.1)
ob.Write()
		
'''
fout=open(path+"/BestBin.txt","w")
fout.truncate()
fout.write("Number of bins: "+str(BestBinning.nbins))
fout.write("\nBinning scheme:\n")
fout.write(ArrayToStr(BestBinning.edges))
fout.write("\nExpected limit: "+str(BestBinning.limit))
fout.close()
'''
