import ROOT
import collections
import fnmatch
import array

def makeHistoName(sample, sel, var, syst=""):
    name = sample +  "_" + sel + "_" + var
    if syst:
        name += "_"
        name += syst
    return name;

def makeHisto2DName(sample, sel, var1, var2, syst=""):
    name = sample +  "_" + sel + "_" + var1 + "_" + var2
    if syst:
        name += "_"
        name += syst
    return name;

def matchInDictionary(diction, pattern):
    matches = []
    for key in diction:
        if fnmatch.fnmatch(key, pattern):
            matches.append(key)
    return matches

def checkBinningCompatibility (newbinning, oldbinning):
    """ oldbinning must include newbinning boundaries """
    for x in newbinning:
        if not x in oldbinning: return False
    return True

class OutputManager:
    """ Handles the input from AnalysisHelper and manages the output to a file
    to be used for datacards and analysis"""
    def __init__(self):
        self.histos      = collections.OrderedDict()
        self.histos2D    = collections.OrderedDict()
        
        self.selections  = []
        self.sel_def     = []
        self.sel_regions = []
        self.variables   = []
        self.variables2D = []
        # self.samples     = []
        self.data        = []
        self.bkgs        = []
        self.sigs        = []

    def readAll (self, rootfile):
        """ read all histograms from rootfile """

        for key in rootfile.GetListOfKeys():
            kname = key.GetName()
            # print key, kname
            obj = rootfile.Get(kname)

            ## 1D plots
            if isinstance(obj, ROOT.TH1):
                self.histos[kname] = obj

            ## 2D plots
            if isinstance(obj, ROOT.TH2):
                self.histos2D[kname] = obj

        # if check:
        #     print '.... checking reading'
        #     for sample in self.samples:
        #         for sel in self.selections:
        #             for var in self.variables:
        #                 pass

    def groupTogether (self, newName, partsList):
        """ regroup the samples in partsList into a new sample that replaces old ones"""
        
        print "... replacing" , newName , " <== " , partsList

        if len(partsList) == 0:
            print "** ERROR: while grouping together histos, input list has size 0"
            raise ValueError

        ## 1D
        for sel in self.selections:
            for var in self.variables:
                newHistoName = makeHistoName(newName, sel, var)
                
                if newHistoName in self.histos:
                    print "** ERROR: while grouping together histos, name", newName, "already exists in input as", newHistoName, "cannot continue..."
                    raise ValueError

                for idx, oldName in enumerate(partsList):
                    oldHistoName = makeHistoName(oldName, sel, var)
                    if idx == 0:
                        self.histos[newHistoName] = self.histos[oldHistoName].Clone(newHistoName)
                        self.histos[newHistoName].SetTitle(newHistoName)
                    else:
                        self.histos[newHistoName].Add(self.histos[oldHistoName])
                    del self.histos[oldHistoName]

                ### now do systs - use the 1st one of the list to get the syst prototype
                ### NOTE: it supposes that all the histos grouped together have the same syst
                protoName = makeHistoName(partsList[0], sel, var)
                allSysts = matchInDictionary(self.histos, protoName+'_*')
                allSysts = [x.replace(protoName+'_', '') for x in allSysts]
                # if 'TT' in partsList[0] and var == 'MT2' and 'defaultBtagLLNoIsoBBTTCut_SR' in sel: print protoName, allSysts
                for syst in allSysts:
                    newHistoName = makeHistoName(newName, sel, var, syst)
                    for idx, oldName in enumerate(partsList):
                        oldHistoName = makeHistoName(oldName, sel, var, syst)
                        if idx == 0:
                            self.histos[newHistoName] = self.histos[oldHistoName].Clone(newHistoName)
                            self.histos[newHistoName].SetTitle(newHistoName)
                        else:
                            self.histos[newHistoName].Add(self.histos[oldHistoName])
                        del self.histos[oldHistoName]

        ## 2D
        for sel in self.selections:
            for var in self.variables2D:
                var1 = var.rsplit(':', 1)[0]
                var2 = var.rsplit(':', 1)[1]
                newHistoName = makeHisto2DName(newName, sel, var1, var2)
                
                if newHistoName in self.histos2D:
                    print "** ERROR: while grouping together histos, name", newName, "already exists in input as", newHistoName, "cannot continue..."
                    raise ValueError

                for idx, oldName in enumerate(partsList):
                    oldHistoName = makeHisto2DName(oldName, sel, var1, var2)
                    if idx == 0:
                        self.histos2D[newHistoName] = self.histos2D[oldHistoName].Clone(newHistoName)
                        self.histos2D[newHistoName].SetTitle(newHistoName)
                    else:
                        self.histos2D[newHistoName].Add(self.histos2D[oldHistoName])
                    del self.histos2D[oldHistoName]

                ### now do systs - FIXME
                ### now do systs - use the 1st one of the list to get the syst prototype
                ### NOTE: it supposes that all the histos grouped together have the same syst
                protoName = makeHisto2DName(partsList[0], sel, var1, var2)
                allSysts = matchInDictionary(self.histos2D, protoName)
                allSysts = [x.replace(protoName+'_', '') for x in allSysts]
                # print allSysts
        
        ## replace entries in data/sig/bkg with their sum
        theList = None
        if partsList[0]   in self.bkgs: theList = self.bkgs
        elif partsList[0] in self.data: theList = self.data
        elif partsList[0] in self.sigs: theList = self.sigs

        for sample in partsList:
            theList.remove(sample) ## ok since all occurences are unique
        theList.append(newName)


    # def setShape(self, fromSel, toSel, fromSample, toSample=None, fromVar='*', toVar=None, forSysts='*', do2D=False):
    #     """ take the shape from a (sample, selection, variable) and use it in another (sample, selection, variable)
    #     None values mean that the same as fromX field is used"""

    #     if not toSample:
    #         toSample = fromSample
    #     if not toVar:
    #         toVar = fromVar

    #     ## 1D
    #     for sel in self.selections:
    #         if not fnmatch.fnmatch(sel, toSel): continue
            
    #         for var in self.variables:
    #             if not fnmatch.fnmatch(var, toVar): continue
            
    #             for sample in self.samples:
    #                 if not fnmatch.fnmatch(sample, toSample): continue

    #                 hTo = self.histos[makeHistoName(sample, sel, var)]
    #                 hFrom = self.histos[makeHistoName(sample, sel, var)]


    ## FIXME: how to treat systematics properly ? Do we need to do ann alternative QCD histo for every syst?
    def makeQCD (self, SR, yieldSB, shapeSB, SBtoSRfactor, doFitIf='False', fitFunc='[0] + [1]*x', QCDname='QCD', removeNegBins = True):
        
        print "... building QCD w/ name:", QCDname, ". SR:" , SR, " yieldSB:", yieldSB, " shapeSB:", shapeSB, " SBtoSRfactor:", SBtoSRfactor
        print "    >> doFitIf:", doFitIf , "fitFunction:", fitFunc , '\n'
        
        for var in self.variables:
            for sel in self.sel_def:
                
                # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' : print "DOING ", var, sel

                ## make shape hist
                for idx, data in enumerate(self.data):
                    hname = makeHistoName(data, sel+'_'+shapeSB, var)
                    if idx == 0:
                        hQCD = self.histos[hname].Clone(makeHistoName(QCDname, sel+'_'+SR, var)) ## use SR name as this is where QCD refers to
                        hQCD.SetTitle(hQCD.GetName())
                    else:
                        hQCD.Add(self.histos[hname])
                    # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' : print ">> data - SHAPE: " , hname, hQCD.Integral()

                # subtract bkg
                for bkg in self.bkgs:
                    hname = makeHistoName(bkg, sel+'_'+shapeSB, var)
                    hQCD.Add(self.histos[hname], -1.)
                    # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' : print ">> -- bkg - SHAPE: " , hname, hQCD.Integral()

                ## remove negative bins if needed
                if removeNegBins:
                    for ibin in range(1, hQCD.GetNbinsX()+1):
                        if hQCD.GetBinContent(ibin) < 0:
                            hQCD.SetBinContent(ibin, 1.e-6)

                ## now compute yield to be set
                for idx, data in enumerate(self.data):
                    hname = makeHistoName(data, sel+'_'+yieldSB, var)
                    if idx == 0:
                        hyieldQCD = self.histos[hname].Clone(makeHistoName(QCDname+'yield', sel+'_'+yieldSB, var))
                    else:
                        hyieldQCD.Add(self.histos[hname])
                    # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' : print ">> data: " , hname, qcdYield

                for bkg in self.bkgs:
                    hname = makeHistoName(bkg, sel+'_'+yieldSB, var)
                    hyieldQCD.Add(self.histos[hname], -1)
                    # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' : print ">> -- bkg: " , hname, qcdYield

                # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' :  print qcdYield
                ## now scale
                qcdYield = hyieldQCD.Integral()
                sc = SBtoSRfactor*qcdYield/hQCD.Integral() if hQCD.Integral() > 0 else 0.0
                hQCD.Scale(sc)

                ## add to collection
                # if var == 'MT2' and sel == 'defaultBtagLLNoIsoBBTTCut' :  print 'saving as ' , hQCD.GetName()

                if eval(doFitIf):
                    try:
                        self.QCDfits
                    except AttributeError:
                        self.QCDfits = collections.OrderedDict()
                        self.QCDfitresults = collections.OrderedDict()
                    
                    # the previous QCD histogram is still stored but as 'uncorr', and the new one gets QCDname
                    hQCD.SetName('uncorr'+hQCD.GetName())
                    hQCDCorr = hQCD.Clone(makeHistoName(QCDname, sel+'_'+SR, var))
                    hQCDCorr.Scale(1./hQCDCorr.Integral())
                    hQCDNum = hyieldQCD.Clone(makeHistoName(QCDname+'FIT', sel+'_'+SR, var))
                    hQCDNum.Scale(1./hQCDNum.Integral()) ## both num and denom are normalized to 1
                    
                    hQCDNum.Divide(hQCDCorr)
                    fFitFunc = ROOT.TF1('QCDFitFunc', fitFunc, hQCD.GetXaxis().GetXmin(), hQCD.GetXaxis().GetXmax())
                    fitresult = hQCDNum.Fit(fFitFunc, "S")
                    fitresult.SetName(makeHistoName(QCDname+'fitresult', sel+'_'+SR, var))

                    self.QCDfits[hQCDNum.GetName()] = hQCDNum
                    self.QCDfitresults[fitresult.GetName()] = fitresult
                    # the fit gets attached to the histo, so that one can retrieve parameters as TF1* f = histo->GetFunction("QCDFitFunc")

                    ## scale the QCD template according to the fit function at the bin center
                    normaliz = hQCD.Integral()
                    hQCDCorr.Multiply(fFitFunc) # NOTE: multiplication is done in the function range, so it's important to set this properly before. Errors are propagated.
                    hQCDCorr.Scale(normaliz/hQCDCorr.Integral())
                    self.histos[hQCDCorr.GetName()] = hQCDCorr

                ## store hQCD - is either 'QCD' if no fit was done or uncorrQCD if fit was done, in any case is the final one to plot
                self.histos[hQCD.GetName()] = hQCD


        ### FIXME: now do 2D histos

    def rebin(self, var, sel, newbinning):        
        print '... rebinning histos for var:' , var, 'sel:', sel
        newbinning_a = array.array('d', newbinning)
        for idx, s in enumerate(self.data + self.bkgs + self.sigs):
            htorebin_name = makeHistoName(s, sel, var)
            h = self.histos[htorebin_name]
            if idx == 0: # all histos have the same binning, don't waste time
                # for i in range(1, h.GetNbinsX()+2): print var, i, h.GetBinLowEdge(i)
                oldbinning = [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)]
            if not checkBinningCompatibility (newbinning, oldbinning):
                print "** OutputManager : rebin : warning: binnings are not compatible, won't rebin", var, sel
                print "old:", oldbinning, "new:", newbinning
                continue
            h.SetName('prerebin_'+htorebin_name)
            hnew = h.Rebin(len(newbinning)-1, htorebin_name, newbinning_a)
            # print sel, var, hnew.GetNbinsX(), hnew.GetName()
            self.histos[hnew.GetName()] = hnew

            ## the following is the correct one, but slow
            # protoName = makeHistoName(s, sel, var)
            # allSysts = matchInDictionary(self.histos, protoName+'_*')
            # allSysts = [x.replace(protoName+'_', '') for x in allSysts]

            ## I should build this once
            # # now systematics
            # if not hasattr(self, 'systMap'):
            #     self.buildSystMap()

            # for the moment, only top pt rew has a systematic, let's speed up and use this information
            #### FIXME: to be generalized! (with the previous function)

            allSysts = []
            if s == 'TT':
                allSysts = ['topUp', 'topDown']
                print '.. this is a TT sample, hence the only one with systs', allSysts
            
            for syst in allSysts:
                htorebin_name = makeHistoName(s, sel, var, syst)
                h = self.histos[htorebin_name]
                h.SetName('prerebin_'+htorebin_name)
                hnew = h.Rebin(len(newbinning)-1, htorebin_name, newbinning_a)
                self.histos[hnew.GetName()] = hnew

    def saveToFile(self, fOut, saveQCDFit=True):
        fOut.cd()
        for elem in self.histos:
            self.histos[elem].Write()
        for elem in self.histos2D:
            self.histos2D[elem].Write()
        if hasattr(self, 'QCDfits') and saveQCDFit:
            for elem in self.QCDfits:
                self.QCDfits[elem].Write()
            for elem in self.QCDfitresults:
                self.QCDfitresults[elem].Write()

    # def buildSystMap(self):
    #     """ make a dictionary with all the possible systematics as suffix to speed up repeated histo lookup """

    #     ### note: there are (currently!) no syst associated to vars, and to specific regions
    #     ### so it is enough to check only the syst for every sample and for a nominal selection (sel + 'SR')

    #     self.systMap = {}

