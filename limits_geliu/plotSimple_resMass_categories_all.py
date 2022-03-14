import re, optparse
import os.path
from math import *
# from ROOT import *
import ROOT

#####

def redrawBorder():
   # this little macro redraws the axis tick marks and the pad border lines.
   ROOT.gPad.Update();
   ROOT.gPad.RedrawAxis();
   l = ROOT.TLine()
   l.SetLineWidth(3)
   l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax());
   l.DrawLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax());
   l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax());
   l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin());

def getExpValue( kl,  yt): 
    BR =1 
    return (2.09*yt*yt*yt*yt +   0.28*yt*yt*kl*kl  -1.37*yt*yt*yt*kl)*2.44477/BR;


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
        return -1.0
    else:
        return matches[-1]


c1 = ROOT.TCanvas("c1", "c1", 650, 500)
c1.SetFrameLineWidth(3)
c1.SetBottomMargin (0.15)
c1.SetRightMargin (0.05)
c1.SetLeftMargin (0.15)
c1.SetGridx()
c1.SetGridy()

mg = ROOT.TMultiGraph()

Grexp=[]

colors=[ROOT.kBlue+2,ROOT.kAzure+10,ROOT.kViolet]

var = 'DNNoutSM_kl_1'

tags = ["cards_TauTauLegacy2018_25Mar2021_m","cards_MuTauLegacy2018_25Mar2021_m","cards_ETauLegacy2018_25Mar2021_m"]
#tags = ['cards_ETauLegacy2018_test','cards_MuTauLegacy2018_test','cards_TauTauLegacy2018_test']
#tags = ['cards_ETauLegacy2018_test_noDRcut','cards_MuTauLegacy2018_test_noDRcut','cards_TauTauLegacy2018_test_noDRcut']
#tag='CombChan2018_res_17Feb2021'

channels = ["TauTau","MuTau","ETau"]
#selections = ["sboostedLLMcut"]
selections = ["s1b1jresolvedMcut", "s2b0jresolvedMcut", "sboostedLLMcut", "VBFloose"]
#selections = ["comb_cat"]

massval = [250, 260, 270, 300, 320, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
print massval

legend = ROOT.TLegend(0,0,0,0)         

legend.SetX1(0.5)
legend.SetY1(0.2)
legend.SetX2(0.88)
legend.SetY2(0.4)

legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.035)
legend.SetBorderSize(0)
# legend
legend.SetHeader('95% CL upper limits, Median expected')

### read the scan with normal width
for i in range(len(tags)):
    tag=tags[i]
   
    gr2sigma = ROOT.TGraphAsymmErrors()
    gr1sigma = ROOT.TGraphAsymmErrors()
    grexp = ROOT.TGraph()
    grobs = ROOT.TGraph()

    ptsList = [] # (x, obs, exp, p2s, p1s, m1s, m2s)

    for mass in massval:
        fName = tag+'/comb_cat/combined_out/comb.Radion{0}.log'.format(mass)
       #fName = 'cards_'+tag+'/'+sel+var+'/comb.{1}.Radion{0}.log'.format(mass,sel)


        exp   = 10.*parseFile(fName)            
        obs   = exp#parseFile(fName, exp=False) 
        m1s_t = 10.*parseFile(fName, CL='16.0') 
        p1s_t = 10.*parseFile(fName, CL='84.0') 
        m2s_t = 10.*parseFile(fName, CL=' 2.5') 
        p2s_t = 10.*parseFile(fName, CL='97.5') 

        ## because the other code wants +/ sigma vars as deviations, without sign, from the centeal exp value...
        p2s = p2s_t - exp
        p1s = p1s_t - exp
        m2s = exp - m2s_t
        m1s = exp - m1s_t
        xval = mass

        ptsList.append((xval, obs, exp, p2s, p1s, m1s, m2s))

    # print lambdasfiner

    # grexp.SetPoint(ipt, xval, exp)
    # grobs.SetPoint(ipt, xval, obs)
    # gr1sigma.SetPoint(ipt, xval, exp)
    # gr2sigma.SetPoint(ipt, xval, exp)
    # gr1sigma.SetPointError(ipt, 0,0,m1s,p1s)
    # gr2sigma.SetPointError(ipt, 0,0,m2s,p2s)
    ptsList.sort()
    for ipt, pt in enumerate(ptsList):
        xval = pt[0]
        obs  = pt[1]
        exp  = pt[2]
        p2s  = pt[3]
        p1s  = pt[4]
        m1s  = pt[5]
        m2s  = pt[6]
        print xval, exp
        grexp.SetPoint(ipt, xval, exp)
        grobs.SetPoint(ipt, xval, obs)
        gr1sigma.SetPoint(ipt, xval, exp)
        gr2sigma.SetPoint(ipt, xval, exp)
        gr1sigma.SetPointError(ipt, 0,0,m1s,p1s)
        gr2sigma.SetPointError(ipt, 0,0,m2s,p2s)
    Grexp.append(grexp)
    
    ######## set styles
    grexp.SetMarkerStyle(24)
    grexp.SetMarkerColor(4)
    grexp.SetMarkerSize(0.8)
    grexp.SetLineColor(colors[i])
    grexp.SetLineWidth(3)
    grexp.SetLineStyle(2)
    grexp.SetFillColor(0)

    if 'ETau' in tag:
	ch="bb e#tau_{h}"
    elif 'MuTau' in tag:
	ch="bb #mu_{}#tau_{h}"
    elif 'TauTau' in tag:
	ch="bb #tau_{h}#tau_{h}"
    legend.AddEntry(Grexp[i], ch, "l")


##### text
# pt = ROOT.TPaveText(0.1663218,0.886316,0.3045977,0.978947,"brNDC")
pt = ROOT.TPaveText(0.1663218-0.02,0.886316,0.3045977-0.02,0.978947,"brNDC")
pt.SetBorderSize(0)
pt.SetTextAlign(12)
pt.SetTextFont(62)
pt.SetTextSize(0.05)
pt.SetFillColor(0)
pt.SetFillStyle(0)
pt.AddText("CMS #font[52]{Internal}" )
#pt.AddText("CMS" )
# pt.AddText("#font[52]{preliminary}")
pt2 = ROOT.TPaveText(0.736111,0.9066667,0.847222,0.954641,"brNDC")
pt2.SetBorderSize(0)
pt2.SetFillColor(0)
pt2.SetTextSize(0.040)
pt2.SetTextFont(42)
pt2.SetFillStyle(0)
if '2016' in tag:
   pt2.AddText("2016 - 35.9 fb^{-1} (13 TeV)")
elif '2017' in tag:
   pt2.AddText("2017 - 41.6 fb^{-1} (13 TeV)")
elif '2018' in tag:
   pt2.AddText("2018 - 59.7 fb^{-1} (13 TeV)")
else:
   pt2.AddText("Run2 - 137.1 fb^{-1} (13 TeV)")
   
pt4 = ROOT.TPaveText(0.4819196+0.036,0.7780357+0.015+0.02,0.9008929+0.036,0.8675595+0.015,"brNDC")
pt4.SetTextAlign(12)
pt4.SetFillColor(ROOT.kWhite)
pt4.SetFillStyle(1001)
pt4.SetTextFont(42)
pt4.SetTextSize(0.05)
pt4.SetBorderSize(0)
pt4.SetTextAlign(32)
if 'ETau' in tag:
   pt4.AddText("bb e#tau_{h}") 
elif 'MuTau' in tag:
   pt4.AddText("bb #mu_{}#tau_{h}")
elif 'TauTau' in tag:
   pt4.AddText("bb #tau_{h}#tau_{h}")
   #pt4.AddText("bb #tau_{#mu}#tau_{h} + bb #tau_{e}#tau_{h} + bb #tau_{h}#tau_{h}") 
   #pt4.AddText("bb #tau_{e}#tau_{h}") 
   #pt4.AddText(sel + ", bb #tau_{h}#tau_{h}")
   #pt4.AddText("GGF-HH Comb. cat.")
   #pt4.AddText("HH #rightarrow bb#tau#tau")

offs = 0.020
height = 0.05
pt5 = ROOT.TPaveText(0.4819196+0.036+0.10,0.7780357+0.015-offs,0.9008929+0.036,0.7780357+0.015-offs+height,"brNDC")
pt5.SetTextAlign(12)
pt5.SetFillColor(ROOT.kWhite)
pt5.SetFillStyle(1001)
pt5.SetTextFont(42)
pt5.SetTextSize(0.05)
pt5.SetBorderSize(0)
pt5.SetTextAlign(32)
# pt5.AddText("bb #mu_{}#tau_{h} + bb e#tau_{h} + bb #tau_{h}#tau_{h}") 
# pt5.AddText("bb #tau_{#mu}#tau_{h} + bb #tau_{e}#tau_{h} + bb #tau_{h}#tau_{h}") 
#pt5.AddText('#scale[0.8]{{0}}'.format(sel))
sel="Combine categories"
pt5.AddText(sel)


xmin=200
xmax=3100

   


hframe = ROOT.TH1F('hframe', '', 100, 250, 3100)
hframe.SetMinimum(0.1)
if '2016' in tag:
   hframe.SetMaximum(10000)
   hframe.SetMinimum(0.1)#5000)
elif '2017' in tag:
   hframe.SetMaximum(10000)
   hframe.SetMinimum(0.1)#5000)
elif '2018' in tag:
   hframe.SetMinimum(0.1)
   hframe.SetMaximum(10000)#4000)
else:
   hframe.SetMaximum(10000)
   hframe.SetMinimum(0.1)#5000)
   
hframe.GetYaxis().SetTitleSize(0.047)
hframe.GetXaxis().SetTitleSize(0.055)
hframe.GetYaxis().SetLabelSize(0.045)
hframe.GetXaxis().SetLabelSize(0.045)
hframe.GetXaxis().SetLabelOffset(0.012)
hframe.GetYaxis().SetTitleOffset(1.2)
hframe.GetXaxis().SetTitleOffset(1.1)

hframe.GetYaxis().SetTitle("95% CL on #sigma #times #bf{#it{#Beta}}(S#rightarrowHH#rightarrow bb#tau#tau) [fb]")
#hframe.GetYaxis().SetTitle("95% CL on #sigma#times B (S#rightarrow HH#rightarrow#tau#tau) [fb]")
hframe.GetXaxis().SetTitle("m_{S} [GeV]")

hframe.SetStats(0)
ROOT.gPad.SetTicky()
hframe.Draw()

# mg.Draw("pmc plc same")
#gr2sigma.Draw("3same")
#gr1sigma.Draw("3same")
for i in range(3):
   Grexp[i].Draw("Lsame")
#grobs.Draw("Lsame")

#graph.Draw("l same") #LP NOTE: THIS SHOULD BE THE THEORY GRAPH

#graph2.Draw("l same")
#Graph_syst_Scale.Draw("e3 same");
#Graph_syst_Scale2.Draw("e3 same");

pt.Draw()
pt2.Draw()
# pt4.Draw()
#txt_kt1.Draw()
#txt_kt2.Draw()
#redrawBorder()
c1.Update()
c1.RedrawAxis("g")
c1.SetLogy()
legend.Draw()
#pt4.Draw()
pt5.Draw()

c1.Update()
#raw_input()

print 'Limit for ', sel
smXS = 31.05*0.073



#c1.Print("plots/klscan_"+tag+".pdf", 'pdf')
c1.Print("plots/massScan_Legacy2018_25Mar2021_combinedcategories.pdf",'pdf')
