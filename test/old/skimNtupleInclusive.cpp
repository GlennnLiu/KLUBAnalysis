#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBranch.h"
#include "TString.h"
#include "TLorentzVector.h"

// bigTree is produced on an existing ntuple as follows (see at the end of the file)
#include "bigTree.h"
#include "smallTree.h"
#include "OfflineProducerHelper.h"
#include "PUReweight.h"
#include "triggerReader.h"
#include "bJetRegrVars.h"
#include "bTagSF.h"
//#include "../../HHKinFit/interface/HHKinFitMaster.h"
#include "../../HHKinFit2/include/HHKinFitMasterHeavyHiggs.h"
//#include "../../HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "ScaleFactor.h"
#include "ConfigParser.h"
//#include "../../HHKinFit2/include/exceptions/HHInvMConstraintException.h"
//#include "../../HHKinFit2/include/exceptions/HHEnergyRangeException.h"
#include "exceptions/HHInvMConstraintException.h"
#include "exceptions/HHEnergyRangeException.h"
#include "exceptions/HHEnergyConstraintException.h"

#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

//
using namespace std ;

const double aTopRW=0.156;
const double bTopRW= -0.00137;
const float DYscale_LL[3] = {1.17835, 1.80015, 0.809161} ; // computed from fit for LL and MM b tag
const float DYscale_MM[3] = {1.20859, 1.0445 , 1.54734 } ;

// const float stitchWeights [5] = {1.11179e-7, 3.04659e-9, 3.28633e-9, 3.48951e-9, 2.5776e-9} ; // weights DY stitch in njets, to be updated at each production (depend on n evtsn processed)
const float stitchWeights [5] = {11.55916, 0.316751, 0.341677, 0.362801, 0.267991} ; // weights DY stitch in njets, to be updated at each production (depend on n evtsn processed)
/*  the modern way of making lorenzvectors (no warnings)
    here is an example of usage:

    https://github.com/govoni/LHEAnalysis/blob/master/epem_01.cpp

    here are the headers needed:

    #include "Math/Vector3D.h"
    #include "Math/Vector4D.h"

    using namespace ROOT::Math;

    typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;

    I am missing about how to set one of those things up with Pt, eta, phi coordinates

*/

TFile *fileTO = TFile::Open("/home/llr/cms/cadamuro/TauTagAndProbe/CMSSW_7_6_3/src/TauTagAndProbe/TauTagAndProbe/test/turnOn_2GevIso.root");
TH1F *hTurnOn = (TH1F*)fileTO->Get("hTurnOn");
TFile *fileTOTrig = TFile::Open("/home/llr/cms/ortona/diHiggs/CMSSW_7_4_7/src/KLUBAnalysis/studies/MakeTurnOn/turnOnTrig05.root");
TH1F *hTurnOnTrig = (TH1F*)fileTOTrig->Get("hTurnOn");

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
// open input txt file and append all the files it contains to TChain
void appendFromFileList (TChain* chain, TString filename)
{
  //cout << "=== inizio parser ===" << endl;
  std::ifstream infile(filename.Data());
  std::string line;
  while (std::getline(infile, line))
  {
    line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
    while (line.find(" ") != std::string::npos) line = line.erase(line.find(" "), 1); // remove white spaces
    while (line.find("\n") != std::string::npos) line = line.erase(line.find("\n"), 1); // remove new line characters
    while (line.find("\r") != std::string::npos) line = line.erase(line.find("\r"), 1); // remove carriage return characters
    if (!line.empty()) // skip empty lines
      chain->Add(line.c_str());
  }
  return;
}


float turnOnCB(float x, float m0, float sigma, float alpha, float n, float norm){
  float sqrtPiOver2 = TMath::Sqrt(TMath::PiOver2());
  float sqrt2 = TMath::Sqrt(2.);
  float sig = abs(sigma);
  float t = (x - m0)/sig * alpha / abs(alpha);
  float absAlpha = abs(alpha/sig);
  float a = TMath::Power(n/absAlpha, n) * TMath::Exp(-0.5 * absAlpha * absAlpha);
  float b = absAlpha - n/absAlpha;
  float arg = absAlpha / sqrt2;
  float ApproxErf = TMath::Erf(arg);
  if (arg > 5.) ApproxErf = 1.;
  else if(arg < -5.) ApproxErf = -1.;
  float leftArea = (1. + ApproxErf) * sqrtPiOver2;
  float rightArea = ( a * 1./TMath::Power(absAlpha-b, n-1) ) / (n-1);
  float area = leftArea + rightArea;
  if (t <= absAlpha){
    arg = t / sqrt2;
    if (arg > 5.) ApproxErf = 1.;
    else {
      if(arg < -5.) ApproxErf = -1.;
      else ApproxErf = TMath::Erf(arg);
    }
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1-n)) / area;
}


float turnOnSF(float pt,int genInd=66615){
  //return 1.0/turnOnCB(pt,3.60274e+01,5.89434e+00,5.82870e+00,1.83737e+00,9.58000e-01)*turnOnCB(pt,3.45412e+01,5.63353e+00,2.49242e+00,3.35896e+00,1);
  //return turnOnCB(pt,3.45412e+01,5.63353e+00,2.49242e+00,3.35896e+00,1);
  //if(genInd>=0)
  float threshold = 120;
  if(pt>threshold) pt = threshold;

//Riccardo
  if(false){
    if(genInd!=66615)
      return turnOnCB(pt,3.94747E+01, 7.23546E+00 ,1.08089E+01 ,1.33930E+00, 1.00000E+00);
    else
      return turnOnCB(pt,3.85953E+01,5.74632E+00 ,5.08553E+00, 5.45593E+00, 9.42168E-01);
  }
  //return turnOnCB(pt,3.60274e+01,5.89434e+00,5.82870e+00,1.83737e+00,9.58000e-01);
//Luca
  return hTurnOn->GetBinContent(hTurnOn->FindBin(pt));//*TMath::Sqrt(hTurnOnTrig->GetBinContent(hTurnOnTrig->FindBin(pt)));

  //TF1 func("fakefunc","expo(0)",20,400);
  //func.SetParameter(0,-0.68459);
  //func.SetParameter(1,-0.000944352);
  /*
    double p[5];
    p[0]                        =   0.00742864;
    p[1]                        =  0.000167935;
    p[2]                        = -7.51164e-07;
    p[3]                        =    1.417e-08;
    p[4]                        = -6.09917e-11;
    TF1 func("fakefunc","pol4(0)",20,400);
    func.SetParameters(p);
    return func.Eval(pt);*/
}
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -
// open the first file in the input list, retrieve the histogram "Counters" for the trigger names and return a copy of it
TH1F* getFirstFileHisto (TString filename, bool isForTriggers=true)
{
  std::ifstream infile(filename.Data());
  std::string line;
  while (std::getline(infile, line))
  {
    line = line.substr(0, line.find("#", 0)); // remove comments introduced by #
    while (line.find(" ") != std::string::npos) line = line.erase(line.find(" "), 1); // remove white spaces
    while (line.find("\n") != std::string::npos) line = line.erase(line.find("\n"), 1); // remove new line characters
    while (line.find("\r") != std::string::npos) line = line.erase(line.find("\r"), 1); // remove carriage return characters
    if (!line.empty()) // skip empty lines
      break;
  }

  TFile* fIn = TFile::Open (line.c_str());
  TH1F* dummy = (TH1F*) fIn->Get ("HTauTauTree/Counters");
  TString name = "Counters_perTrigger";
  if(!isForTriggers) {
    dummy = (TH1F*) fIn->Get ("HTauTauTree/TauIDs");
    name = "Counters_pertauID";
  }
  TH1F* histo = new TH1F (*dummy);
  histo-> SetDirectory(0);
  histo->SetName (name.Data());
  fIn->Close();
  return histo;
}
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -


bool isDegenere (TLorentzVector & first, TLorentzVector & second)
{
  //if ( fabs(first.Px()-second.Px())+fabs(first.Py()-second.Py())+fabs(first.Pz()-second.Pz())<0.03 * (fabs(first.Px())+fabs(first.Py())+fabs(first.Pz()) )) return true;
  if (first.DeltaR(second) < 0.03) return true;

  return false ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- -


float deltaPhi (float phi1, float phi2)
{
  float deltaphi = fabs (phi1 - phi2) ;
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308 ;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi ;
  return deltaphi ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


struct leptSort: public std::binary_function<pair<TLorentzVector, float> &, pair<TLorentzVector, float> &, bool>
{
  bool operator() (pair<TLorentzVector, float> & x, pair<TLorentzVector, float> & y)
    {
      return x.first.Pt () < y.first.Pt () ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


template <class T>
struct scoreSortSingle: public std::binary_function<pair <T, float> &, pair <T, float> &, bool>
{
  bool operator() (pair <T, float> & x, pair <T, float> & y)
    {
      return x.second < y.second ;
    }
} ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


pair<int, int>
chooseHighestPtJets (vector <pair <int, float> > & jets_and_btag)
{
  // leap of faith
  return pair<int, int> (jets_and_btag.at (0).first, jets_and_btag.at (1).first) ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


pair<int, int>
chooseHighestBtagJets (vector <pair <int, float> > & jets_and_btag)
{
  sort (jets_and_btag.rbegin (), jets_and_btag.rend (), scoreSortSingle<int> ()) ;

  int firstBjetIndex = jets_and_btag.at (0).first ;
  int secondBjetIndex = jets_and_btag.at (1).first ;
  return pair<int, int> (firstBjetIndex, secondBjetIndex) ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


/**
   muons = 0
   electrons = 1
   taus = 2
*/
float getIso (unsigned int iDau, float pt, bigTree & theBigTree)
{
  int type = theBigTree.particleType->at (iDau) ;
  // is tauH
  if (type == 2)
    return theBigTree.daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(iDau) ;
  // muon
  if (type == 1 || type == 0)
    return theBigTree.combreliso->at(iDau);

  return -1 ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char** argv)
{
  // read input file and cfg
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  TRandom3 *rnd = new TRandom3(0);

  if (argc < 7)
  {
    cerr << "missing input parameters" << endl ;
    cerr << "usage: " << argv[0]
	 << " inputFileNameList outputFileName crossSection isData configFile runHHKinFit" << endl
	 << "OPTIONAL: xsecScale(stitch) HTMax(stitch) isTTBar DY_Nbs" << endl ;
    exit (1) ;
  }
  TString inputFile = argv[1] ;
  TString outputFile = argv[2] ;
  float XS = atof (argv[3]) ;
  bool isMC = true;
  int isDatabuf = atoi (argv[4]);
  if (isDatabuf == 1)
  {
    cout << "RUNNING ON DATA\n" ;
    isMC = false;
    XS = 1.;
  }

  if (gConfigParser) return 1 ;
  gConfigParser = new ConfigParser () ;

  TString config ;
  config.Form ("%s",argv[5]) ;

  bool runHHKinFit = false;
  string opt7 (argv[6]);
  if (opt7 == "1") runHHKinFit = true;

  cout << "Running HH in fit: " << runHHKinFit << " " << argv[6] << endl;

  // optional parameters for stitching
  float xsecScale = 1.0;
  float HTMax = -999.0;
  bool isTTBar = false;
  bool DY_Nbs = false; // run on genjets to count in DY samples the number of b jets
  bool DY_tostitch = false;

  if (argc >= 8)
  {
    xsecScale = atof (argv[7]);
    cout << " ** INFO: doing stitching, xsec scaled by: " << xsecScale << endl;
    XS = XS*xsecScale;

    if (argc >= 9)
    {
      HTMax = atof(argv[8]);
      cout << " ** INFO: in stitching, removing HT < " << HTMax << endl;

      if (argc >= 10)
      {
        int isTTBarI = atoi(argv[9]);
        if (isTTBarI == 1) isTTBar = true;
        cout << " ** INFO: is this a TTbar sample? : " << isTTBar << endl;
        if (argc >= 11)
        {
          int I_DY_Nbs = atoi(argv[10]);
          if (I_DY_Nbs == 1)
          {
            DY_Nbs = true;
            DY_tostitch = true; // FIXME!! this is ok only if we use jet binned samples
          }
          // if (argc >= 12)
          // {
          //   int I_DY_tostitch = atoi(argv[11]);
          //   if (I_DY_tostitch == 1) DY_tostitch = true;
          // }
        }
      }
    }
  }

  // force to !TTbar if isData -- you never know...
  if (!isMC) isTTBar = false;
  cout << " ** INFO: going to loop on genjets to find number of b? " << DY_Nbs << endl;


  cout << " ** INFO: running on file list : " << inputFile << endl;
  cout << " ** INFO: saving output in     : " << outputFile << endl;


  if (! (gConfigParser->init (config)))
  {
    cout << ">>> parseConfigFile::Could not open configuration file " << config << endl ;
    return -1 ;
  }

  bool   beInclusive         = gConfigParser->readBoolOption   ("selections::beInclusive") ;
  float  PUjetID_minCut      = gConfigParser->readFloatOption  ("parameters::PUjetIDminCut") ;
  int    saveOS              = gConfigParser->readIntOption    ("parameters::saveOS") ;
  float  lepCleaningCone     = gConfigParser->readFloatOption  ("parameters::lepCleaningCone") ;
  int    bChoiceFlag         = gConfigParser->readFloatOption  ("parameters::bChoiceFlag") ;
  int    PUReweight_MC       = gConfigParser->readFloatOption  ("parameters::PUReweightMC") ;
  int    PUReweight_target   = gConfigParser->readFloatOption  ("parameters::PUReweighttarget") ;
  string leptonSelectionFlag = gConfigParser->readStringOption ("parameters::lepSelections") ;
  int maxNjetsSaved          = gConfigParser->readIntOption    ("parameters::maxNjetsSaved") ;

  vector<string> trigMuTau   =  (isMC ? gConfigParser->readStringListOption ("triggersMC::MuTau")  : gConfigParser->readStringListOption ("triggersData::MuTau")) ;
  vector<string> trigTauTau   = (isMC ? gConfigParser->readStringListOption ("triggersMC::TauTau") : gConfigParser->readStringListOption ("triggersData::TauTau")) ;
  vector<string> trigEleTau   = (isMC ? gConfigParser->readStringListOption ("triggersMC::EleTau") : gConfigParser->readStringListOption ("triggersData::EleTau")) ;
  vector<string> trigEleMu   =  (isMC ? gConfigParser->readStringListOption ("triggersMC::EleMu")  : gConfigParser->readStringListOption ("triggersData::EleMu")) ;
  //I didn't store MuMu and I don't care for eleele
  vector<string> trigEleEle   =  (isMC ? gConfigParser->readStringListOption ("triggersMC::EleEle")  : gConfigParser->readStringListOption ("triggersData::EleEle")) ;
  vector<string> trigMuMu   =  (isMC ? gConfigParser->readStringListOption ("triggersMC::MuMu")  : gConfigParser->readStringListOption ("triggersData::MuMu")) ;

  cout << "TRIGGERS: " << endl;

  cout << "MuTau" << endl; cout << "  --> ";
  for (unsigned int i = 0 ; i < trigMuTau.size(); i++) cout << "  " << trigMuTau.at(i);
  cout << endl;

  cout << "EleTau" << endl; cout << "  --> ";
  for (unsigned int i = 0 ; i < trigEleTau.size(); i++) cout << "  " << trigEleTau.at(i);
  cout << endl;

  cout << "TauTau" << endl; cout << "  --> ";
  for (unsigned int i = 0 ; i < trigTauTau.size(); i++) cout << "  " << trigTauTau.at(i);
  cout << endl;


  int skipTriggers = 0; //1: reweight, 0: skip, -1: dont skip
  if (gConfigParser->isDefined ("debug::skipTriggers"))
  {
    skipTriggers = gConfigParser->readIntOption ("debug::skipTriggers");
  }
  cout << "skipTriggers? " << skipTriggers << endl;
  TFile *fTurnOn;
  TGraph * gTurnOn;
  if(skipTriggers>0){
    fTurnOn = TFile::Open("config/tauLegEff.root");
    gTurnOn = (TGraph*)fTurnOn->Get("tauLegEff");
  }
  bool skipIso = true;
  if (gConfigParser->isDefined ("debug::skipIso"))
  {
    skipIso = gConfigParser->readBoolOption ("debug::skipIso");
  }

  int strategy = 2016;
  if (gConfigParser->isDefined ("debug::strategy"))
  {
    strategy = gConfigParser->readIntOption ("debug::strategy");
  }


  TH1F *hFakeSF_num = new TH1F("hFSF_num","hFSF_num",40,20,200);
  TH1F *hFakeSF_den = new TH1F("hFSF_den","hFSF_den",40,20,200);

  string bRegrWeights("");
  bool computeBregr = gConfigParser->readBoolOption ("bRegression::computeBregr");
  if (computeBregr) bRegrWeights = gConfigParser->readStringOption("bRegression::weights");

  // input and output setup
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  TChain * bigChain = new TChain ("HTauTauTree/HTauTauTree") ;
  //bigChain->Add (inputFile) ;
  appendFromFileList (bigChain, inputFile);
  bigTree theBigTree (bigChain) ;

  //Create a new file + a clone of old tree header. Do not copy events
  TFile * smallFile = new TFile (outputFile, "recreate") ;
  smallFile->cd () ;
  smallTree theSmallTree ("HTauTauTree") ;

  // these are needed for the HHKinFit
  //vector<Int_t> hypo_mh1 ; //FIXME why is this an integer?!
  //hypo_mh1.push_back (125) ;
  //vector<Int_t> hypo_mh2 ;
  //hypo_mh2.push_back (125) ;
  int hypo_mh1=125,hypo_mh2=125;

  //int eventsNumber = theBigTree.fChain->GetEntries () ;

  float totalEvents = 0. ;
  float selectedEvents = 0. ;

  int totalNoWeightsEventsNum = 0 ;
  int selectedNoWeightsEventsNum = 0 ;

  int selectionsNumber = 3 ;
  vector<float> counter (selectionsNumber + 1, 0.) ;

  // ------------------------------

  TH1F* hTriggers = getFirstFileHisto (inputFile);
  TH1F* hTauIDS = getFirstFileHisto (inputFile,false);
  triggerReader trigReader (hTriggers);
  trigReader.addTauTauTrigs (trigTauTau);
  trigReader.addMuTauTrigs  (trigMuTau);
  trigReader.addEleTauTrigs (trigEleTau);
  trigReader.addMuEleTrigs  (trigEleMu);
  trigReader.addMuMuTrigs   (trigMuMu);
  trigReader.addEleEleTrigs (trigEleEle);

  // ------------------------------

  OfflineProducerHelper oph (hTriggers, hTauIDS) ;

  // ------------------------------

  bJetRegrVars bjrv;
  TMVA::Reader *bRreader = new TMVA::Reader();
  bjrv.setReader (bRreader);
  string bRegrMethodName = "BDTG method";
  if (computeBregr)
    bRreader->BookMVA( bRegrMethodName.c_str(), bRegrWeights.c_str() );
  // ------------------------------

  PUReweight reweight (PUReweight::NONE); // none : no PU reweight (always returns 1)

  // ------------------------------
  string bTag_SFFile = gConfigParser->readStringOption("bTagScaleFactors::SFFile") ;
  string bTag_effFile = gConfigParser->readStringOption("bTagScaleFactors::effFile") ;
  cout << "B Tag SF file: " << bTag_SFFile << endl;
  bTagSF bTagSFHelper (bTag_SFFile, bTag_effFile, ""); // third field unused, but could be needed to select efficiencies for different selection levels

  // ------------------------------
  ScaleFactor * myScaleFactor[2][2]; // [0: mu, 1: ele] [0: trigger, 1: ID]
  for (int i = 0 ; i < 2; i++)
    for (int j = 0; j < 2; j++)
      myScaleFactor[i][j]= new ScaleFactor();

  myScaleFactor[0][0] -> init_ScaleFactor("weights/data/Muon/Muon_IsoMu18_fall15.root");
  myScaleFactor[0][1] -> init_ScaleFactor("weights/data/Muon/Muon_IdIso0p1_fall15.root");
  myScaleFactor[1][0] -> init_ScaleFactor("weights/data/Electron/Electron_Ele23_fall15.root");
  myScaleFactor[1][1] -> init_ScaleFactor("weights/data/Electron/Electron_IdIso0p1_fall15.root");

  // loop over events
  //for (Long64_t iEvent = 0 ; iEvent < eventsNumber ; ++iEvent)
  for (Long64_t iEvent = 0 ; true ; ++iEvent)
  {
    if (iEvent % 10000 == 0)  cout << "reading event " << iEvent << endl ;
    int selID = 0 ;

    theSmallTree.clearVars () ;
    int got = theBigTree.fChain->GetEntry(iEvent);
    if (got == 0) break;
    //theBigTree.GetEntry (iEvent) ;

    // directly reject events outside HT range in case of stitching of inclusive sample-- they should not count in weights
    if (HTMax > 0)
    {
      if (theBigTree.lheHt > HTMax) continue;
    }

    float stitchWeight = 1.0;
    if (DY_tostitch)
    {
      int njets = theBigTree.lheNOutPartons;
      stitchWeight = stitchWeights[njets];
    }

    // gen info -- fetch tt pair and compute top PT reweight
    float topPtReweight = 1.0; // 1 for all the other samples
    if (isTTBar)
    {
      float ptTop1 = -1.0;
      float ptTop2 = -1.0;

      for (unsigned int igen = 0; igen < theBigTree.genpart_pdg->size(); igen++)
      {
	int pdg = theBigTree.genpart_pdg->at(igen);
	//int genflags = theBigTree.genpart_flags->at(igen);
	//int topDM = theBigTree.genpart_TopDecayMode->at(igen);

	if (abs(pdg) == 6) // top -- from LLRNtuples, I am guaranteed to have only 2 tops
	{
	  TLorentzVector TopV;
	  TopV.SetPxPyPzE (theBigTree.genpart_px->at(igen), theBigTree.genpart_py->at(igen), theBigTree.genpart_pz->at(igen), theBigTree.genpart_e->at(igen) ) ;
	  if (ptTop1 < 0) ptTop1 = TopV.Pt();
	  else if (ptTop2 < 0) { ptTop2 = TopV.Pt(); break; } // check done in paralles shows that I never have > 2 top .  break is safe .
	  //else cout << " !! skim warning: sample is declared as as ttbar, but I have > 2 gen top in the event! " << endl;
	}
      }
      if (ptTop1 < 0 || ptTop2 < 0)
      {
	cout << "Warning: event is TTbar but I didn't find 2 tops (1,2) :" << ptTop1 << " " << ptTop2 << endl;
      }
      else
      {
	float SFTop1 = TMath::Exp(aTopRW+bTopRW*ptTop1);
	float SFTop2 = TMath::Exp(aTopRW+bTopRW*ptTop2);
	topPtReweight = TMath::Sqrt (SFTop1*SFTop2); // save later together with other weights
	theSmallTree.m_topReweight = topPtReweight ;
      }
    }

    // For Drell-Yan only: loop on genjets and count how many are there with 0,1,2 b
    // 0: 0bjet, 1: 1 b jet, 2: >= 2 b jet
    theSmallTree.m_DYscale_LL = 1.0; // all the other MC samples + data have no weight
    theSmallTree.m_DYscale_MM = 1.0;

    if (isMC && DY_Nbs)
    {
      TLorentzVector vgj;
      int nbs = 0;
      for (unsigned int igj = 0; igj < theBigTree.genjet_px->size(); igj++)
      {
	vgj.SetPxPyPzE(theBigTree.genjet_px->at(igj), theBigTree.genjet_py->at(igj), theBigTree.genjet_pz->at(igj), theBigTree.genjet_e->at(igj));
	if (vgj.Pt() > 20 && TMath::Abs(vgj.Eta()) < 2.5)
	{
	  int theFlav = theBigTree.genjet_hadronFlavour->at(igj);
	  if (abs(theFlav) == 5) nbs++;

	  // about 2% of DY events print the following message :-(
	  // if (theFlav == -999) cout << "** warning: gen jet with flav = -999 of pt: " << vgj.Pt() << " eta: " << vgj.Eta() << endl;
	}
      }
      if (nbs > 2) nbs = 2;
      theSmallTree.m_nBhadrons = nbs;

      theSmallTree.m_DYscale_LL = DYscale_LL[nbs];
      theSmallTree.m_DYscale_MM = DYscale_MM[nbs];

    }

    if (isMC)
    {
      totalEvents += theBigTree.aMCatNLOweight * reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) * topPtReweight;
      counter.at (selID++) += theBigTree.aMCatNLOweight * reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) * topPtReweight;
    }
    else
    {
      totalEvents += 1 ;
      counter.at (selID++) += 1 ;
    }
    ++totalNoWeightsEventsNum ;

    if (theBigTree.indexDau1->size () == 0) continue ;
    int metbit = theBigTree.metfilterbit;
    int metpass = metbit & (1 << 0);
    metpass += metbit & (1 << 2);
    metpass += metbit & (1 << 5);
    metpass += metbit & (1 << 6);
    //if(metpass > 0) continue ;
    //WARNING!!!!!

    if (isMC) counter.at (selID++) += theBigTree.aMCatNLOweight * reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) * topPtReweight;
    else      counter.at (selID++) += 1 ;

    // assume that the ordering of the pair numbering
    // is the same as the priority we want to give:
    //  0 : muTau
    //  1 : eTau
    //  2 : TauTau
    //  3 : mumu
    //  4 : ee
    //  5 : emu

    map<int, int> foundPairs ; // pairType, pairIndexInVectors

    // determine kind of pairs according to trigger
    int trigPairType = -99;
    Long64_t triggerbit = theBigTree.triggerbit;

    // // for mumu skim only
    // if (trigReader.checkOR (3, triggerbit) ) trigPairType = 3;

    // for all the other skims
    if (trigReader.checkOR (0, triggerbit) ) trigPairType = 0; // muTau
    //if (false) trigPairType = 0; // muTau
    else if (trigReader.checkOR (1, triggerbit) ) trigPairType = 1; // etau
    else if (trigReader.checkOR (2, triggerbit) ) trigPairType = 2; // tautau
    else if (trigReader.checkOR (3, triggerbit) ) trigPairType = 3; // mumu
    else if (trigReader.checkOR (4, triggerbit) ) trigPairType = 4; // ee
    else if (trigReader.checkOR (5, triggerbit) ) trigPairType = 5; // emu

    theSmallTree.m_RunNumber = theBigTree.RunNumber ;
    float turnOnReweight = 1.0;
    vector<int> tau_daus;
    //int nmu=0,ne=0;ntau=0;
    int iPair =-1;
    if(strategy == 2016){
      int pairType =-1;
      for(int iType=0; iType<3 && iPair<0 ;iType++){
	for (unsigned int ipair = 0 ; ipair < theBigTree.indexDau1->size () ; ++ipair){
	  if (!oph.pairPassBaseline (&theBigTree, ipair, (leptonSelectionFlag + string("Iso")).c_str ())) {
	    continue ;
	  }
	  /*
	  //Gen Matching
	  int firstDaughterIndex = theBigTree.indexDau1->at (ipair) ;
	  int secondDaughterIndex = theBigTree.indexDau2->at (ipair) ;
	  int type1 = theBigTree.particleType->at (firstDaughterIndex) ;
	  int type2 = theBigTree.particleType->at (secondDaughterIndex) ;
	  pairType = oph.getPairType (type1, type2) ;
	  if(pairType != iType) continue;

	  TLorentzVector tlv_firstLepton
	  (
	  theBigTree.daughters_px->at (firstDaughterIndex),
	  theBigTree.daughters_py->at (firstDaughterIndex),
	  theBigTree.daughters_pz->at (firstDaughterIndex),
	  theBigTree.daughters_e->at (firstDaughterIndex)
	  ) ;
	  TLorentzVector tlv_secondLepton
	  (
	  theBigTree.daughters_px->at (secondDaughterIndex),
	  theBigTree.daughters_py->at (secondDaughterIndex),
	  theBigTree.daughters_pz->at (secondDaughterIndex),
	  theBigTree.daughters_e->at (secondDaughterIndex)
	  ) ;
	  bool matched1 = false, matched2=false;
	  for(int igen=0; igen<theBigTree.genpart_pdg->size(); igen++){
	  if(abs(theBigTree.genpart_pdg->at(igen))!=66615 || pairType != 2) continue;
	  TLorentzVector tlv_gen
	  (
	  theBigTree.genpart_px->at (igen),
	  theBigTree.genpart_py->at (igen),
	  theBigTree.genpart_pz->at (igen),
	  theBigTree.genpart_e->at (igen)
	  ) ;
	  if(tlv_gen.DeltaR(tlv_firstLepton)<0.5)matched1=true;
	  if(tlv_gen.DeltaR(tlv_secondLepton)<0.5)matched2=true;
	  }
	  if(!(matched1 && matched2)) continue;
	  */
	  iPair = (int)ipair;
	  break;
	}
      }
      if (iPair<0)continue;
    }
    //int itau=0,ie=0,imu=0;
    // FIXME!! ee and mumu are "eaten" by the muTau / eTau final states when we use single lepton triggers
    for (unsigned int iiPair = 0 ; iiPair < theBigTree.indexDau1->size () ; ++iiPair)
      //for(int iiiii=0;iiiii<1;iiiii++)
    {
      if(strategy==2015) {
	iPair = (int)iiPair;
	if (!oph.pairPassBaseline (&theBigTree, iPair, leptonSelectionFlag.c_str ())) continue ;
      }
      else iiPair = theBigTree.indexDau1->size ()+10;

      // FIXME need to implement here the choice of iso / anti-iso

      int firstDaughterIndex = theBigTree.indexDau1->at (iPair) ;
      int secondDaughterIndex = theBigTree.indexDau2->at (iPair) ;
      TLorentzVector tlv_firstLepton
        (
          theBigTree.daughters_px->at (firstDaughterIndex),
          theBigTree.daughters_py->at (firstDaughterIndex),
          theBigTree.daughters_pz->at (firstDaughterIndex),
          theBigTree.daughters_e->at (firstDaughterIndex)
          ) ;
      TLorentzVector tlv_secondLepton
        (
          theBigTree.daughters_px->at (secondDaughterIndex),
          theBigTree.daughters_py->at (secondDaughterIndex),
          theBigTree.daughters_pz->at (secondDaughterIndex),
          theBigTree.daughters_e->at (secondDaughterIndex)
          ) ;
      //if (isDegenere (tlv_firstLepton, tlv_secondLepton)) continue; // antiEle and antiMu should already do this dirty job
      int type1 = theBigTree.particleType->at (firstDaughterIndex) ;
      int type2 = theBigTree.particleType->at (secondDaughterIndex) ;
      int pairType = oph.getPairType (type1, type2) ;
      /*
      //Gen Matching of the selected pair
      bool matched1 = false, matched2=false;
      for(int igen=0; igen<theBigTree.genpart_pdg->size(); igen++){
      if(abs(theBigTree.genpart_pdg->at(igen))!=66615 || pairType != 2) continue;
      TLorentzVector tlv_gen
      (
      theBigTree.genpart_px->at (igen),
      theBigTree.genpart_py->at (igen),
      theBigTree.genpart_pz->at (igen),
      theBigTree.genpart_e->at (igen)
      ) ;
      if(tlv_gen.DeltaR(tlv_firstLepton)<0.5)matched1=true;
      if(tlv_gen.DeltaR(tlv_secondLepton)<0.5)matched2=true;
      }
      if(!(matched1 && matched2)) continue;
      */
      if(isMC){
	if (type1==2) tau_daus.push_back(firstDaughterIndex);
	if (type2==2) tau_daus.push_back(secondDaughterIndex);
      }

      // for tautau pairs, need to check the match to the L1 seed excluded in run D
      if (pairType == 2 && isMC && strategy == 2015)
      {

	bool hasL1Match_1 = theBigTree.daughters_isL1IsoTau28Matched->at (firstDaughterIndex);
	bool hasL1Match_2 = theBigTree.daughters_isL1IsoTau28Matched->at (secondDaughterIndex);
	bool goodL1 = (hasL1Match_1 && hasL1Match_2);
	bool goodHLT = trigReader.checkOR (2, triggerbit);
	if (!goodL1  && skipTriggers < 0) continue;
      }

      if (foundPairs.find (pairType) != foundPairs.end ()) continue ;

      if (skipTriggers < 0)
      {
	//trigger matching
	if (pairType != trigPairType) continue;
	if (!trigReader.checkOR (pairType,theBigTree.daughters_L3FilterFired->at(firstDaughterIndex)))continue;
	if (pairType == 2 && !trigReader.checkOR (pairType,theBigTree.daughters_L3FilterFired->at(secondDaughterIndex)))continue;
      }
      foundPairs[pairType] = iPair ;
    }

    if (foundPairs.size () == 0) continue ;
    // by now take the MOST ISOLATED pair (no OS/SS)
    int chosenTauPair = foundPairs.begin ()->second ;
    int isOS = theBigTree.isOSCand->at (chosenTauPair) ;

    if (saveOS == 1 && !isOS) continue ;
    if (saveOS == 0 &&  isOS) continue ;

    int firstDaughterIndex = theBigTree.indexDau1->at (chosenTauPair) ;
    const TLorentzVector tlv_firstLepton
      (
        theBigTree.daughters_px->at (firstDaughterIndex),
        theBigTree.daughters_py->at (firstDaughterIndex),
        theBigTree.daughters_pz->at (firstDaughterIndex),
        theBigTree.daughters_e->at (firstDaughterIndex)
        ) ;
    int secondDaughterIndex = theBigTree.indexDau2->at (chosenTauPair) ;
    const TLorentzVector tlv_secondLepton
      (
        theBigTree.daughters_px->at (secondDaughterIndex),
        theBigTree.daughters_py->at (secondDaughterIndex),
        theBigTree.daughters_pz->at (secondDaughterIndex),
        theBigTree.daughters_e->at (secondDaughterIndex)
        ) ;

    if(isMC && skipTriggers>0){
      int type1 = theBigTree.particleType->at (firstDaughterIndex) ;
      int type2 = theBigTree.particleType->at (secondDaughterIndex) ;
      int pairType = oph.getPairType (type1, type2) ;

      turnOnReweight = 0;
      std::sort(tau_daus.begin(), tau_daus.end());
      tau_daus.erase(unique(tau_daus.begin(), tau_daus.end()), tau_daus.end());
      float tau_probs[tau_daus.size()];
      int found =0;
      theSmallTree.m_RunNumber = tau_daus.size();

      //now I have a list of unique indexes
      for(int it=0; it<(int)tau_daus.size();it++){
	//cout<<tau_daus.at(it)<<endl;;
	int tau_ind = tau_daus.at(it);
	float x = theBigTree.daughters_px->at(tau_ind);
	float y = theBigTree.daughters_py->at(tau_ind);
	float pt_t = TMath::Sqrt(x*x+y*y);
	float iso1 = getIso (tau_ind, pt_t, theBigTree) ;

	float prob = turnOnSF(pt_t);
	if(skipTriggers==1){
	  tau_probs[it] = prob;
	} else {
	  if(rnd->Rndm()<= prob)found++;
	}
      }
      if(skipTriggers==1){
	if(strategy==2015){
	  for (int bit=0; bit<TMath::Power(2,(int)tau_daus.size());bit++){
	    int nBits = 0;
	    float tmpprob =1.0;
	    for(int it=0; it<(int)tau_daus.size();it++){
	      int tau_ind = tau_daus.at(it);
	      if(oph.checkBit(bit, it)){
		nBits++;
		tmpprob *= tau_probs[it];
	      }else{
		tmpprob *= (1-tau_probs[it]);
	      }
	    }
	    if(nBits>1) turnOnReweight += tmpprob;
	  }

	}else{ //strategy 2016
	  int gen = theBigTree.daughters_genindex->at(firstDaughterIndex);
	  int pdg1=-1,pdg2=-1;
	  if (gen>=0)pdg1=theBigTree.genpart_pdg->at(gen);
	  gen = theBigTree.daughters_genindex->at(secondDaughterIndex);
	  if (gen>=0)pdg2=theBigTree.genpart_pdg->at(gen);
	  if(pairType == 2 )turnOnReweight = turnOnSF(tlv_firstLepton.Pt(),pdg1)*turnOnSF(tlv_secondLepton.Pt(),pdg2);
	  else turnOnReweight = 0.98;//*turnOnSF(tlv_secondLepton.Pt(),theBigTree.daughters_genindex->at(secondDaughterIndex));
	}
      }
      else if(skipTriggers==2) {
	if(strategy == 2016){
	  if(pairType == 2) {
	    if(rnd->Rndm() > turnOnSF(tlv_firstLepton.Pt(),theBigTree.daughters_genindex->at(firstDaughterIndex))) continue;
	    if(rnd->Rndm() > turnOnSF(tlv_secondLepton.Pt(),theBigTree.daughters_genindex->at(secondDaughterIndex))) continue;
	  }else if(rnd->Rndm() > 0.98) continue;
	}
	else if(found<2)continue;
      }
    }

    if (isMC){
      counter.at (selID++) += theBigTree.aMCatNLOweight * reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) * topPtReweight * turnOnReweight;
    } else      counter.at (selID++) += 1 ;

    TLorentzVector tlv_tauH = tlv_firstLepton + tlv_secondLepton ;
    TLorentzVector tlv_tauH_SVFIT ;

    // in case the SVFIT mass is calculated
    if (theBigTree.SVfitMass->at (chosenTauPair) > -900.)
    {
      theSmallTree.m_tauH_SVFIT_mass = theBigTree.SVfitMass->at (chosenTauPair) ;

      theSmallTree.m_tauH_SVFIT_pt = theBigTree.SVfit_pt->at (chosenTauPair) ;
      theSmallTree.m_tauH_SVFIT_eta = theBigTree.SVfit_eta->at (chosenTauPair) ;
      theSmallTree.m_tauH_SVFIT_phi = theBigTree.SVfit_phi->at (chosenTauPair) ;
      theSmallTree.m_tauH_SVFIT_METphi = theBigTree.SVfit_fitMETPhi->at (chosenTauPair) ;
      theSmallTree.m_tauH_SVFIT_METrho = theBigTree.SVfit_fitMETRho->at (chosenTauPair) ;

      tlv_tauH_SVFIT.SetPtEtaPhiM (
	theBigTree.SVfit_pt->at (chosenTauPair),
	theBigTree.SVfit_eta->at (chosenTauPair),
	theBigTree.SVfit_phi->at (chosenTauPair),
	theBigTree.SVfitMass->at (chosenTauPair)
	) ;
    }

    //else cout << "Evt: " << iEvent << " - DEBUG: SVfit for the pair: " << chosenTauPair << " type=" << oph.getPairType (theBigTree.particleType->at (firstDaughterIndex), theBigTree.particleType->at (secondDaughterIndex)) << " is: " << theBigTree.SVfitMass->at (chosenTauPair) << endl;

    vector <pair <int, float> > jets_and_btag ;
    // loop over jets
    for (unsigned int iJet = 0 ; iJet < theBigTree.jets_px->size () ; ++iJet)
    {
      // PG filter jets at will
      if (theBigTree.jets_PUJetID->at (iJet) < PUjetID_minCut) continue ;

      TLorentzVector tlv_jet
	(
	  theBigTree.jets_px->at (iJet),
	  theBigTree.jets_py->at (iJet),
	  theBigTree.jets_pz->at (iJet),
	  theBigTree.jets_e->at (iJet)
	  ) ;
      if (tlv_jet.Pt () < 20. /*GeV*/) continue ;
      if (tlv_jet.DeltaR (tlv_firstLepton) < lepCleaningCone) continue ;
      if (tlv_jet.DeltaR (tlv_secondLepton) < lepCleaningCone) continue ;
      if (TMath::Abs(tlv_jet.Eta()) > 2.4) continue;
      if (tlv_jet.Pt() < 30) continue;

      jets_and_btag.push_back (std::pair <int, float> (
				 iJet, theBigTree.bCSVscore->at (iJet)
				 )) ;

    } // loop over jets

      // apply some selections here
      // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
      // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
      // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
      // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

    if (!beInclusive && jets_and_btag.size () < 2) continue ;
    if (isMC) counter.at (selID++) += theBigTree.aMCatNLOweight * reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) * topPtReweight * turnOnReweight;
    else      counter.at (selID++) += 1 ;

    // fill the variables of interest
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

    const int type1 = theBigTree.particleType->at (firstDaughterIndex) ;
    const int type2 = theBigTree.particleType->at (secondDaughterIndex) ;
    const int pType = oph.getPairType (type1, type2) ;
    theSmallTree.m_pairType = pType ;

    theSmallTree.m_PUReweight = (isMC ? reweight.weight(PUReweight_MC,PUReweight_target,theBigTree.npu) : 1) ;
    theSmallTree.m_MC_weight = (isMC ? theBigTree.aMCatNLOweight * XS * stitchWeight : 1) ;
    vector<float> bTagWeight = bTagSFHelper.getEvtWeight (jets_and_btag, theBigTree.jets_px, theBigTree.jets_py, theBigTree.jets_pz, theBigTree.jets_e, theBigTree.jets_HadronFlavour, pType) ;
    theSmallTree.m_bTagweightL = (isMC ? bTagWeight.at(0) : 1.0) ;
    theSmallTree.m_bTagweightM = (isMC ? bTagWeight.at(1) : 1.0) ;
    theSmallTree.m_bTagweightT = (isMC ? bTagWeight.at(2) : 1.0) ;
    theSmallTree.m_TTtopPtreweight = (isTTBar ? topPtReweight : 1.0) ;
    theSmallTree.m_turnOnreweight = (isMC ? turnOnReweight : 1.0) ;
    theSmallTree.m_lheht = (isMC ? theBigTree.lheHt : 0) ;
    theSmallTree.m_EventNumber = theBigTree.EventNumber ;
    theSmallTree.m_npv = theBigTree.npv ;
    theSmallTree.m_npu = theBigTree.npu ;
    theSmallTree.m_lumi = theBigTree.lumi ;
    theSmallTree.m_triggerbit = theBigTree.triggerbit ;
    theSmallTree.m_rho = theBigTree.rho ;
    theSmallTree.m_isMC = isMC ;
    theSmallTree.m_isOS = theBigTree.isOSCand->at (chosenTauPair) ;
    theSmallTree.m_met_phi = theBigTree.metphi ;
    theSmallTree.m_met_et = theBigTree.met ;
    theSmallTree.m_mT = theBigTree.mT_Dau1->at (chosenTauPair) ;



    theSmallTree.m_tauH_pt = tlv_tauH.Pt () ;
    theSmallTree.m_tauH_eta = tlv_tauH.Eta () ;
    theSmallTree.m_tauH_phi = tlv_tauH.Phi () ;
    theSmallTree.m_tauH_e = tlv_tauH.E () ;
    theSmallTree.m_tauH_mass = tlv_tauH.M () ;

    theSmallTree.m_L3filter1 = theBigTree.daughters_L3FilterFiredLast->at (firstDaughterIndex) ;
    theSmallTree.m_L3filterlast1 = theBigTree.daughters_L3FilterFiredLast->at (firstDaughterIndex) ;
    theSmallTree.m_L3filter2 = theBigTree.daughters_L3FilterFiredLast->at (secondDaughterIndex) ;
    theSmallTree.m_L3filterlast2 = theBigTree.daughters_L3FilterFiredLast->at (secondDaughterIndex) ;


    theSmallTree.m_dau1_iso = getIso (firstDaughterIndex, tlv_firstLepton.Pt (), theBigTree) ;
    theSmallTree.m_dau1_photonPtSumOutsideSignalCone = theBigTree.photonPtSumOutsideSignalCone->at (firstDaughterIndex) ;
    for(int i=0;i<hTauIDS->GetNbinsX();i++){
      if(string(hTauIDS->GetXaxis()->GetBinLabel(i+1))==string("byLooseCombinedIsolationDeltaBetaCorr3Hits")){
	theSmallTree.m_dau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (firstDaughterIndex) & (1 << i);
	theSmallTree.m_dau2_byLooseCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (secondDaughterIndex) & (1 << i);
      }
      if(string(hTauIDS->GetXaxis()->GetBinLabel(i+1))==string("byMediumCombinedIsolationDeltaBetaCorr3Hits")){
	theSmallTree.m_dau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (firstDaughterIndex) & (1 << i);
	theSmallTree.m_dau2_byMediumCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (secondDaughterIndex) & (1 << i);
      }
      if(string(hTauIDS->GetXaxis()->GetBinLabel(i+1))==string("byTightCombinedIsolationDeltaBetaCorr3Hits")){
	theSmallTree.m_dau1_byTightCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (firstDaughterIndex) & (1 << i);
	theSmallTree.m_dau2_byTightCombinedIsolationDeltaBetaCorr3Hits = theBigTree.tauID->at (secondDaughterIndex) & (1 << i);
      }
    }
    theSmallTree.m_dau1_pt = tlv_firstLepton.Pt () ;
    theSmallTree.m_dau1_eta = tlv_firstLepton.Eta () ;
    theSmallTree.m_dau1_phi = tlv_firstLepton.Phi () ;
    theSmallTree.m_dau1_e = theBigTree.daughters_e->at (firstDaughterIndex) ;
    theSmallTree.m_dau1_flav = theBigTree.daughters_charge->at (firstDaughterIndex) *
      (theBigTree.particleType->at (firstDaughterIndex) + 1) ;
    // 1 = from muons collection
    // 2 = from electrons collection
    // 3 = from tauH collection

    theSmallTree.m_dau2_iso = getIso (secondDaughterIndex, tlv_secondLepton.Pt (), theBigTree) ;
    theSmallTree.m_dau2_photonPtSumOutsideSignalCone = theBigTree.photonPtSumOutsideSignalCone->at (secondDaughterIndex) ;
    theSmallTree.m_dau2_pt = tlv_secondLepton.Pt () ;
    theSmallTree.m_dau2_eta = tlv_secondLepton.Eta () ;
    theSmallTree.m_dau2_phi = tlv_secondLepton.Phi () ;
    theSmallTree.m_dau2_e = theBigTree.daughters_e->at (secondDaughterIndex) ;
    theSmallTree.m_dau2_flav = theBigTree.daughters_charge->at (secondDaughterIndex) *
      (theBigTree.particleType->at (secondDaughterIndex) + 1) ;

    float trigSF = 1.0;
    float idAndIsoSF = 1.0;
    // particle 2 is always a TAU --  FIXME: not good for emu
    if (type1 < 2 && isMC) // mu
    {
      trigSF = myScaleFactor[type1][0]->get_ScaleFactor(tlv_firstLepton.Pt(),tlv_firstLepton.Eta());
      idAndIsoSF = myScaleFactor[type1][1]->get_ScaleFactor(tlv_firstLepton.Pt(),tlv_firstLepton.Eta());
    }

    theSmallTree.m_trigSF     = (isMC ? trigSF : 1.0);
    theSmallTree.m_IdAndIsoSF = (isMC ? idAndIsoSF : 1.0);

    // loop over leptons
    vector<pair<TLorentzVector, float> > dummyLeptCollection ;
    for (unsigned int iLep = 0 ;
	 (iLep < theBigTree.daughters_px->size ()) ;
	 ++iLep)
    {
      // skip the H decay candiates
      if (int (iLep) == firstDaughterIndex ||
	  int (iLep) == secondDaughterIndex) continue ;

      // remove taus
      if (theBigTree.particleType->at (iLep) == 2)
      {
	continue ;
      }
      else if (theBigTree.particleType->at (iLep) == 0) // muons
      {
	if (!oph.muBaseline (&theBigTree, iLep, 10., 2.4, 0.3)) continue ;
      }
      else if (theBigTree.particleType->at (iLep) == 1) // electrons
      {
	if (!oph.eleBaseline (&theBigTree, iLep, 10., 0.3, 1)) continue ;
//              if (!oph.eleBaseline (&theBigTree, iLep, 10., 2.5, 0.3, 1)) continue ;
      }
      TLorentzVector tlv_dummyLepton
	(
	  theBigTree.daughters_px->at (iLep),
	  theBigTree.daughters_py->at (iLep),
	  theBigTree.daughters_pz->at (iLep),
	  theBigTree.daughters_e->at (iLep)
	  ) ;
      dummyLeptCollection.push_back (pair<TLorentzVector, float> (
				       tlv_dummyLepton,
				       theBigTree.daughters_charge->at (iLep) *
				       (theBigTree.particleType->at (iLep) + 1)
				       )) ;

    } // loop over leptons

    sort (dummyLeptCollection.rbegin (), dummyLeptCollection.rend (), leptSort ()) ;
    for (unsigned int iLep = 0 ;
	 (iLep < dummyLeptCollection.size ()) && (theSmallTree.m_nleps < 2) ;
	 ++iLep)
    {
      theSmallTree.m_leps_pt.push_back (dummyLeptCollection.at (iLep).first.Pt ()) ;
      theSmallTree.m_leps_eta.push_back (dummyLeptCollection.at (iLep).first.Eta ()) ;
      theSmallTree.m_leps_phi.push_back (dummyLeptCollection.at (iLep).first.Phi ()) ;
      theSmallTree.m_leps_e.push_back (dummyLeptCollection.at (iLep).first.E ()) ;
      theSmallTree.m_leps_flav.push_back (dummyLeptCollection.at (iLep).second) ;
      ++theSmallTree.m_nleps ;
    }

    if (jets_and_btag.size () >= 2)
    {
      pair<int, int> eventJets ;
      if (bChoiceFlag == 1)       eventJets = chooseHighestBtagJets (jets_and_btag) ;
      else if (bChoiceFlag == 2)  eventJets = chooseHighestPtJets (jets_and_btag) ;

      TLorentzVector tlv_firstBjet
	(
	  theBigTree.jets_px->at (eventJets.first),
	  theBigTree.jets_py->at (eventJets.first),
	  theBigTree.jets_pz->at (eventJets.first),
	  theBigTree.jets_e->at (eventJets.first)
	  ) ;
      TLorentzVector tlv_secondBjet
	(
	  theBigTree.jets_px->at (eventJets.second),
	  theBigTree.jets_py->at (eventJets.second),
	  theBigTree.jets_pz->at (eventJets.second),
	  theBigTree.jets_e->at (eventJets.second)
	  ) ;

      // compute b jet energy regression
      // eventJets.first , eventJets.second

      double ptRegr[2] = {tlv_firstBjet.Pt(), tlv_secondBjet.Pt()};
      if (computeBregr)
      {
	for (int iBJet = 0; iBJet <=1; iBJet++)
	{
	  int bidx = (iBJet == 0 ? eventJets.first : eventJets.second);
	  bjrv.Jet_pt     = (iBJet == 0 ? tlv_firstBjet.Pt()  : tlv_secondBjet.Pt());
	  bjrv.Jet_eta    = (iBJet == 0 ? tlv_firstBjet.Eta() : tlv_secondBjet.Eta());
	  //bjrv.Jet_corr         = theBigTree.jets_rawPt->at(bidx);
	  bjrv.Jet_corr   = theBigTree.jetRawf->at(bidx); // should be 1./jetrawf ??
	  bjrv.rho              = theBigTree.rho;
	  bjrv.Jet_mt           = theBigTree.jets_mT->at(bidx);
	  bjrv.Jet_leadTrackPt  = theBigTree.jets_leadTrackPt->at(bidx);
	  bjrv.Jet_leptonPtRel  = theBigTree.jets_leptonPtRel->at(bidx);
	  bjrv.Jet_leptonPt     = theBigTree.jets_leptonPt->at(bidx);
	  bjrv.Jet_leptonDeltaR = theBigTree.jets_leptonDeltaR->at(bidx);
	  bjrv.Jet_neHEF   = theBigTree.jets_nHEF->at(bidx);
	  bjrv.Jet_neEmEF  = theBigTree.jets_nEmEF->at(bidx);
	  bjrv.Jet_chMult  = theBigTree.jets_chMult->at(bidx);
	  bjrv.Jet_vtxPt   = theBigTree.jets_vtxPt->at(bidx);
	  bjrv.Jet_vtxMass = theBigTree.jets_vtxMass->at(bidx);
	  bjrv.Jet_vtx3dL  = theBigTree.jets_vtx3dL->at(bidx);
	  bjrv.Jet_vtxNtrk = theBigTree.jets_vtxNtrk->at(bidx);
	  bjrv.Jet_vtx3deL = theBigTree.jets_vtx3deL->at(bidx);

	  ptRegr[iBJet] = (bRreader->EvaluateRegression (bRegrMethodName.c_str()))[0];
	}
      }

      TLorentzVector tlv_firstBjet_raw = tlv_firstBjet;
      TLorentzVector tlv_secondBjet_raw = tlv_secondBjet;

      // SET raw pt 1 and 2 ordered in pt for a test of b jet regression
      // float FIXME1 = tlv_firstBjet_raw.Pt() ;
      // float FIXME2 = tlv_secondBjet_raw.Pt() ;
      // float theMax, theMin;
      // if (FIXME1 > FIXME2)
      // {
      //   theMax = FIXME1;
      //   theMin = FIXME2;
      // }
      // else
      // {
      //   theMax = FIXME2;
      //   theMin = FIXME1;
      // }

      // theSmallTree.m_bjet1_pt_raw = theMax;
      // theSmallTree.m_bjet2_pt_raw = theMin;

      theSmallTree.m_bjet1_pt_raw = tlv_firstBjet_raw.Pt();
      theSmallTree.m_bjet2_pt_raw = tlv_secondBjet_raw.Pt();

      TLorentzVector tlv_bH_raw = tlv_firstBjet + tlv_secondBjet ;
      theSmallTree.m_bH_mass_raw = tlv_bH_raw.M();
      // FIXME : here mass is manually set to 0, should we change it?
      float ptScale1 = ptRegr[0] / tlv_firstBjet.Pt() ;
      float ptScale2 = ptRegr[1] / tlv_secondBjet.Pt() ;
      tlv_firstBjet.SetPtEtaPhiE (ptRegr[0], tlv_firstBjet.Eta(), tlv_firstBjet.Phi(), ptScale1*tlv_firstBjet.Energy());
      tlv_secondBjet.SetPtEtaPhiE (ptRegr[1], tlv_secondBjet.Eta(), tlv_secondBjet.Phi(), ptScale2*tlv_secondBjet.Energy());

      theSmallTree.m_bjet1_pt  = tlv_firstBjet.Pt () ;
      theSmallTree.m_bjet1_eta = tlv_firstBjet.Eta () ;
      theSmallTree.m_bjet1_phi = tlv_firstBjet.Phi () ;
      theSmallTree.m_bjet1_e = theBigTree.jets_e->at (eventJets.first) ;
      theSmallTree.m_bjet1_bID = theBigTree.bCSVscore->at (eventJets.first) ;
      theSmallTree.m_bjet1_flav = theBigTree.jets_HadronFlavour->at (eventJets.first) ;

      theSmallTree.m_bjet2_pt  = tlv_secondBjet.Pt () ;
      theSmallTree.m_bjet2_eta = tlv_secondBjet.Eta () ;
      theSmallTree.m_bjet2_phi = tlv_secondBjet.Phi () ;
      theSmallTree.m_bjet2_e = theBigTree.jets_e->at (eventJets.second) ;
      theSmallTree.m_bjet2_bID = theBigTree.bCSVscore->at (eventJets.second) ;
      theSmallTree.m_bjet2_flav = theBigTree.jets_HadronFlavour->at (eventJets.second) ;

      bool hasgj1 = false;
      bool hasgj2 = false;
      if (isMC){

	int mcind = theBigTree.jets_genjetIndex->at(eventJets.first);
	if (mcind>=0){
	  TLorentzVector gen(theBigTree.genjet_px->at(mcind),theBigTree.genjet_py->at(mcind),theBigTree.genjet_pz->at(mcind),theBigTree.genjet_e->at(mcind));
	  theSmallTree.m_genjet1_pt = gen.Pt();
	  theSmallTree.m_genjet1_eta = gen.Eta();
	  theSmallTree.m_genjet1_phi = gen.Phi();
	  theSmallTree.m_genjet1_e = gen.E();
	  if (gen.Pt() > 8) hasgj1 = true;
	}

	mcind = theBigTree.jets_genjetIndex->at(eventJets.second);
	if (mcind>=0){
	  TLorentzVector gen(theBigTree.genjet_px->at(mcind),theBigTree.genjet_py->at(mcind),theBigTree.genjet_pz->at(mcind),theBigTree.genjet_e->at(mcind));
	  theSmallTree.m_genjet2_pt = gen.Pt();
	  theSmallTree.m_genjet2_eta = gen.Eta();
	  theSmallTree.m_genjet2_phi = gen.Phi();
	  theSmallTree.m_genjet2_e = gen.E();
	  if (gen.Pt() > 8) hasgj2 = true;
	}
      }
      theSmallTree.m_bjet1_hasgenjet = hasgj1 ;
      theSmallTree.m_bjet2_hasgenjet = hasgj2 ;

      float METx = theBigTree.METx->at (chosenTauPair) ;
      float METy = theBigTree.METy->at (chosenTauPair) ;
      //float METpt = 0;//TMath::Sqrt (METx*METx + METy*METy) ;

      TLorentzVector tlv_bH = tlv_firstBjet + tlv_secondBjet ;
      TLorentzVector tlv_neutrinos =  tlv_bH - tlv_bH_raw;
      theSmallTree.m_met_et_corr = theBigTree.met - tlv_neutrinos.Et() ;

      const TVector2 ptmiss = TVector2(METx, METy) ;
      //TVector2 ptmiss = TVector2(METx,METy);
      TMatrixD metcov (2, 2) ;
      metcov (0,0) = theBigTree.MET_cov00->at (chosenTauPair) ;
      metcov (1,0) = theBigTree.MET_cov10->at (chosenTauPair) ;
      metcov (0,1) = theBigTree.MET_cov01->at (chosenTauPair) ;
      metcov (1,1) = theBigTree.MET_cov11->at (chosenTauPair) ;

      const TMatrixD stableMetCov = metcov;
      theSmallTree.m_bH_pt = tlv_bH.Pt () ;
      theSmallTree.m_bH_eta = tlv_bH.Eta () ;
      theSmallTree.m_bH_phi = tlv_bH.Phi () ;
      theSmallTree.m_bH_e = tlv_bH.E () ;
      theSmallTree.m_bH_mass = tlv_bH.M () ;

      TLorentzVector tlv_HH = tlv_bH + tlv_tauH ;
      TLorentzVector tlv_HH_raw = tlv_bH_raw + tlv_tauH ;
      theSmallTree.m_HH_pt = tlv_HH.Pt () ;
      theSmallTree.m_HH_eta = tlv_HH.Eta () ;
      theSmallTree.m_HH_phi = tlv_HH.Phi () ;
      theSmallTree.m_HH_e = tlv_HH.E () ;
      theSmallTree.m_HH_mass = tlv_HH.M () ;
      theSmallTree.m_HH_mass_raw = tlv_HH_raw.M () ;
      // in case the SVFIT mass is calculated
      if (theBigTree.SVfitMass->at (chosenTauPair) > -900.)
      {
	TLorentzVector tlv_HHsvfit = tlv_bH + tlv_tauH_SVFIT ;
	theSmallTree.m_HHsvfit_pt = tlv_HHsvfit.Pt () ;
	theSmallTree.m_HHsvfit_eta = tlv_HHsvfit.Eta () ;
	theSmallTree.m_HHsvfit_phi = tlv_HHsvfit.Phi () ;
	theSmallTree.m_HHsvfit_e = tlv_HHsvfit.E () ;
	theSmallTree.m_HHsvfit_mass = tlv_HHsvfit.M () ;
      } // in case the SVFIT mass is calculated

//           //intance of fitter master class
//           HHKinFitMaster kinFits = HHKinFitMaster (&tlv_firstBjet, &tlv_secondBjet,

      bool wrongHHK=false;
      float HHKmass = -999;
      float HHKChi2 = -999;
      if (runHHKinFit)
      {
	HHKinFit2::HHKinFitMasterHeavyHiggs kinFits = HHKinFit2::HHKinFitMasterHeavyHiggs ( tlv_firstBjet, tlv_secondBjet,
											    tlv_firstLepton, tlv_secondLepton,  ptmiss, stableMetCov) ;
	HHKinFit2::HHKinFitMasterHeavyHiggs kinFitsraw = HHKinFit2::HHKinFitMasterHeavyHiggs ( tlv_firstBjet_raw, tlv_secondBjet_raw,
											       tlv_firstLepton, tlv_secondLepton,  ptmiss, stableMetCov) ;
	//           kinFits.setAdvancedBalance (&ptmiss, metcov) ;
	//           kinFits.setSimpleBalance (ptmiss.Pt (),10) ; //alternative which uses only the absolute value of ptmiss in the fit
	//
	//           kinFits.addMh1Hypothesis (hypo_mh1) ;
	//           kinFits.addMh2Hypothesis (hypo_mh2) ;
	kinFits.addHypo(hypo_mh1,hypo_mh2);
	kinFitsraw.addHypo(hypo_mh1,hypo_mh2);
	try{
	  kinFits.fit();//doFit () ;
	}
	catch(HHKinFit2::HHInvMConstraintException e){
	  cout<<"INVME THIS EVENT WAS WRONG, INV MASS CONSTRAIN EXCEPTION"<<endl;
	  cout<<"INVME masshypo1 = 125,    masshypo2 = 125"<<endl;
	  cout<<"INVME Tau1"<<endl;
	  cout<<"INVME (E,Px,Py,Pz,M) "<<tlv_firstLepton.E()<<","<<tlv_firstLepton.Px()<<","<<tlv_firstLepton.Py()<<","<<tlv_firstLepton.Pz()<<","<<tlv_firstLepton.M()<<endl;//tlv_firstLepton.Print();
	  cout<<"INVME Tau2"<<endl;
	  cout<<"INVME (E,Px,Py,Pz,M) "<<tlv_secondLepton.E()<<","<<tlv_secondLepton.Px()<<","<<tlv_secondLepton.Py()<<","<<tlv_secondLepton.Pz()<<","<<tlv_secondLepton.M()<<endl;
	  cout<<"INVME B1"<<endl;
	  cout<<"INVME (E,Px,Py,Pz,M) "<<tlv_firstBjet.E()<<","<<tlv_firstBjet.Px()<<","<<tlv_firstBjet.Py()<<","<<tlv_firstBjet.Pz()<<","<<tlv_firstBjet.M()<<endl;
	  cout<<"INVME B2"<<endl;
	  cout<<"INVME (E,Px,Py,Pz,M) "<<tlv_secondBjet.E()<<","<<tlv_secondBjet.Px()<<","<<tlv_secondBjet.Py()<<","<<tlv_secondBjet.Pz()<<","<<tlv_secondBjet.M()<<endl;
	  cout<<"INVME MET"<<endl;
	  cout<<"INVME (E,Px,Py,Pz,M) "<<","<<ptmiss.Px()<<","<<ptmiss.Py()<<endl;
	  cout<<"INVME METCOV "<<endl;
	  cout<<"INVME "<<metcov (0,0)<<"  "<<metcov (0,1)<<endl;// = theBigTree.MET_cov00->at (chosenTauPair) ;
	  cout<<"INVME "<<metcov (1,0)<<"  "<<metcov (1,1)<<endl;// = theBigTree.MET_cov10->at (chosenTauPair) ;
	  cout<<"INVME tau1, tau2, b1, b2"<<endl;
	  cout<<"INVME ";
	  tlv_firstLepton.Print();
	  cout<<"INVME ";
	  tlv_secondLepton.Print();
	  cout<<"INVME ";
	  tlv_firstBjet.Print();
	  cout<<"INVME ";
	  tlv_secondBjet.Print();
	  wrongHHK=true;
	}
	catch (HHKinFit2::HHEnergyRangeException e){
	  cout<<"ERANGE THIS EVENT WAS WRONG, ENERGY RANGE EXCEPTION"<<endl;
	  cout<<"ERANGE masshypo1 = 125,    masshypo2 = 125"<<endl;
	  cout<<"ERANGE Tau1"<<endl;
	  cout<<"ERANGE (E,Px,Py,Pz,M) "<<tlv_firstLepton.E()<<","<<tlv_firstLepton.Px()<<","<<tlv_firstLepton.Py()<<","<<tlv_firstLepton.Pz()<<","<<tlv_firstLepton.M()<<endl;//tlv_firstLepton.Print();
	  cout<<"ERANGE Tau2"<<endl;
	  cout<<"ERANGE (E,Px,Py,Pz,M) "<<tlv_secondLepton.E()<<","<<tlv_secondLepton.Px()<<","<<tlv_secondLepton.Py()<<","<<tlv_secondLepton.Pz()<<","<<tlv_secondLepton.M()<<endl;
	  cout<<"ERANGE B1"<<endl;
	  cout<<"ERANGE (E,Px,Py,Pz,M) "<<tlv_firstBjet.E()<<","<<tlv_firstBjet.Px()<<","<<tlv_firstBjet.Py()<<","<<tlv_firstBjet.Pz()<<","<<tlv_firstBjet.M()<<endl;
	  cout<<"ERANGE B2"<<endl;
	  cout<<"ERANGE (E,Px,Py,Pz,M) "<<tlv_secondBjet.E()<<","<<tlv_secondBjet.Px()<<","<<tlv_secondBjet.Py()<<","<<tlv_secondBjet.Pz()<<","<<tlv_secondBjet.M()<<endl;
	  cout<<"ERANGE MET"<<endl;
	  cout<<"ERANGE (E,Px,Py,Pz,M) "<<","<<ptmiss.Px()<<","<<ptmiss.Py()<<endl;
	  cout<<"ERANGE METCOV "<<endl;
	  cout<<"ERANGE "<<metcov (0,0)<<"  "<<metcov (0,1)<<endl;// = theBigTree.MET_cov00->at (chosenTauPair) ;
	  cout<<"ERANGE "<<metcov (1,0)<<"  "<<metcov (1,1)<<endl;// = theBigTree.MET_cov10->at (chosenTauPair) ;
	  cout<<"ERANGE tau1, tau2, b1, b2"<<endl;
	  cout<<"ERANGE ";
	  tlv_firstLepton.Print();
	  cout<<"ERANGE ";
	  tlv_secondLepton.Print();
	  cout<<"ERANGE ";
	  tlv_firstBjet.Print();
	  cout<<"ERANGE ";
	  tlv_secondBjet.Print();
	  wrongHHK=true;
	}
	catch(HHKinFit2::HHEnergyConstraintException e){
	  cout<<"ECON THIS EVENT WAS WRONG, ENERGY CONSTRAIN EXCEPTION"<<endl;
	  cout<<"ECON masshypo1 = 125,    masshypo2 = 125"<<endl;
	  cout<<"ECON Tau1"<<endl;
	  cout<<"ECON (E,Px,Py,Pz,M) "<<tlv_firstLepton.E()<<","<<tlv_firstLepton.Px()<<","<<tlv_firstLepton.Py()<<","<<tlv_firstLepton.Pz()<<","<<tlv_firstLepton.M()<<endl;//tlv_firstLepton.Print();
	  cout<<"ECON Tau2"<<endl;
	  cout<<"ECON (E,Px,Py,Pz,M) "<<tlv_secondLepton.E()<<","<<tlv_secondLepton.Px()<<","<<tlv_secondLepton.Py()<<","<<tlv_secondLepton.Pz()<<","<<tlv_secondLepton.M()<<endl;
	  cout<<"ECON B1"<<endl;
	  cout<<"ECON (E,Px,Py,Pz,M) "<<tlv_firstBjet.E()<<","<<tlv_firstBjet.Px()<<","<<tlv_firstBjet.Py()<<","<<tlv_firstBjet.Pz()<<","<<tlv_firstBjet.M()<<endl;
	  cout<<"ECON B2"<<endl;
	  cout<<"ECON (E,Px,Py,Pz,M) "<<tlv_secondBjet.E()<<","<<tlv_secondBjet.Px()<<","<<tlv_secondBjet.Py()<<","<<tlv_secondBjet.Pz()<<","<<tlv_secondBjet.M()<<endl;
	  cout<<"ECON MET"<<endl;
	  cout<<"ECON (E,Px,Py,Pz,M) "<<","<<ptmiss.Px()<<","<<ptmiss.Py()<<endl;
	  cout<<"ECON METCOV "<<endl;
	  cout<<"ECON "<<metcov (0,0)<<"  "<<metcov (0,1)<<endl;// = theBigTree.MET_cov00->at (chosenTauPair) ;
	  cout<<"ECON "<<metcov (1,0)<<"  "<<metcov (1,1)<<endl;// = theBigTree.MET_cov10->at (chosenTauPair) ;
	  cout<<"ECON tau1, tau2, b1, b2"<<endl;
	  cout<<"ECON ";
	  tlv_firstLepton.Print();
	  cout<<"ECON ";
	  tlv_secondLepton.Print();
	  cout<<"ECON ";
	  tlv_firstBjet.Print();
	  cout<<"ECON ";
	  tlv_secondBjet.Print();
	  wrongHHK=true;
	}
	if(!wrongHHK){
	  HHKmass = kinFits.getMH () ;
	  HHKChi2 = kinFits.getChi2 () ;
	}else{
	  if(isOS)HHKmass = -333;
	}
	bool wrongHHKraw =false;
	try {
	  kinFitsraw.fit();
	}
	catch(HHKinFit2::HHInvMConstraintException e){wrongHHKraw=true;}
	catch(HHKinFit2::HHEnergyConstraintException e){wrongHHKraw=true;}
	catch (HHKinFit2::HHEnergyRangeException e){wrongHHKraw=true;}
	if(!wrongHHKraw){
	  theSmallTree.m_HHKin_mass_raw = kinFitsraw.getMH();
	}
	else theSmallTree.m_HHKin_mass_raw = -100 ;
	if (theBigTree.SVfitMass->at (chosenTauPair) > -900. && !wrongHHK)
	{
	  TLorentzVector b1 = kinFits.getFittedBJet1();
	  TLorentzVector b2 = kinFits.getFittedBJet2();
	  TLorentzVector bH_HKin = b1 + b2;
	  TLorentzVector tlv_HHsvfit = bH_HKin + tlv_tauH_SVFIT ;

	  theSmallTree.m_HHkinsvfit_bHmass = bH_HKin.M();
	  theSmallTree.m_HHkinsvfit_pt = tlv_HHsvfit.Pt () ;
	  theSmallTree.m_HHkinsvfit_eta = tlv_HHsvfit.Eta () ;
	  theSmallTree.m_HHkinsvfit_phi = tlv_HHsvfit.Phi () ;
	  theSmallTree.m_HHkinsvfit_e = tlv_HHsvfit.E () ;
	  theSmallTree.m_HHkinsvfit_m = tlv_HHsvfit.M () ;
	} // in case the SVFIT mass is calculated

      } // end if doing HHKinFit
      theSmallTree.m_HHKin_mass = HHKmass;//kinFits.getMH () ;
      theSmallTree.m_HHKin_chi2 = HHKChi2;//kinFits.getChi2 () ;


      theSmallTree.m_HH_deltaPhi = deltaPhi (tlv_bH.Phi (), tlv_tauH.Phi ()) ;
      theSmallTree.m_tauHMet_deltaPhi = deltaPhi (theBigTree.metphi, tlv_tauH.Phi ()) ;
      theSmallTree.m_bHMet_deltaPhi = deltaPhi (theBigTree.metphi, tlv_bH.Phi ()) ;

      theSmallTree.m_ditau_deltaPhi = deltaPhi (tlv_firstLepton.Phi (), tlv_secondLepton.Phi ()) ;
      theSmallTree.m_dib_deltaPhi = deltaPhi (tlv_firstBjet.Phi (), tlv_secondBjet.Phi ()) ;

      theSmallTree.m_ditau_deltaR = tlv_firstLepton.DeltaR(tlv_secondLepton) ;
      theSmallTree.m_dib_deltaR = tlv_firstBjet.DeltaR(tlv_secondBjet) ;

      vector <float> dRBTau;
      dRBTau.push_back (tlv_firstLepton.DeltaR(tlv_firstBjet));
      dRBTau.push_back (tlv_firstLepton.DeltaR(tlv_secondBjet));
      dRBTau.push_back (tlv_secondLepton.DeltaR(tlv_firstBjet));
      dRBTau.push_back (tlv_secondLepton.DeltaR(tlv_secondBjet));
      theSmallTree.m_btau_deltaRmin = *std::min_element(dRBTau.begin(), dRBTau.end());
      theSmallTree.m_btau_deltaRmax = *std::max_element(dRBTau.begin(), dRBTau.end());

      // loop over jets
      for (unsigned int iJet = 0 ;
	   (iJet < theBigTree.jets_px->size ()) && (theSmallTree.m_njets < maxNjetsSaved) ;
	   ++iJet)
      {
	// PG filter jets at will
	if (theBigTree.jets_PUJetID->at (iJet) < PUjetID_minCut) continue ;

	// skip the H decay candiates
	if (int (iJet) == eventJets.first ){
	  theSmallTree.m_bjet1_jecUnc = theBigTree.jets_jecUnc->at(iJet);
	  continue;
	}else if(int (iJet) == eventJets.second){
	  theSmallTree.m_bjet2_jecUnc = theBigTree.jets_jecUnc->at(iJet);
	  continue ;
	}
	TLorentzVector tlv_dummyJet (
	  theBigTree.jets_px->at (iJet),
	  theBigTree.jets_py->at (iJet),
	  theBigTree.jets_pz->at (iJet),
	  theBigTree.jets_e->at (iJet)
	  ) ;

	// remove jets that overlap with the tau selected in the leg 1 and 2
	if (type1 == 2) {
	  if (tlv_firstLepton.DeltaR(tlv_dummyJet) < 0.5){
	    theSmallTree.m_dau1_jecUnc = theBigTree.jets_jecUnc->at(iJet);
	    continue;
	  }
	}
	if (type2 == 2) {
	  if (tlv_secondLepton.DeltaR(tlv_dummyJet) < 0.5){
	    theSmallTree.m_dau2_jecUnc = theBigTree.jets_jecUnc->at(iJet);
	    continue;
	  }
	}

	// find matching gen jet
	bool hasgj = false;
	if (isMC)
	{
	  int mcind = theBigTree.jets_genjetIndex->at(iJet);
	  if (mcind>=0)
	  {
	    TLorentzVector thisGenJet(theBigTree.genjet_px->at(mcind),theBigTree.genjet_py->at(mcind),theBigTree.genjet_pz->at(mcind),theBigTree.genjet_e->at(mcind));
	    if (thisGenJet.Pt() > 8) hasgj = true;
	  }
	}

	theSmallTree.m_jets_pt.push_back (tlv_dummyJet.Pt ()) ;
	theSmallTree.m_jets_eta.push_back (tlv_dummyJet.Eta ()) ;
	theSmallTree.m_jets_phi.push_back (tlv_dummyJet.Phi ()) ;
	theSmallTree.m_jets_e.push_back (theBigTree.jets_e->at (iJet)) ;
	theSmallTree.m_jets_btag.push_back (theBigTree.bCSVscore->at (iJet)) ;
	theSmallTree.m_jets_flav.push_back (theBigTree.jets_HadronFlavour->at (iJet)) ;
	theSmallTree.m_jets_jecUnc.push_back (theBigTree.jets_jecUnc->at (iJet)) ;
	theSmallTree.m_jets_hasgenjet.push_back (hasgj) ;
	++theSmallTree.m_njets ;
      } // loop over jets

    } // if there's two jets in the event, at least

    if (isMC) selectedEvents += theBigTree.aMCatNLOweight ;
    else selectedEvents += 1 ;
    ++selectedNoWeightsEventsNum ;

    theSmallTree.Fill () ;
  } // loop over events

  cout << "1: " << totalEvents << endl ;
  cout << "2: " << selectedEvents << endl ;
  cout << "3: " << totalNoWeightsEventsNum << endl ;
  cout << "4: " << selectedNoWeightsEventsNum << endl ;
  for (unsigned int i = 0 ; i < counter.size () ; ++i)
    cout << "5 + i: " << counter.at (i) << endl ;

  if (totalEvents != 0) cout << "efficiency = " << selectedEvents / totalEvents << endl ;
  else                  cout << "NO events found\n" ;
  TH1F h_eff ("h_eff", "h_eff", 4 + counter.size (), 0., float (4 + counter.size ())) ;
  h_eff.SetBinContent (1, totalEvents) ;
  h_eff.SetBinContent (2, selectedEvents) ;
  h_eff.SetBinContent (3, totalNoWeightsEventsNum) ;
  h_eff.SetBinContent (4, selectedNoWeightsEventsNum) ;
  h_eff.GetXaxis()->SetBinLabel (1, "totalEvents") ;
  h_eff.GetXaxis()->SetBinLabel (2, "selectedEvents") ;
  h_eff.GetXaxis()->SetBinLabel (3, "totalNoWeightsEventsNum") ;
  h_eff.GetXaxis()->SetBinLabel (4, "selectedNoWeightsEventsNum") ;
  for (unsigned int i = 0 ; i < counter.size () ; ++i){
    h_eff.SetBinContent (5 + i, counter.at (i)) ;
    TString lab; lab.Form("step %d",i);
    h_eff.GetXaxis()->SetBinLabel (5 + i, lab.Data()) ;
  }
  smallFile->cd() ;
  h_eff.Write () ;
  hFakeSF_den->Write();
  hFakeSF_num->Write();
  smallFile->Write () ;
  smallFile->Close () ;

  bool computeMVA = gConfigParser->readBoolOption ("TMVA::computeMVA");

  if (computeMVA)
  {
    bool doMuTau  = gConfigParser->isDefined("TMVA::weightsMuTau");
    bool doETau   = gConfigParser->isDefined("TMVA::weightsETau");
    bool doTauTau = gConfigParser->isDefined("TMVA::weightsTauTau");
    bool doLepTau = gConfigParser->isDefined("TMVA::weightsLepTau");

    string TMVAweightsTauTau  = "";
    string TMVAweightsMuTau   = "";
    string TMVAweightsETau    = "";
    string TMVAweightsLepTau    = "";

    if (doMuTau)  TMVAweightsMuTau  = gConfigParser->readStringOption ("TMVA::weightsMuTau");
    if (doETau)   TMVAweightsETau   = gConfigParser->readStringOption ("TMVA::weightsETau");
    if (doTauTau) TMVAweightsTauTau = gConfigParser->readStringOption ("TMVA::weightsTauTau");
    if (doLepTau) TMVAweightsLepTau = gConfigParser->readStringOption ("TMVA::weightsLepTau");

    bool TMVAspectatorsIn      = gConfigParser->readBoolOption   ("TMVA::spectatorsPresent");
    vector<string> TMVAspectators = gConfigParser->readStringListOption   ("TMVA::spectators");
    vector<string> TMVAvariables  = gConfigParser->readStringListOption   ("TMVA::variables");

    TFile *outFile = TFile::Open(outputFile,"UPDATE");
    TTree *treenew = (TTree*)outFile->Get("HTauTauTree");

    TMVA::Reader * reader = new TMVA::Reader () ;
    Float_t mvatautau,mvamutau, mvaetau, mvaleptau;
    TBranch *mvaBranchmutau;
    TBranch *mvaBranchtautau;
    TBranch *mvaBranchetau;
    TBranch *mvaBranchleptau;

    vector<float> address (TMVAvariables.size () + TMVAspectators.size () * TMVAspectatorsIn, 0.) ;
    for (unsigned int iv = 0 ; iv < TMVAvariables.size () ; ++iv)
    {
      treenew->SetBranchAddress (TMVAvariables.at (iv).c_str (), &(address.at (iv))) ;
      reader->AddVariable (TMVAvariables.at (iv), &(address.at (iv))) ;
    }

    for (unsigned int iv = 0 ; iv < TMVAspectators.size () && TMVAspectatorsIn ; ++iv)
    {
      int addressIndex = iv + TMVAvariables.size () ;
      treenew->SetBranchAddress (TMVAspectators.at (iv).c_str (), &(address.at (addressIndex))) ;
      reader->AddSpectator (TMVAspectators.at (iv), &(address.at (addressIndex))) ;
    }

    //if (treenew->GetListOfBranches ()->FindObject (mvaName.c_str ())) {
    //  treenew->SetBranchAddress ("MuTauKine", &mvamutau, &mvaBranchmutau) ;
    //  treenew->SetBranchAddress ("TauTauKine", &mvatautau, &mvaBranchtautau) ;
    //}
    //else{

    if (doMuTau)  mvaBranchmutau = treenew->Branch ("MuTauKine", &mvamutau, "MuTauKine/F") ;
    if (doETau)   mvaBranchetau = treenew->Branch ("ETauKine", &mvaetau, "ETauKine/F") ;
    if (doTauTau) mvaBranchtautau = treenew->Branch ("TauTauKine", &mvatautau, "TauTauKine/F") ;
    if (doLepTau) mvaBranchleptau = treenew->Branch ("LepTauKine", &mvaleptau, "LepTauKine/F") ;
    //}
    if (doMuTau)   reader->BookMVA ("MuTauKine",  TMVAweightsMuTau.c_str ()) ;
    if (doETau)    reader->BookMVA ("ETauKine",  TMVAweightsETau.c_str ()) ;
    if (doTauTau)  reader->BookMVA ("TauTauKine",  TMVAweightsTauTau.c_str ()) ;
    if (doLepTau)  reader->BookMVA ("LepTauKine",  TMVAweightsLepTau.c_str ()) ;

    int nentries = treenew->GetEntries();
    for(int i=0;i<nentries;i++){
      treenew->GetEntry(i);

      if (doMuTau)   mvamutau= reader->EvaluateMVA ("MuTauKine") ;
      if (doETau)    mvaetau= reader->EvaluateMVA ("ETauKine") ;
      if (doTauTau)  mvatautau= reader->EvaluateMVA ("TauTauKine") ;
      if (doLepTau)  mvaleptau= reader->EvaluateMVA ("LepTauKine") ;

      if (doMuTau)    mvaBranchmutau->Fill();
      if (doETau)     mvaBranchetau->Fill();
      if (doTauTau)   mvaBranchtautau->Fill();
      if (doLepTau)   mvaBranchleptau->Fill();
    }

    outFile->cd () ;
    h_eff.Write () ;
    treenew->Write ("", TObject::kOverwrite) ;

    delete reader;
  }

  cout << "... SKIM finished, exiting." << endl;
  return 0 ;
}



/* how to prepare bigTree:

   prepare the base class with ROOT
   ------------------------------------------------

   root /data_CMS/cms/govoni/test_submit_to_tier3/HiggsTauTauOutput_DY_-1Events_0Skipped_1436202480.82/output_9.root
   TTree * t = (TTree *) _file0->Get ("HTauTauTree/HTauTauTree")
   t->MakeSelector ("bigTree")

   remove the following functions from the .h
   ------------------------------------------------

   virtual void    Init(TTree *tree);
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   move the Init implementation in the class definition
   ------------------------------------------------

   change the following line
   ------------------------------------------------

   bigTree(TTree ETC....  --->
   bigTree (TChain * inputChain) : fChain (inputChain) { Init(fChain) ; }

   remove all the rest of the implementations
   ------------------------------------------------

   delete bigTree.C
   ------------------------------------------------

   remove the following lines to simplify things
   ------------------------------------------------

   #endif

   #ifdef bigTree_cxx

   remove the inheritance:
   ------------------------------------------------

   class bigTree : public TSelector { --->
   class bigTree { --->

   remove the following lines
   ------------------------------------------------

   #include <TSelector.h>

   remove the following lines from Init
   ------------------------------------------------

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   remove the ClassDef call
   ------------------------------------------------

*/
