////////////////////////////////////////////////////////////////
// Macro to make histogram .root files from trees             //
// Selection implemented here, and kinematic regions defined  //
// Run separately for each sample, channel                    //
////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include "TSystem.h"
#include "RooUnfold/src/RooUnfold.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldDagostini.h"
#include "RooUnfold/src/RooUnfoldErrors.h"
#include "RooUnfold/src/RooUnfoldInvert.h"
#include "RooUnfold/src/RooUnfoldParms.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
#include "RooUnfold/src/RooUnfoldTUnfold.h"
#include "BTagCalibrationStandalone.h"

#include <iostream>
#include <string>
#include <vector>

#include "TLorentzVector.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);

float getMuonTrkSF(float eta, bool returnSF){
  // Taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonWorkInProgressAndPagResults#Results_on_the_full_2016_data
  float edges[16] = {-2.4,-2.1,-1.6,-1.2,-0.9,-0.6,-0.3,-0.2,0.2,0.3,0.6,0.9,1.2,1.6,2.1,2.4};
  float SFs[15] = {0.991237,0.994853,0.996413,0.997157,0.997512,0.99756,0.996745,0.996996,0.99772,0.998604,0.998321,0.997682,0.995252,0.994919,0.987334};
  float errUp[15] = {0.000748552,0.000200175,0.000173235,0.000186386,0.0000895796,0.0000866951,0.000163849,0.0000729397,0.000165896,0.0000835331,0.0000921935,0.0001819,0.00017926,0.000184303,0.000729069};
  float errDn[15] = {0.000754067,0.000200328,0.000174266,0.000185342,0.0000897593,0.0000873434,0.000164963,0.0000727173,0.000166809,0.000084481,0.0000922173,0.000180505,0.000179773,0.000184726,0.000736922};
  float SF = 0.0;
  float err = 0.0;
  for (int ii = 0; ii < 15; ii++){
    if (edges[ii] <= eta && edges[ii+1] > eta){
      SF = SFs[ii];
      err = max(errUp[ii],errDn[ii]);
    }
  }
  if (returnSF) return SF;
  else return err / SF;
}

float getMuonSF( double eta, double pt, TH2F* h_muID, TH1F* h_muTrig, TString syst, bool usePost){
  
  // Tracking: eta [-2.4,2.4]
  // ID: x: pt [20,120] y: |eta| [0.0,2.4]
  // Trigger: pt [50,500]
  
  float myPt = pt;
  if (pt > 499.9) myPt = 499.9;
  
  int ibin1 = h_muTrig->FindBin(myPt);
  float SF_muTrig = h_muTrig->GetBinContent(ibin1);
  float err_muTrig = abs(h_muTrig->GetBinError(ibin1) / SF_muTrig);
  
  if (pt > 119.9) myPt = 119.9;
  
  int ibin2X = h_muID->GetXaxis()->FindBin(myPt);
  int ibin2Y = h_muID->GetYaxis()->FindBin(abs(eta));
  int ibin2 = h_muID->GetBin(ibin2X,ibin2Y);
  float SF_muID = h_muID->GetBinContent(ibin2);
  float err_muID = abs(h_muID->GetBinError(ibin2) / SF_muID);
  
  float SF_muReco = getMuonTrkSF(eta,true);
  float err_muReco = getMuonTrkSF(eta,false);
  
  // TODO: add isolation SF
  // For now, take muon iso SF as ~1.0 (https://indico.cern.ch/event/605620/contributions/2441087/attachments/1398025/2132153/VHFMeeting_X53_01.18.17.pdf)
  
  // Total uncertainty is statistical uncertainties from each SF, plus 1% for ID, plus 0.5% each for trigger and isolation
  float err_tot = sqrt(pow(err_muTrig,2) + pow(err_muID,2) + pow(err_muReco,2) + pow(0.01,2) + pow(0.005,2) + pow(0.005,2));
  
  float SF_tot = SF_muReco * SF_muID * SF_muTrig;
  
  if (usePost){
    if (syst == "up") return (SF_tot - 0.078 * err_tot) * (1 + 0.994 * err_tot); //TODO: modify to have correct postfit values
    else if (syst == "down") return (SF_tot - 0.078 * err_tot) * (1 - 0.994 * err_tot);
    else return SF_tot - 0.078 * err_tot;
  }
  else{
    if (syst == "up") return SF_tot * (1 + err_tot);
    else if (syst == "down") return SF_tot * (1 - err_tot);
    else return SF_tot;
  }
};

double getElectronSF(double eta, double pt, TH2F* h_elReco, TH2F* h_elID, TH2F* h_elIso, TH1F* h_elTrig, TString syst, bool usePost){

  //Tracking: x: eta [-2.5,2.5] y: pt [35,500]
  //ID: x: pt [10,200] y: |eta| [0,2.5]
  //Iso: x: pt [10,200] y: |eta| [0,2.5]
  //Trig: pt [50,500]
  
  float myPt = pt;
  if (pt > 499.9) myPt = 499.9;

  int ibinX1 = h_elReco->GetXaxis()->FindBin(eta);
  int ibinY1 = h_elReco->GetYaxis()->FindBin(myPt);
  int ibin1 = h_elReco->GetBin(ibinX1,ibinY1);
  float SF_elReco = h_elReco->GetBinContent(ibin1);
  float err_elReco = abs(h_elReco->GetBinError(ibin1) / SF_elReco);

  int ibin2 = h_elTrig->FindBin(myPt);
  float SF_elTrig = h_elTrig->GetBinContent(ibin2);
  float err_elTrig = abs(h_elTrig->GetBinError(ibin2) / SF_elTrig);
  
  if (myPt > 199.9) myPt = 199.9;  

  int ibinX3 = h_elID->GetXaxis()->FindBin(myPt);
  int ibinY3 = h_elID->GetYaxis()->FindBin(abs(eta));
  int ibin3 = h_elID->GetBin(ibinX3,ibinY3);
  float SF_elID = h_elID->GetBinContent(ibin3);
  float err_elID = abs(h_elID->GetBinError(ibin3) / SF_elID);

  int ibinX4 = h_elIso->GetXaxis()->FindBin(myPt);
  int ibinY4 = h_elIso->GetYaxis()->FindBin(abs(eta));
  int ibin4 = h_elIso->GetBin(ibinX4,ibinY4);
  float SF_elIso = h_elIso->GetBinContent(ibin4);
  float err_elIso = abs(h_elIso->GetBinError(ibin4) / SF_elIso);

  float err_tot = sqrt(pow(err_elReco,2)+pow(err_elID,2)+pow(err_elIso,2)+pow(err_elTrig,2));
  
  // Add 1% systematic for high-pt electrons (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors)
  if (myPt > 80.0) err_tot = sqrt(pow(err_tot,2) + pow(0.01,2));
  
  float SF_tot = SF_elReco * SF_elID * SF_elIso * SF_elTrig;

  if (usePost){
    if (syst == "up") return (SF_tot + 1.77 * err_tot) * (1 + 0.47 * err_tot); //TODO: modify to have correct postfit values
    else if (syst == "down") return (SF_tot + 1.77 * err_tot) * (1 - 0.47 * err_tot);
    else return SF_tot + 1.77 * err_tot;
  }
  else{
    if (syst == "up") return SF_tot * (1 + err_tot);
    else if (syst == "down") return SF_tot * (1 - err_tot);
    else return SF_tot;
  }
};

// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------


void makeHists(TString INDIR, TString OUTDIR, TString sample, TString channel, bool isData = false, bool isSignal = false, TString lepID = "Medium", TString iso = "None", bool doHemiCuts = false, float metCut = 0.0, bool doTriangular = false, bool isQCD = false, TString systematic = "nom", int oddOrEven = 0, bool usePost = false) {
  
  SetPlotStyle();
  TH1::AddDirectory(kFALSE);
  
  // minimum top pt cut
  float MINTOPPT = 400.0;

  //--------------------------
  // Setup for response matrix
  //--------------------------
  
  gSystem->Load("RooUnfold/libRooUnfold.so");

  const int nptbins = 8;
  float ptbins[nptbins+1] = {0.0,200.0,400.0,500.0,600.0,700.0,800.0,1200.0,2000.0};
  TH1D* h_bins = new TH1D("bins", ";;", nptbins, ptbins);

  // additional binnings for unfolding checks 
  const int nptbins2 = 10;
  float ptbins2[nptbins2+1] = {0.0,300.0,400.0,450.0,525.0,600.0,700.0,800.0,950.0,1200.0,2000.0};
  const int nptbins3 = 9;
  float ptbins3[nptbins3+1] = {300.0,400.0,450.0,525.0,600.0,700.0,800.0,950.0,1200.0,2000.0};
  const int nptbins4 = 7;
  float ptbins4[nptbins4+1] = {300.0,400.0,500.0,600.0,700.0,800.0,1200.0,2000.0};
  const int nptbins5 = 7;
  float ptbins5[nptbins5+1] = {300.0,400.0,475.0,575.0,700.0,850.0,1200.0,2000.0};
  const int nptbins5fine = 14;
  float ptbins5fine[nptbins5fine+1] = {300.0,350.0,400.0,430.0,475.0,525.0,575.0,625.0,700.0,775.0,850.0,1025.0,1200.0,1600.0,2000.0};
  const int nptbins6 = 7;
  float ptbins6[nptbins6+1] = {300.0,400.0,475.0,550.0,650.0,800.0,1200.0,2000.0};
  const int nptbins6fine = 14;
  float ptbins6fine[nptbins6fine+1] = {300.0,350.0,400.0,430.0,475.0,505.0,550.0,600.0,650.0,725.0,800.0,1000.0,1200.0,1600.0,2000.0};

  TH1D* h_bins2 = new TH1D("bins2", ";;", nptbins2, ptbins2);
  TH1D* h_bins3 = new TH1D("bins3", ";;", nptbins3, ptbins3);
  TH1D* h_bins4 = new TH1D("bins4", ";;", nptbins4, ptbins4);
  TH1D* h_bins5 = new TH1D("bins5", ";;", nptbins5, ptbins5);
  TH1D* h_bins5fine = new TH1D("bins5fine", ";;", nptbins5fine, ptbins5fine);
  TH1D* h_bins6 = new TH1D("bins6", ";;", nptbins6, ptbins6);
  TH1D* h_bins6fine = new TH1D("bins6fine", ";;", nptbins6fine, ptbins6fine);

  RooUnfoldResponse response(h_bins, h_bins);
  response.SetName("response_pt");
  TH1D* h_ptGenTop = new TH1D("ptGenTop", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins, ptbins);
  TH1D* h_ptRecoTop = new TH1D("ptRecoTop", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins, ptbins);

  // additional binnings for unfolding checks 
  RooUnfoldResponse response2(h_bins2, h_bins2);
  RooUnfoldResponse response3(h_bins3, h_bins3);
  RooUnfoldResponse response4(h_bins4, h_bins4);
  RooUnfoldResponse response5(h_bins5, h_bins5);
  RooUnfoldResponse response5fine(h_bins5fine, h_bins5);
  RooUnfoldResponse response6(h_bins6, h_bins6);
  RooUnfoldResponse response6fine(h_bins6fine, h_bins6);
  response2.SetName("response_pt2");
  response3.SetName("response_pt3");
  response4.SetName("response_pt4");
  response5.SetName("response_pt5");
  response5fine.SetName("response_pt5fine");
  response6.SetName("response_pt6");
  response6fine.SetName("response_pt6fine");

  TH1D* h_ptGenTop2 = new TH1D("ptGenTop2", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptGenTop3 = new TH1D("ptGenTop3", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptGenTop4 = new TH1D("ptGenTop4", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptGenTop5 = new TH1D("ptGenTop5", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptGenTop6 = new TH1D("ptGenTop6", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);

  TH1D* h_ptRecoTop2 = new TH1D("ptRecoTop2", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptRecoTop3 = new TH1D("ptRecoTop3", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptRecoTop4 = new TH1D("ptRecoTop4", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptRecoTop5 = new TH1D("ptRecoTop5", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptRecoTop5fine = new TH1D("ptRecoTop5fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5fine, ptbins5fine);
  TH1D* h_ptRecoTop6 = new TH1D("ptRecoTop6", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);
  TH1D* h_ptRecoTop6fine = new TH1D("ptRecoTop6fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6fine, ptbins6fine);

  // reweighted top pt spectra at reco/truth level (using weights based on parton-level top quark pt)
  TH1D* h_ptGenTopMod  = new TH1D("ptGenTopMod",  ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins, ptbins);
  TH1D* h_ptGenTopMod2 = new TH1D("ptGenTopMod2", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptGenTopMod3 = new TH1D("ptGenTopMod3", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptGenTopMod4 = new TH1D("ptGenTopMod4", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptGenTopMod5 = new TH1D("ptGenTopMod5", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptGenTopMod6 = new TH1D("ptGenTopMod6", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);
  TH1D* h_ptRecoTopMod  = new TH1D("ptRecoTopMod",  ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins, ptbins);
  TH1D* h_ptRecoTopMod2 = new TH1D("ptRecoTopMod2", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptRecoTopMod3 = new TH1D("ptRecoTopMod3", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptRecoTopMod4 = new TH1D("ptRecoTopMod4", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptRecoTopMod5 = new TH1D("ptRecoTopMod5", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptRecoTopMod5fine = new TH1D("ptRecoTopMod5fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5fine, ptbins5fine);
  TH1D* h_ptRecoTopMod6 = new TH1D("ptRecoTopMod6", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);
  TH1D* h_ptRecoTopMod6fine = new TH1D("ptRecoTopMod6fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6fine, ptbins6fine);

  TH1D* h_ptGenTopModDown  = new TH1D("ptGenTopModDown",  ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins, ptbins);
  TH1D* h_ptGenTopModDown2 = new TH1D("ptGenTopModDown2", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptGenTopModDown3 = new TH1D("ptGenTopModDown3", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptGenTopModDown4 = new TH1D("ptGenTopModDown4", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptGenTopModDown5 = new TH1D("ptGenTopModDown5", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptGenTopModDown6 = new TH1D("ptGenTopModDown6", ";p_{T}(generated top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);
  TH1D* h_ptRecoTopModDown  = new TH1D("ptRecoTopModDown",  ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins, ptbins);
  TH1D* h_ptRecoTopModDown2 = new TH1D("ptRecoTopModDown2", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins2, ptbins2);
  TH1D* h_ptRecoTopModDown3 = new TH1D("ptRecoTopModDown3", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins3, ptbins3);
  TH1D* h_ptRecoTopModDown4 = new TH1D("ptRecoTopModDown4", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins4, ptbins4);
  TH1D* h_ptRecoTopModDown5 = new TH1D("ptRecoTopModDown5", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5, ptbins5);
  TH1D* h_ptRecoTopModDown5fine = new TH1D("ptRecoTopModDown5fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins5fine, ptbins5fine);
  TH1D* h_ptRecoTopModDown6 = new TH1D("ptRecoTopModDown6", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6, ptbins6);
  TH1D* h_ptRecoTopModDown6fine = new TH1D("ptRecoTopModDown6fine", ";p_{T}(reconstructed top) [GeV]; Events / 10 GeV", nptbins6fine, ptbins6fine);

  float LUM = 35867.0; //From https://hypernews.cern.ch/HyperNews/CMS/get/physics-announcements/4495/1.html -- not filtered for SingleMuon or SingleElectron yet

  double weight_response = LUM * 831.76 / 77229341.; //lum * xsec / Nevents for PowhegPythia8

  // ----------------------------------------------------------------------------------------------------------
  // If running on signal, load truth information for events not passing reco, in order to fill response matrix
  // ----------------------------------------------------------------------------------------------------------

  if (sample.Contains("fullTruth")){
    TChain* treeTO = new TChain("trueTree");
    treeTO->Add(INDIR + sample + ".root");
    
    if (treeTO->GetEntries() == 0) {
      cout << "File doesn't exist or is empty, returning..." << endl;
      return;
    }
    
    vector<int>*   truthChannel_TO = 0;
    vector<float>* genTopPt_TO     = 0;
    vector<float>* genTopEta_TO    = 0;
    vector<float>* genTopPhi_TO    = 0;
    vector<float>* genTTbarMass_TO = 0;
    vector<float>* eventWeight_nom_TO = 0;
    vector<float>* eventWeight_puUp_TO = 0;
    vector<float>* eventWeight_puDown_TO = 0;
    vector<float>* eventWeight_Q2Up_TO = 0;
    vector<float>* eventWeight_Q2Down_TO = 0;
    vector<float>* eventWeight_PDF_TO = 0;
    vector<float>* eventWeight_alphaUp_TO = 0;
    vector<float>* eventWeight_alphaDown_TO = 0;
    
    TBranch* b_truthChannel_TO;
    TBranch* b_genTopPt_TO;
    TBranch* b_genTopEta_TO;
    TBranch* b_genTopPhi_TO;
    TBranch* b_genTTbarMass_TO;
    TBranch* b_eventWeight_nom_TO;
    TBranch* b_eventWeight_puUp_TO;
    TBranch* b_eventWeight_puDown_TO;
    TBranch* b_eventWeight_Q2Up_TO;
    TBranch* b_eventWeight_Q2Down_TO;
    TBranch* b_eventWeight_PDF_TO;
    TBranch* b_eventWeight_alphaUp_TO;
    TBranch* b_eventWeight_alphaDown_TO;
    
    treeTO->SetBranchAddress("truthChannel"          , &truthChannel_TO          , &b_truthChannel_TO          );     
    treeTO->SetBranchAddress("genTopPt"              , &genTopPt_TO              , &b_genTopPt_TO              );
    treeTO->SetBranchAddress("genTopEta"             , &genTopEta_TO             , &b_genTopEta_TO             );
    treeTO->SetBranchAddress("genTopPhi"             , &genTopPhi_TO             , &b_genTopPhi_TO             );
    treeTO->SetBranchAddress("genTTbarMass"          , &genTTbarMass_TO          , &b_genTTbarMass_TO          );
    treeTO->SetBranchAddress("eventWeight_nom"       , &eventWeight_nom_TO       , &b_eventWeight_nom_TO       );
    treeTO->SetBranchAddress("eventWeight_puUp"      , &eventWeight_puUp_TO      , &b_eventWeight_puUp_TO      );
    treeTO->SetBranchAddress("eventWeight_puDown"    , &eventWeight_puDown_TO    , &b_eventWeight_puDown_TO    );
    treeTO->SetBranchAddress("eventWeight_Q2Up"      , &eventWeight_Q2Up_TO      , &b_eventWeight_Q2Up_TO      );
    treeTO->SetBranchAddress("eventWeight_Q2Down"    , &eventWeight_Q2Down_TO    , &b_eventWeight_Q2Down_TO    );
    treeTO->SetBranchAddress("eventWeight_PDF"       , &eventWeight_PDF_TO       , &b_eventWeight_PDF_TO       );
    treeTO->SetBranchAddress("eventWeight_alphaUp"   , &eventWeight_alphaUp_TO   , &b_eventWeight_alphaUp_TO   );
    treeTO->SetBranchAddress("eventWeight_alphaDown" , &eventWeight_alphaDown_TO , &b_eventWeight_alphaDown_TO );
    
    for (int i=0; i<treeTO->GetEntries(); i++) {
      
      treeTO->GetEntry(i,0);
      
      if (i%1000000 == 0) cout << "Event " << i << " (trueTree)" << endl;
      
      if (oddOrEven != 0){
	if (i % 2 == oddOrEven - 1) continue;
      }

      float weight = eventWeight_nom_TO->at(0);

      if (systematic == "puUp")    weight = eventWeight_puUp_TO->at(0);
      if (systematic == "puDown")  weight = eventWeight_puDown_TO->at(0);
      if (systematic == "Q2Up")    weight *= eventWeight_Q2Up_TO->at(0);
      if (systematic == "Q2Down")  weight *= eventWeight_Q2Down_TO->at(0);
      if (systematic == "PDFUp")   weight *= (1.0 + eventWeight_PDF_TO->at(0));
      if (systematic == "PDFDown") weight *= (1.0 - eventWeight_PDF_TO->at(0));
      if (systematic == "ASUp")    weight *= eventWeight_alphaUp_TO->at(0);
      if (systematic == "ASDown")  weight *= eventWeight_alphaDown_TO->at(0);

      // Do channel selection at parton level
      if ((int)truthChannel_TO->size() == 0){
	cout << "Error: no truthChannel information!" << endl;
	continue;
      }
      if (truthChannel_TO->at(0) == 2) {
	cout << "Error: signal event is not semileptonic at parton level!" << endl;
	continue;
      }
      if (truthChannel_TO->at(0) == 0 && channel == "el") continue;
      if (truthChannel_TO->at(0) == 1 && channel == "mu") continue;

      h_ptGenTop->Fill(genTopPt_TO->at(0),weight);
      response.Miss(genTopPt_TO->at(0),weight*weight_response);

      h_ptGenTop2->Fill(genTopPt_TO->at(0),weight);
      h_ptGenTop3->Fill(genTopPt_TO->at(0),weight);
      h_ptGenTop4->Fill(genTopPt_TO->at(0),weight);
      h_ptGenTop5->Fill(genTopPt_TO->at(0),weight);
      h_ptGenTop6->Fill(genTopPt_TO->at(0),weight);
      response2.Miss(genTopPt_TO->at(0),weight*weight_response);
      response3.Miss(genTopPt_TO->at(0),weight*weight_response);
      response4.Miss(genTopPt_TO->at(0),weight*weight_response);
      response5.Miss(genTopPt_TO->at(0),weight*weight_response);
      response5fine.Miss(genTopPt_TO->at(0),weight*weight_response);
      response6.Miss(genTopPt_TO->at(0),weight*weight_response);
      response6fine.Miss(genTopPt_TO->at(0),weight*weight_response);

      float oldpt = genTopPt_TO->at(0);
      float w_ptup = (1.0 + 0.0008*oldpt);
      float w_ptdn = (1.0 - 0.0008*oldpt);

      h_ptGenTopMod->Fill(oldpt,weight*w_ptup);
      h_ptGenTopMod2->Fill(oldpt,weight*w_ptup);
      h_ptGenTopMod3->Fill(oldpt,weight*w_ptup);
      h_ptGenTopMod4->Fill(oldpt,weight*w_ptup);
      h_ptGenTopMod5->Fill(oldpt,weight*w_ptup);
      h_ptGenTopMod6->Fill(oldpt,weight*w_ptup);

      h_ptGenTopModDown->Fill(oldpt,weight*w_ptdn);
      h_ptGenTopModDown2->Fill(oldpt,weight*w_ptdn);
      h_ptGenTopModDown3->Fill(oldpt,weight*w_ptdn);
      h_ptGenTopModDown4->Fill(oldpt,weight*w_ptdn);
      h_ptGenTopModDown5->Fill(oldpt,weight*w_ptdn);
      h_ptGenTopModDown6->Fill(oldpt,weight*w_ptdn);

    }
    treeTO->Delete();
  }
  
  // -------------------------
  // Load information for reco
  // -------------------------
  
  BTagCalibration calib("CSVv2", "CSVv2_Moriond17_B_H.csv");
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, "central", {"up","down"});

  reader.load(calib, BTagEntry::FLAV_B, "comb"); //b
  reader.load(calib, BTagEntry::FLAV_C, "comb"); //c
  reader.load(calib, BTagEntry::FLAV_UDSG, "incl"); //light

  
  // ---------------------------------------------------------------------------------------------------------------
  // Get SFs

  // Muon ID
  // Taken from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults#Results_on_the_full_2016_data
  TFile* f_muID_BCDEF = TFile::Open("EfficienciesAndSF_BCDEF.root","read");
  TFile* f_muID_GH = TFile::Open("EfficienciesAndSF_GH.root","read");
  TH2F* h_muID;

  if (lepID == "Medium") {
    TH2F* h_muID_BCDEF = (TH2F*) f_muID_BCDEF->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
    TH2F* h_muID_GH = (TH2F*) f_muID_GH->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
    h_muID_BCDEF->Sumw2();
    h_muID_GH->Sumw2();
    h_muID_BCDEF->Scale(20.9 / 38.0); // Combine SFs on basis of luminosity: 20.9 (BCDEF), 17.1 (GH)
    h_muID_GH->Scale(17.1 / 38.0);    // TODO: fix to match latest lumi?
    h_muID = (TH2F*) h_muID_BCDEF->Clone();
    h_muID->Add(h_muID_GH);
    h_muID_BCDEF->Delete();
    h_muID_GH->Delete();
  }
  else {
    TH2F* h_muID_BCDEF = (TH2F*) f_muID_BCDEF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
    TH2F* h_muID_GH = (TH2F*) f_muID_GH->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone();
    h_muID_BCDEF->Sumw2();
    h_muID_GH->Sumw2();
    h_muID_BCDEF->Scale(20.9 / 38.0);
    h_muID_GH->Scale(17.1 / 38.0);
    h_muID = (TH2F*) h_muID_BCDEF->Clone();
    h_muID->Add(h_muID_GH);
    h_muID_BCDEF->Delete();
    h_muID_GH->Delete();
  }
  f_muID_BCDEF->Close();
  f_muID_GH->Close();
  delete f_muID_BCDEF;
  delete f_muID_GH;
  if (h_muID == 0) cout << "Muon ID SFs are missing!" << endl;

  // Muon isolation
  // We are currently assuming this to be ~1
  
  // Muon trigger
  TFile* f_muTrig = TFile::Open("muTrigEff.root","read");
  TH1F* h_muTrig = (TH1F*) f_muTrig->Get("lepPt")->Clone();
  f_muTrig->Close();
  delete f_muTrig;
  if (h_muTrig == 0) cout << "Muon trigger SFs are missing!" << endl;

  // Electron reconstruction
  // Taken from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
  TFile* f_elReco = TFile::Open("eleRecoSF.root","read");
  TH2F* h_elReco = (TH2F*) f_elReco->Get("EGamma_SF2D")->Clone();
  f_elReco->Close();
  delete f_elReco;
  if (h_elReco == 0) cout << "Electron reco SFs are missing!" << endl;

  // Electron ID and isolation
  // Taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Scale_Factors_for_Moriond2017
  TFile* f_elSFs = TFile::Open("scaleFactors.root","read"); 
  TH2F* h_elID;
  if (lepID == "Medium") h_elID = (TH2F*) f_elSFs->Get("GsfElectronToCutBasedSpring15M")->Clone();
  else h_elID = (TH2F*) f_elSFs->Get("GsfElectronToCutBasedSpring15T")->Clone();
  TH2F* h_elIso = (TH2F*) f_elSFs->Get("MVAVLooseElectronToMini")->Clone();
  f_elSFs->Close();
  delete f_elSFs;
  if (h_elID == 0) cout << "Electron ID SFs are missing!" << endl;
  if (h_elIso == 0) cout << "Electron Iso SFs are missing!" << endl;

  // Electron trigger
  TFile* f_elTrig = TFile::Open("elTrigEff.root","read");
  TH1F* h_elTrig = (TH1F*) f_elTrig->Get("lepPt")->Clone();
  f_elTrig->Close();
  delete f_elTrig;
  if (h_elTrig == 0) cout << "Electron trigger SFs are missing!" << endl;

  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("recoTree");
  tree->Add(INDIR + sample + ".root");
  
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
  
  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches
  
  // reco level
  vector<float>* metPt                 = 0;
  vector<float>* metPhi                = 0;
  vector<float>* ht                    = 0;
  vector<float>* muPt                  = 0;
  vector<float>* muEta                 = 0;
  vector<float>* muPhi                 = 0;
  vector<float>* muMiniIso             = 0;
  vector<float>* muPtRelPt15           = 0;
  vector<float>* muPtRelPt20           = 0;
  vector<float>* muPtRelPt25           = 0;
  vector<float>* muPtRelPt30           = 0;
  vector<float>* muPtRelPt35           = 0;
  vector<float>* muPtRelPt40           = 0;
  vector<float>* muPtRelPt45           = 0;
  vector<float>* mudRPt15              = 0;
  vector<float>* mudRPt20              = 0;
  vector<float>* mudRPt25              = 0;
  vector<float>* mudRPt30              = 0;
  vector<float>* mudRPt35              = 0;
  vector<float>* mudRPt40              = 0;
  vector<float>* mudRPt45              = 0;
  vector<int>*   muTight               = 0;
  vector<float>* muCharge              = 0;
  vector<float>* elPt                  = 0;
  vector<float>* elEta                 = 0;
  vector<float>* elPhi                 = 0;
  vector<float>* elMiniIso             = 0;
  vector<float>* elPtRelPt15           = 0;
  vector<float>* elPtRelPt20           = 0;
  vector<float>* elPtRelPt25           = 0;
  vector<float>* elPtRelPt30           = 0;
  vector<float>* elPtRelPt35           = 0;
  vector<float>* elPtRelPt40           = 0;
  vector<float>* elPtRelPt45           = 0;
  vector<float>* eldRPt15              = 0;
  vector<float>* eldRPt20              = 0;
  vector<float>* eldRPt25              = 0;
  vector<float>* eldRPt30              = 0;
  vector<float>* eldRPt35              = 0;
  vector<float>* eldRPt40              = 0;
  vector<float>* eldRPt45              = 0;
  vector<int>*   elTight               = 0;
  vector<float>* elCharge              = 0;
  vector<float>* ak4jetPt              = 0;
  vector<float>* ak4jetEta             = 0;
  vector<float>* ak4jetPhi             = 0;
  vector<float>* ak4jetMass            = 0;
  vector<float>* ak4jetCSV             = 0;
  vector<float>* ak4jetVtxMass         = 0;
  vector<float>* ak4jetJECunc          = 0;
  vector<float>* ak4jetPtJERup         = 0;
  vector<float>* ak4jetPtJERdown       = 0;
  vector<float>* ak4jetEtaJERup        = 0;
  vector<float>* ak4jetEtaJERdown      = 0;
  vector<float>* ak4jetPhiJERup        = 0;
  vector<float>* ak4jetPhiJERdown      = 0;
  vector<float>* ak4jetMassJERup       = 0;
  vector<float>* ak4jetMassJERdown     = 0;
  vector<float>* ak4jetHadronFlavour   = 0;
  vector<float>* ak8jetPt              = 0;
  vector<float>* ak8jetEta             = 0;
  vector<float>* ak8jetPhi             = 0;
  vector<float>* ak8jetMass            = 0;
  vector<float>* ak8jetMassPruned      = 0;
  vector<float>* ak8jetMassFiltered    = 0;
  vector<float>* ak8jetMassTrimmed     = 0;   
  vector<float>* ak8jetTau1            = 0;
  vector<float>* ak8jetTau2            = 0;
  vector<float>* ak8jetTau3            = 0;
  vector<float>* ak8jetCSV             = 0;
  vector<float>* ak8jetSDmass          = 0;
  vector<float>* ak8jetSDsubjet0pt     = 0;
  vector<float>* ak8jetSDsubjet0eta    = 0;
  vector<float>* ak8jetSDsubjet0phi    = 0;
  vector<float>* ak8jetSDsubjet0mass   = 0;
  vector<float>* ak8jetSDsubjet0CSV    = 0;
  vector<float>* ak8jetSDsubjet1pt     = 0;
  vector<float>* ak8jetSDsubjet1eta    = 0;
  vector<float>* ak8jetSDsubjet1phi    = 0;
  vector<float>* ak8jetSDsubjet1mass   = 0;
  vector<float>* ak8jetSDsubjet1CSV    = 0;
  vector<float>* ak8jetJECunc          = 0;
  vector<float>* ak8jetPtJERup         = 0;
  vector<float>* ak8jetPtJERdown       = 0;
  vector<float>* ak8jetEtaJERup        = 0;
  vector<float>* ak8jetEtaJERdown      = 0;
  vector<float>* ak8jetPhiJERup        = 0;
  vector<float>* ak8jetPhiJERdown      = 0;
  vector<float>* ak8jetMassJERup       = 0;
  vector<float>* ak8jetMassJERdown     = 0;

  vector<float>* eventWeight_nom       = 0;
  vector<float>* eventWeight_puUp      = 0;
  vector<float>* eventWeight_puDown    = 0;
  vector<float>* eventWeight_Q2Up      = 0;
  vector<float>* eventWeight_Q2Down    = 0;
  vector<float>* eventWeight_PDF       = 0;
  vector<float>* eventWeight_alphaUp   = 0;
  vector<float>* eventWeight_alphaDown = 0;
  vector<int>*   truthChannel          = 0;

  // BRANCHES
  // reco level
  TBranch* b_metPt                 ;
  TBranch* b_metPhi                ;
  TBranch* b_ht                    ;
  TBranch* b_muPt                  ;
  TBranch* b_muEta                 ;
  TBranch* b_muPhi                 ;
  TBranch* b_muMiniIso             ;
  TBranch* b_muPtRelPt15           ;
  TBranch* b_muPtRelPt20           ;
  TBranch* b_muPtRelPt25           ;
  TBranch* b_muPtRelPt30           ;
  TBranch* b_muPtRelPt35           ;
  TBranch* b_muPtRelPt40           ;
  TBranch* b_muPtRelPt45           ;
  TBranch* b_mudRPt15              ;
  TBranch* b_mudRPt20              ;
  TBranch* b_mudRPt25              ;
  TBranch* b_mudRPt30              ;
  TBranch* b_mudRPt35              ;
  TBranch* b_mudRPt40              ;
  TBranch* b_mudRPt45              ;
  TBranch* b_muTight               ;
  TBranch* b_muCharge              ;
  TBranch* b_elPt                  ;
  TBranch* b_elEta                 ;
  TBranch* b_elPhi                 ;
  TBranch* b_elMiniIso             ;
  TBranch* b_elPtRelPt15           ;
  TBranch* b_elPtRelPt20           ;
  TBranch* b_elPtRelPt25           ;
  TBranch* b_elPtRelPt30           ;
  TBranch* b_elPtRelPt35           ;
  TBranch* b_elPtRelPt40           ;
  TBranch* b_elPtRelPt45           ;
  TBranch* b_eldRPt15              ;
  TBranch* b_eldRPt20              ;
  TBranch* b_eldRPt25              ;
  TBranch* b_eldRPt30              ;
  TBranch* b_eldRPt35              ;
  TBranch* b_eldRPt40              ;
  TBranch* b_eldRPt45              ;
  TBranch* b_elTight               ;
  TBranch* b_elCharge              ;
  TBranch* b_ak4jetPt              ;
  TBranch* b_ak4jetEta             ;
  TBranch* b_ak4jetPhi             ;
  TBranch* b_ak4jetMass            ;
  TBranch* b_ak4jetCSV             ;
  TBranch* b_ak4jetVtxMass         ;
  TBranch* b_ak4jetJECunc          ;
  TBranch* b_ak4jetPtJERup         ;
  TBranch* b_ak4jetPtJERdown       ;
  TBranch* b_ak4jetEtaJERup        ;
  TBranch* b_ak4jetEtaJERdown      ;
  TBranch* b_ak4jetPhiJERup        ;
  TBranch* b_ak4jetPhiJERdown      ;
  TBranch* b_ak4jetMassJERup       ;
  TBranch* b_ak4jetMassJERdown     ;
  TBranch* b_ak4jetHadronFlavour   ;
  TBranch* b_ak8jetPt              ;
  TBranch* b_ak8jetEta             ;
  TBranch* b_ak8jetPhi             ;
  TBranch* b_ak8jetMass            ;
  TBranch* b_ak8jetMassPruned      ;
  TBranch* b_ak8jetMassFiltered    ;
  TBranch* b_ak8jetMassTrimmed     ;   
  TBranch* b_ak8jetTau1            ;
  TBranch* b_ak8jetTau2            ;
  TBranch* b_ak8jetTau3            ;
  TBranch* b_ak8jetCSV             ;
  TBranch* b_ak8jetSDmass          ;
  TBranch* b_ak8jetSDsubjet0pt     ;
  TBranch* b_ak8jetSDsubjet0eta    ;
  TBranch* b_ak8jetSDsubjet0phi    ;
  TBranch* b_ak8jetSDsubjet0mass   ;
  TBranch* b_ak8jetSDsubjet0CSV    ;
  TBranch* b_ak8jetSDsubjet1pt     ;
  TBranch* b_ak8jetSDsubjet1eta    ;
  TBranch* b_ak8jetSDsubjet1phi    ;
  TBranch* b_ak8jetSDsubjet1mass   ;
  TBranch* b_ak8jetSDsubjet1CSV    ;
  TBranch* b_ak8jetJECunc          ;
  TBranch* b_ak8jetPtJERup         ;
  TBranch* b_ak8jetPtJERdown       ;
  TBranch* b_ak8jetEtaJERup        ;
  TBranch* b_ak8jetEtaJERdown      ;
  TBranch* b_ak8jetPhiJERup        ;
  TBranch* b_ak8jetPhiJERdown      ;
  TBranch* b_ak8jetMassJERup       ;
  TBranch* b_ak8jetMassJERdown     ;

  TBranch* b_eventWeight_nom       ;
  TBranch* b_eventWeight_puUp      ;
  TBranch* b_eventWeight_puDown    ;
  TBranch* b_eventWeight_Q2Up      ;
  TBranch* b_eventWeight_Q2Down    ;
  TBranch* b_eventWeight_PDF       ;
  TBranch* b_eventWeight_alphaUp   ;
  TBranch* b_eventWeight_alphaDown ;
  TBranch* b_truthChannel          ;

  tree->SetBranchAddress("metPt"                  , &metPt               , &b_metPt               );
  tree->SetBranchAddress("metPhi"                 , &metPhi              , &b_metPhi              );
  tree->SetBranchAddress("ht"                     , &ht                  , &b_ht                  );
  tree->SetBranchAddress("muPt"                   , &muPt                , &b_muPt                );
  tree->SetBranchAddress("muEta"                  , &muEta               , &b_muEta               );
  tree->SetBranchAddress("muPhi"                  , &muPhi               , &b_muPhi               );
  tree->SetBranchAddress("muMiniIso"              , &muMiniIso           , &b_muMiniIso           );
  tree->SetBranchAddress("muPtRelPt15"            , &muPtRelPt15         , &b_muPtRelPt15         );
  tree->SetBranchAddress("muPtRelPt20"            , &muPtRelPt20         , &b_muPtRelPt20         );
  tree->SetBranchAddress("muPtRelPt25"            , &muPtRelPt25         , &b_muPtRelPt25         );
  tree->SetBranchAddress("muPtRelPt30"            , &muPtRelPt30         , &b_muPtRelPt30         );
  tree->SetBranchAddress("muPtRelPt35"            , &muPtRelPt35         , &b_muPtRelPt35         );
  tree->SetBranchAddress("muPtRelPt40"            , &muPtRelPt40         , &b_muPtRelPt40         );
  tree->SetBranchAddress("muPtRelPt45"            , &muPtRelPt45         , &b_muPtRelPt45         );
  tree->SetBranchAddress("mudRPt15"               , &mudRPt15            , &b_mudRPt15            );
  tree->SetBranchAddress("mudRPt20"               , &mudRPt20            , &b_mudRPt20            );
  tree->SetBranchAddress("mudRPt25"               , &mudRPt25            , &b_mudRPt25            );
  tree->SetBranchAddress("mudRPt30"               , &mudRPt30            , &b_mudRPt30            );
  tree->SetBranchAddress("mudRPt35"               , &mudRPt35            , &b_mudRPt35            );
  tree->SetBranchAddress("mudRPt40"               , &mudRPt40            , &b_mudRPt40            );
  tree->SetBranchAddress("mudRPt45"               , &mudRPt45            , &b_mudRPt45            );
  tree->SetBranchAddress("muTight"                , &muTight             , &b_muTight             );
  tree->SetBranchAddress("muCharge"               , &muCharge            , &b_muCharge            );
  tree->SetBranchAddress("elPt"                   , &elPt                , &b_elPt                );
  tree->SetBranchAddress("elEta"                  , &elEta               , &b_elEta               );
  tree->SetBranchAddress("elPhi"                  , &elPhi               , &b_elPhi               );
  tree->SetBranchAddress("elMiniIso"              , &elMiniIso           , &b_elMiniIso           );
  tree->SetBranchAddress("elPtRelPt15"            , &elPtRelPt15         , &b_elPtRelPt15         );
  tree->SetBranchAddress("elPtRelPt20"            , &elPtRelPt20         , &b_elPtRelPt20         );
  tree->SetBranchAddress("elPtRelPt25"            , &elPtRelPt25         , &b_elPtRelPt25         );
  tree->SetBranchAddress("elPtRelPt30"            , &elPtRelPt30         , &b_elPtRelPt30         );
  tree->SetBranchAddress("elPtRelPt35"            , &elPtRelPt35         , &b_elPtRelPt35         );
  tree->SetBranchAddress("elPtRelPt40"            , &elPtRelPt40         , &b_elPtRelPt40         );
  tree->SetBranchAddress("elPtRelPt45"            , &elPtRelPt45         , &b_elPtRelPt45         );
  tree->SetBranchAddress("eldRPt15"               , &eldRPt15            , &b_eldRPt15            );
  tree->SetBranchAddress("eldRPt20"               , &eldRPt20            , &b_eldRPt20            );
  tree->SetBranchAddress("eldRPt25"               , &eldRPt25            , &b_eldRPt25            );
  tree->SetBranchAddress("eldRPt30"               , &eldRPt30            , &b_eldRPt30            );
  tree->SetBranchAddress("eldRPt35"               , &eldRPt35            , &b_eldRPt35            );
  tree->SetBranchAddress("eldRPt40"               , &eldRPt40            , &b_eldRPt40            );
  tree->SetBranchAddress("eldRPt45"               , &eldRPt45            , &b_eldRPt45            );
  tree->SetBranchAddress("elTight"                , &elTight             , &b_elTight             );
  tree->SetBranchAddress("elCharge"               , &elCharge            , &b_elCharge            );
  tree->SetBranchAddress("ak4jetPt"               , &ak4jetPt            , &b_ak4jetPt            );
  tree->SetBranchAddress("ak4jetEta"              , &ak4jetEta           , &b_ak4jetEta           );
  tree->SetBranchAddress("ak4jetPhi"              , &ak4jetPhi           , &b_ak4jetPhi           );
  tree->SetBranchAddress("ak4jetMass"             , &ak4jetMass          , &b_ak4jetMass          );
  tree->SetBranchAddress("ak4jetCSV"              , &ak4jetCSV           , &b_ak4jetCSV           );
  tree->SetBranchAddress("ak4jetVtxMass"          , &ak4jetVtxMass       , &b_ak4jetVtxMass       );
  tree->SetBranchAddress("ak4jetJECunc"           , &ak4jetJECunc        , &b_ak4jetJECunc        );
  tree->SetBranchAddress("ak8jetPt"               , &ak8jetPt            , &b_ak8jetPt            );
  tree->SetBranchAddress("ak8jetEta"              , &ak8jetEta           , &b_ak8jetEta           );
  tree->SetBranchAddress("ak8jetPhi"              , &ak8jetPhi           , &b_ak8jetPhi           );
  tree->SetBranchAddress("ak8jetMass"             , &ak8jetMass          , &b_ak8jetMass          );
  tree->SetBranchAddress("ak8jetMassPruned"       , &ak8jetMassPruned    , &b_ak8jetMassPruned    );
  tree->SetBranchAddress("ak8jetMassFiltered"     , &ak8jetMassFiltered  , &b_ak8jetMassFiltered  );
  tree->SetBranchAddress("ak8jetMassTrimmed"      , &ak8jetMassTrimmed   , &b_ak8jetMassTrimmed   );
  tree->SetBranchAddress("ak8jetTau1"             , &ak8jetTau1          , &b_ak8jetTau1          );
  tree->SetBranchAddress("ak8jetTau2"             , &ak8jetTau2          , &b_ak8jetTau2          );
  tree->SetBranchAddress("ak8jetTau3"             , &ak8jetTau3          , &b_ak8jetTau3          );
  tree->SetBranchAddress("ak8jetCSV"              , &ak8jetCSV           , &b_ak8jetCSV           );
  tree->SetBranchAddress("ak8jetSDmass"           , &ak8jetSDmass        , &b_ak8jetSDmass        );
  tree->SetBranchAddress("ak8jetSDsubjet0pt"      , &ak8jetSDsubjet0pt   , &b_ak8jetSDsubjet0pt   );
  tree->SetBranchAddress("ak8jetSDsubjet0eta"     , &ak8jetSDsubjet0eta  , &b_ak8jetSDsubjet0eta  );
  tree->SetBranchAddress("ak8jetSDsubjet0phi"     , &ak8jetSDsubjet0phi  , &b_ak8jetSDsubjet0phi  );
  tree->SetBranchAddress("ak8jetSDsubjet0mass"    , &ak8jetSDsubjet0mass , &b_ak8jetSDsubjet0mass );
  tree->SetBranchAddress("ak8jetSDsubjet0CSV"     , &ak8jetSDsubjet0CSV  , &b_ak8jetSDsubjet0CSV  );
  tree->SetBranchAddress("ak8jetSDsubjet1pt"      , &ak8jetSDsubjet1pt   , &b_ak8jetSDsubjet1pt   ); 
  tree->SetBranchAddress("ak8jetSDsubjet1eta"     , &ak8jetSDsubjet1eta  , &b_ak8jetSDsubjet1eta  );
  tree->SetBranchAddress("ak8jetSDsubjet1phi"     , &ak8jetSDsubjet1phi  , &b_ak8jetSDsubjet1phi  );
  tree->SetBranchAddress("ak8jetSDsubjet1mass"    , &ak8jetSDsubjet1mass , &b_ak8jetSDsubjet1mass );
  tree->SetBranchAddress("ak8jetSDsubjet1CSV"     , &ak8jetSDsubjet1CSV  , &b_ak8jetSDsubjet1CSV  );
  tree->SetBranchAddress("ak8jetJECunc"           , &ak8jetJECunc        , &b_ak8jetJECunc        );

  tree->SetBranchAddress("eventWeight_nom"            , &eventWeight_nom         , &b_eventWeight_nom         );
  
  if (!isData){
    tree->SetBranchAddress("eventWeight_puUp"       , &eventWeight_puUp      , &b_eventWeight_puUp      );
    tree->SetBranchAddress("eventWeight_puDown"     , &eventWeight_puDown    , &b_eventWeight_puDown    );
    tree->SetBranchAddress("ak4jetPtJERup"          , &ak4jetPtJERup         , &b_ak4jetPtJERup         );
    tree->SetBranchAddress("ak4jetPtJERdown"        , &ak4jetPtJERdown       , &b_ak4jetPtJERdown       );
    tree->SetBranchAddress("ak4jetEtaJERup"         , &ak4jetEtaJERup        , &b_ak4jetEtaJERup        );
    tree->SetBranchAddress("ak4jetEtaJERdown"       , &ak4jetEtaJERdown      , &b_ak4jetEtaJERdown      );
    tree->SetBranchAddress("ak4jetPhiJERup"         , &ak4jetPhiJERup        , &b_ak4jetPhiJERup        );
    tree->SetBranchAddress("ak4jetPhiJERdown"       , &ak4jetPhiJERdown      , &b_ak4jetPhiJERdown      );
    tree->SetBranchAddress("ak4jetMassJERup"        , &ak4jetMassJERup       , &b_ak4jetMassJERup       );
    tree->SetBranchAddress("ak4jetMassJERdown"      , &ak4jetMassJERdown     , &b_ak4jetMassJERdown     );
    tree->SetBranchAddress("ak4jetHadronFlavour"    , &ak4jetHadronFlavour   , &b_ak4jetHadronFlavour   );
    tree->SetBranchAddress("ak8jetPtJERup"          , &ak8jetPtJERup         , &b_ak8jetPtJERup         );
    tree->SetBranchAddress("ak8jetPtJERdown"        , &ak8jetPtJERdown       , &b_ak8jetPtJERdown       );
    tree->SetBranchAddress("ak8jetEtaJERup"         , &ak8jetEtaJERup        , &b_ak8jetEtaJERup        );
    tree->SetBranchAddress("ak8jetEtaJERdown"       , &ak8jetEtaJERdown      , &b_ak8jetEtaJERdown      );
    tree->SetBranchAddress("ak8jetPhiJERup"         , &ak8jetPhiJERup        , &b_ak8jetPhiJERup        );
    tree->SetBranchAddress("ak8jetPhiJERdown"       , &ak8jetPhiJERdown      , &b_ak8jetPhiJERdown      );
    tree->SetBranchAddress("ak8jetMassJERup"        , &ak8jetMassJERup       , &b_ak8jetMassJERup       );
    tree->SetBranchAddress("ak8jetMassJERdown"      , &ak8jetMassJERdown     , &b_ak8jetMassJERdown     );

    if (sample.Contains("PowhegPythia8")){
      if (isSignal) tree->SetBranchAddress("truthChannel"           , &truthChannel        , &b_truthChannel        );     
      tree->SetBranchAddress("eventWeight_Q2Up"       , &eventWeight_Q2Up      , &b_eventWeight_Q2Up      );
      tree->SetBranchAddress("eventWeight_Q2Down"     , &eventWeight_Q2Down    , &b_eventWeight_Q2Down    );
      tree->SetBranchAddress("eventWeight_PDF"        , &eventWeight_PDF       , &b_eventWeight_PDF       );
      tree->SetBranchAddress("eventWeight_alphaUp"    , &eventWeight_alphaUp   , &b_eventWeight_alphaUp   );
      tree->SetBranchAddress("eventWeight_alphaDown"  , &eventWeight_alphaDown , &b_eventWeight_alphaDown );
    }
  }
  
  // parton level
  vector<float>* genTopPt       = 0;
  vector<float>* genTopEta      = 0;
  vector<float>* genTopPhi      = 0;
  vector<float>* genMuPt        = 0;
  vector<float>* genMuEta       = 0;
  vector<float>* genMuPhi       = 0;
  vector<float>* genElPt        = 0;
  vector<float>* genElEta       = 0;
  vector<float>* genElPhi       = 0;
  vector<float>* genTTbarMass   = 0;
  
  // particle level
  //vector<float>* genAK4jetPt    = 0;
  //vector<float>* genAK4jetEta   = 0;
  //vector<float>* genAK4jetPhi   = 0;
  //vector<float>* genAK4jetMass  = 0;
  vector<float>* genAK8jetPt    = 0;
  vector<float>* genAK8jetEta   = 0;
  vector<float>* genAK8jetPhi   = 0;
  vector<float>* genAK8jetMass  = 0;
  vector<float>* partMuPt       = 0;
  vector<float>* partMuEta      = 0;
  vector<float>* partMuPhi      = 0;
  vector<float>* partElPt       = 0;
  vector<float>* partElEta      = 0;
  vector<float>* partElPhi      = 0;
  
  if (isSignal) {
    // parton level
    TBranch* b_genTopPt;
    TBranch* b_genTopEta      ;
    TBranch* b_genTopPhi      ;
    TBranch* b_genMuPt        ;
    TBranch* b_genMuEta       ;
    TBranch* b_genMuPhi       ;
    TBranch* b_genElPt        ;
    TBranch* b_genElEta       ;
    TBranch* b_genElPhi       ;
    TBranch* b_genTTbarMass   ;
    
    // particle level
    //TBranch* b_genAK4jetPt    ;
    //TBranch* b_genAK4jetEta   ;
    //TBranch* b_genAK4jetPhi   ;
    //TBranch* b_genAK4jetMass  ;
    TBranch* b_genAK8jetPt    ;
    TBranch* b_genAK8jetEta   ;
    TBranch* b_genAK8jetPhi   ;
    TBranch* b_genAK8jetMass  ;
    TBranch* b_partMuPt       ;
    TBranch* b_partMuEta      ;
    TBranch* b_partMuPhi      ;
    TBranch* b_partElPt       ;
    TBranch* b_partElEta      ;
    TBranch* b_partElPhi      ;
    
    tree->SetBranchAddress("genTopPt"               , &genTopPt            , &b_genTopPt            );
    tree->SetBranchAddress("genTopEta"              , &genTopEta           , &b_genTopEta           );
    tree->SetBranchAddress("genTopPhi"              , &genTopPhi           , &b_genTopPhi           );
    tree->SetBranchAddress("genMuPt"                , &genMuPt             , &b_genMuPt             );
    tree->SetBranchAddress("genMuEta"               , &genMuEta            , &b_genMuEta            );
    tree->SetBranchAddress("genMuPhi"               , &genMuPhi            , &b_genMuPhi            );
    tree->SetBranchAddress("genElPt"                , &genElPt             , &b_genElPt             );
    tree->SetBranchAddress("genElEta"               , &genElEta            , &b_genElEta            );
    tree->SetBranchAddress("genElPhi"               , &genElPhi            , &b_genElPhi            );
    tree->SetBranchAddress("genTTbarMass"           , &genTTbarMass        , &b_genTTbarMass        );
    
    //tree->SetBranchAddress("genAK4jetPt"            , &genAK4jetPt         , &b_genAK4jetPt         );
    //tree->SetBranchAddress("genAK4jetEta"           , &genAK4jetEta        , &b_genAK4jetEta        );
    //tree->SetBranchAddress("genAK4jetPhi"           , &genAK4jetPhi        , &b_genAK4jetPhi        );
    //tree->SetBranchAddress("genAK4jetMass"          , &genAK4jetMass       , &b_genAK4jetMass       );
    tree->SetBranchAddress("genAK8jetPt"            , &genAK8jetPt         , &b_genAK8jetPt         );
    tree->SetBranchAddress("genAK8jetEta"           , &genAK8jetEta        , &b_genAK8jetEta        );
    tree->SetBranchAddress("genAK8jetPhi"           , &genAK8jetPhi        , &b_genAK8jetPhi        );
    tree->SetBranchAddress("genAK8jetMass"          , &genAK8jetMass       , &b_genAK8jetMass       );
    tree->SetBranchAddress("partMuPt"               , &partMuPt            , &b_partMuPt            );
    tree->SetBranchAddress("partMuEta"              , &partMuEta           , &b_partMuEta           );
    tree->SetBranchAddress("partMuPhi"              , &partMuPhi           , &b_partMuPhi           );
    tree->SetBranchAddress("partElPt"               , &partElPt            , &b_partElPt            );
    tree->SetBranchAddress("partElEta"              , &partElEta           , &b_partElEta           );
    tree->SetBranchAddress("partElPhi"              , &partElPhi           , &b_partElPhi           );

  }

  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  TH1F* h_nTrueLepBJet            = new TH1F("nTrueLepBJet" ,";Number of true b jets in lepton hemisphere;Events", 5,-0.5,4.5);

  TH1F* h_puWeight                = new TH1F("puWeight"     ,";PU weight;Events / 0.1"                    ,20,0.0,2.0);
  TH1F* h_lepSF                   = new TH1F("lepSF"        ,";Lepton SF;Events / 0.01"                   ,60,0.5,1.1);

  TH1F* h_genTopPt                = new TH1F("genTopPt"     ,";top quark p_{T} (GeV);Events / 20 GeV"     ,60,0.0,1200.);
  TH1F* h_genTopEta               = new TH1F("genTopEta"    ,";top quark #eta;Events / 0.1"               ,50,-2.5,2.5);
  TH1F* h_genTopPhi               = new TH1F("genTopPhi"    ,";top quark #phi;Events / 0.1"               ,70,-3.5,3.5);
  TH1F* h_genTTbarMass            = new TH1F("genTTbarMass" ,";top quark pair mass (GeV);Events / 10 GeV" ,50,0.0,500.);
  TH1F* h_genLepPt                = new TH1F("genLepPt"     ,";truth lepton p_{T} (GeV);Events / 10 GeV"  ,50,0.0,500.);
  TH1F* h_genLepEta               = new TH1F("genLepEta"    ,";truth lepton #eta;Events / 0.1"            ,50,-2.5,2.5);
  TH1F* h_genLepPhi               = new TH1F("genLepPhi"    ,";truth lepton #phi;Events / 0.1"            ,70,-3.5,3.5);
  
  TH1F* h_metPtPre                = new TH1F("metPtPre"                ,";Missing E_{T} (GeV);Events / 5 GeV"                       ,50,  0.0, 250.);
  TH1F* h_htPre                   = new TH1F("htPre"                   ,";H_{T} (GeV);Events / 20 GeV"                              ,60,400.0, 1600.);
  TH1F* h_htLepPre                = new TH1F("htLepPre"                ,";H_{T}^{lep} (GeV);Events / 20 GeV"                        ,50,  0.0, 1000.);
  TH1F* h_nAK4jetPre              = new TH1F("nAK4jetPre"              ,";# AK4 jets;Events"                                        ,10, -0.5, 9.5);
  TH1F* h_nBjetPre                = new TH1F("nBjetPre"                ,";# b-jet cands;Events"                                     ,5,  -0.5, 4.5);
  TH1F* h_ak4jetPtPre             = new TH1F("ak4jetPtPre"             ,";p_{T} of AK4 jets (GeV);Jets / 10 GeV"                    ,60,  0.0, 600.);
  TH1F* h_ak4jetEtaPre            = new TH1F("ak4jetEtaPre"            ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetEtaPre_bb         = new TH1F("ak4jetEtaPre_bb"         ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetEtaPre_b          = new TH1F("ak4jetEtaPre_b"          ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetEtaPre_cc         = new TH1F("ak4jetEtaPre_cc"         ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetEtaPre_c          = new TH1F("ak4jetEtaPre_c"          ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetEtaPre_l          = new TH1F("ak4jetEtaPre_l"          ,";#eta of AK4 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak4jetPhiPre            = new TH1F("ak4jetPhiPre"            ,";#phi of AK4 jets;Jets / 0.1"                              ,70, -3.5, 3.5);
  TH1F* h_ak4jetMassPre           = new TH1F("ak4jetMassPre"           ,";Mass of AK4 jets (GeV);Jets / 2 GeV"                      ,50,  0.0, 100.);
  TH1F* h_ak4jetCSVPre            = new TH1F("ak4jetCSVPre"            ,";CSVv2 of AK4 jets;Jets / 0.02"                            ,50,  0.0, 1.0);
  TH1F* h_ak4jetVtxMassPre        = new TH1F("ak4jetVtxMassPre"        ,";Secondary vertex mass of AK4 jets (GeV);Jets / 0.1 GeV"   ,50,  0.0, 5.0);
  TH1F* h_nAK8jetPre              = new TH1F("nAK8jetPre"              ,";# AK8 jets;Events"                                        ,10, -0.5, 9.5);
  TH1F* h_nTjetPre                = new TH1F("nTjetPre"                ,";# t-jet cands;Events"                                     ,5,  -0.5, 4.5);
  TH1F* h_ak8jetPtPre             = new TH1F("ak8jetPtPre"             ,";p_{T} of AK8 jets (GeV);Jets / 20 GeV"                    ,60,  0.0, 1200.);
  TH1F* h_ak8jetEtaPre            = new TH1F("ak8jetEtaPre"            ,";#eta of AK8 jets;Jets / 0.1"                              ,50, -2.5, 2.5);
  TH1F* h_ak8jetPhiPre            = new TH1F("ak8jetPhiPre"            ,";#phi of AK8 jets;Jets / 0.1"                              ,70, -3.5, 3.5);
  TH1F* h_ak8jetYPre              = new TH1F("ak8jetYPre"              ,";Rapidity of AK8 jets;Jets / 0.1"                          ,50, -2.5, 2.5);
  TH1F* h_ak8jetMassPre           = new TH1F("ak8jetMassPre"           ,";Mass of AK8 jets (GeV);Jets / 5 GeV"                      ,50,  0.0, 250.);
  TH1F* h_ak8jetMassPrunedPre     = new TH1F("ak8jetMassPrunedPre"     ,";Pruned mass of AK8 jets (GeV);Jets / 5 GeV"               ,50,  0.0, 250.);
  TH1F* h_ak8jetMassFilteredPre   = new TH1F("ak8jetMassFilteredPre"   ,";Filtered mass of AK8 jets (GeV);Jets / 5 GeV"             ,50,  0.0, 250.);
  TH1F* h_ak8jetMassTrimmedPre    = new TH1F("ak8jetMassTrimmedPre"    ,";Trimmed mass of AK8 jets (GeV);Jets / 5 Gev"              ,50,  0.0, 250.);
  TH1F* h_ak8jetTau1Pre           = new TH1F("ak8jetTau1Pre"           ,";#tau_{1} of AK8 jets;Jets / 0.01"                         ,70,  0.0, 0.7);
  TH1F* h_ak8jetTau2Pre           = new TH1F("ak8jetTau2Pre"           ,";#tau_{2} of AK8 jets;Jets / 0.01"                         ,50,  0.0, 0.5);
  TH1F* h_ak8jetTau3Pre           = new TH1F("ak8jetTau3Pre"           ,";#tau_{3} of AK8 jets;Jets / 0.01"                         ,30,  0.0, 0.3);
  TH1F* h_ak8jetTau32Pre          = new TH1F("ak8jetTau32Pre"          ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau32Pre_bb       = new TH1F("ak8jetTau32Pre_bb"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau32Pre_b        = new TH1F("ak8jetTau32Pre_b"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau32Pre_cc       = new TH1F("ak8jetTau32Pre_cc"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau32Pre_c        = new TH1F("ak8jetTau32Pre_c"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau32Pre_l        = new TH1F("ak8jetTau32Pre_l"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre          = new TH1F("ak8jetTau21Pre"          ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre_bb       = new TH1F("ak8jetTau21Pre_bb"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre_b        = new TH1F("ak8jetTau21Pre_b"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre_cc       = new TH1F("ak8jetTau21Pre_cc"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre_c        = new TH1F("ak8jetTau21Pre_c"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau21Pre_l        = new TH1F("ak8jetTau21Pre_l"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                ,100, 0.0, 1.0);
  TH1F* h_ak8jetCSVPre            = new TH1F("ak8jetCSVPre"            ,";CSVv2 of AK8 jets;Jets / 0.02"                            ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDmassPre         = new TH1F("ak8jetSDmassPre"         ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmassPre_bb      = new TH1F("ak8jetSDmassPre_bb"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmassPre_b       = new TH1F("ak8jetSDmassPre_b"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmassPre_cc      = new TH1F("ak8jetSDmassPre_cc"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmassPre_c       = new TH1F("ak8jetSDmassPre_c"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmassPre_l       = new TH1F("ak8jetSDmassPre_l"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"            ,50,  0.0, 250.);
  TH1F* h_ak8jetSDsubjet0ptPre    = new TH1F("ak8jetSDsubjet0ptPre"    ,";p_{T} of leading Soft Drop subjet (GeV);Jets / 10 GeV"    ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet0massPre  = new TH1F("ak8jetSDsubjet0massPre"  ,";Mass of leading Soft Drop subjet (GeV);Jets / 5 GeV"      ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet0CSVPre   = new TH1F("ak8jetSDsubjet0CSVPre"   ,";CSVv2 of leading Soft Drop subjet;Jets / 0.02"            ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjet1ptPre    = new TH1F("ak8jetSDsubjet1ptPre"    ,";p_{T} of subleading Soft Drop subjet (GeV);Jets / 10 GeV" ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet1massPre  = new TH1F("ak8jetSDsubjet1massPre"  ,";Mass of subleading Soft Drop subjet (GeV);Jets / 5 GeV"   ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet1CSVPre   = new TH1F("ak8jetSDsubjet1CSVPre"   ,";CSVv2 of subleading Soft Drop subjet;Jets / 0.02"         ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjetMaxCSVPre = new TH1F("ak8jetSDsubjetMaxCSVPre" ,";Max CSV of Soft Drop subjets;Jets / 0.02"                 ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDm01Pre          = new TH1F("ak8jetSDm01Pre"          ,";Mass of two hardest Soft Drop subjets (GeV);Jets / 5 GeV" ,50,  0.0, 250.);
  TH1F* h_lepPtPre                = new TH1F("lepPtPre"                ,";Lepton p_{T} (GeV);Leptons / 10 GeV"                      ,50,  0.0, 500.);
  TH1F* h_lepEtaPre               = new TH1F("lepEtaPre"               ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepEtaPre_bb            = new TH1F("lepEtaPre_bb"            ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepEtaPre_b             = new TH1F("lepEtaPre_b"             ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepEtaPre_cc            = new TH1F("lepEtaPre_cc"            ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepEtaPre_c             = new TH1F("lepEtaPre_c"             ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepEtaPre_l             = new TH1F("lepEtaPre_l"             ,";Lepton #eta;Leptons / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_lepAbsEtaPre            = new TH1F("lepAbsEtaPre"            ,";Lepton |#eta|;Leptons / 0.1"                              ,25,  0.0, 2.5);
  TH1F* h_lepSignEtaPre           = new TH1F("lepSignEtaPre"           ,";Lepton |#eta| * Q;Leptons / 0.1"                          ,50, -2.5, 2.5);
  TH1F* h_lepPhiPre               = new TH1F("lepPhiPre"               ,";Lepton #phi;Leptons / 0.1"                                ,70, -3.5, 3.5);
  TH1F* h_lepBJetdRPre            = new TH1F("lepBJetdRPre"            ,";#Delta R(lepton, b-jet cand);Leptons / 0.1"               ,35,  0.0, 3.5);
  TH1F* h_lepTJetdRPre            = new TH1F("lepTJetdRPre"            ,";#Delta R(lepton, t-jet cand);Leptons / 0.1"               ,35,  0.0, 3.5);
  TH1F* h_lepBJetPtRelPre         = new TH1F("lepBJetPtRelPre"         ,";p_{T}^{rel}(lepton, b-jet cand) (GeV);Leptons / 2 GeV"    ,50,  0.0, 100.);
  TH1F* h_lepMiniIsoPre           = new TH1F("lepMiniIsoPre"           ,";lepton miniIso;Leptons / 0.02"                            ,50,  0.0, 1.0);
  TH1F* h_leadLepPtPre            = new TH1F("leadLepPtPre"            ,";Leading lepton p_{T} (GeV);Events / 10 GeV"               ,50,  0.0, 500.);

  TH1F* h_metPt1b                = new TH1F("metPt1b"                ,";Missing E_{T} (GeV);Events / 5 GeV"                              ,50,  0.0, 250.);
  TH1F* h_ht1b                   = new TH1F("ht1b"                   ,";H_{T} (GeV);Events / 20 GeV"                                     ,60,400.0, 1600.);
  TH1F* h_htLep1b                = new TH1F("htLep1b"                ,";H_{T}^{lep} (GeV);Events / 20 GeV"                               ,50,  0.0, 1000.);
  TH1F* h_nAK4jet1b              = new TH1F("nAK4jet1b"              ,";# AK4 jets;Events"                                               ,10, -0.5, 9.5);
  TH1F* h_nBjet1b                = new TH1F("nBjet1b"                ,";# b-jet cands;Events"                                            ,5,  -0.5, 4.5);
  TH1F* h_ak4jetPt1b             = new TH1F("ak4jetPt1b"             ,";p_{T} of b-jet candidate (GeV);Events / 10 GeV"                  ,60,  0.0, 600.);
  TH1F* h_ak4jetEta1b            = new TH1F("ak4jetEta1b"            ,";#eta of b-jet candidate;Events / 0.1"                            ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1b_bb         = new TH1F("ak4jetEta1b_bb"         ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1b_b          = new TH1F("ak4jetEta1b_b"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1b_cc         = new TH1F("ak4jetEta1b_cc"         ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1b_c          = new TH1F("ak4jetEta1b_c"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1b_l          = new TH1F("ak4jetEta1b_l"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetPhi1b            = new TH1F("ak4jetPhi1b"            ,";#phi of b-jet candidate;Events / 0.1"                            ,70, -3.5, 3.5);
  TH1F* h_ak4jetMass1b           = new TH1F("ak4jetMass1b"           ,";Mass of b-jet candidate (GeV);Events / 2 GeV"                    ,50,  0.0, 100.);
  TH1F* h_ak4jetCSV1b            = new TH1F("ak4jetCSV1b"            ,";CSVv2 of b-jet candidate;Events / 0.02"                          ,50,  0.0, 1.0);
  TH1F* h_ak4jetVtxMass1b        = new TH1F("ak4jetVtxMass1b"        ,";Secondary vertex mass of b-jet candidate (GeV);Events / 0.1 GeV" ,50,  0.0, 5.0);
  TH1F* h_nAK8jet1b              = new TH1F("nAK8jet1b"              ,";# AK8 jets;Events"                                               ,10, -0.5, 9.5);
  TH1F* h_nTjet1b                = new TH1F("nTjet1b"                ,";# t-jet cands;Events"                                            ,5,  -0.5, 4.5);
  TH1F* h_ak8jetPt1b             = new TH1F("ak8jetPt1b"             ,";p_{T} of t-jet candidate (GeV);Events / 20 GeV"                  ,60,  0.0, 1200.);
  TH1F* h_ak8jetEta1b            = new TH1F("ak8jetEta1b"            ,";#eta of t-jet candidate;Events / 0.1"                            ,50, -2.5, 2.5);
  TH1F* h_ak8jetPhi1b            = new TH1F("ak8jetPhi1b"            ,";#phi of t-jet candidate;Events / 0.1"                            ,70, -3.5, 3.5);
  TH1F* h_ak8jetY1b              = new TH1F("ak8jetY1b"              ,";Rapidity of t-jet candidate;Events / 0.1"                        ,50, -2.5, 2.5);
  TH1F* h_ak8jetMass1b           = new TH1F("ak8jetMass1b"           ,";Mass of t-jet candidate (GeV);Events / 5 GeV"                    ,50,  0.0, 250.);
  TH1F* h_ak8jetMassPruned1b     = new TH1F("ak8jetMassPruned1b"     ,";Pruned mass of t-jet candidate (GeV);Events / 5 GeV"             ,50,  0.0, 250.);
  TH1F* h_ak8jetMassFiltered1b   = new TH1F("ak8jetMassFiltered1b"   ,";Filtered mass of t-jet candidate (GeV);Events / 5 GeV"           ,50,  0.0, 250.);
  TH1F* h_ak8jetMassTrimmed1b    = new TH1F("ak8jetMassTrimmed1b"    ,";Trimmed mass of t-jet candidate (GeV);Events / 5 Gev"            ,50,  0.0, 250.);
  TH1F* h_ak8jetTau11b           = new TH1F("ak8jetTau11b"           ,";#tau_{1} of t-jet candidate;Events / 0.01"                       ,70,  0.0, 0.7);
  TH1F* h_ak8jetTau21b           = new TH1F("ak8jetTau21b"           ,";#tau_{2} of t-jet candidate;Events / 0.01"                       ,50,  0.0, 0.5);
  TH1F* h_ak8jetTau31b           = new TH1F("ak8jetTau31b"           ,";#tau_{3} of t-jet candidate;Events / 0.01"                       ,30,  0.0, 0.3);
  TH1F* h_ak8jetTau321b          = new TH1F("ak8jetTau321b"          ,";#tau_{3}/#tau_{2} of t-jet candidate;Events / 0.01"              ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321b_bb       = new TH1F("ak8jetTau321b_bb"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321b_b        = new TH1F("ak8jetTau321b_b"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321b_cc       = new TH1F("ak8jetTau321b_cc"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321b_c        = new TH1F("ak8jetTau321b_c"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321b_l        = new TH1F("ak8jetTau321b_l"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b          = new TH1F("ak8jetTau211b"          ,";#tau_{2}/#tau_{1} of t-jet candidate;Events / 0.01"              ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b_bb       = new TH1F("ak8jetTau211b_bb"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b_b        = new TH1F("ak8jetTau211b_b"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b_cc       = new TH1F("ak8jetTau211b_cc"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b_c        = new TH1F("ak8jetTau211b_c"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211b_l        = new TH1F("ak8jetTau211b_l"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetCSV1b            = new TH1F("ak8jetCSV1b"            ,";CSVv2 of t-jet candidate;Events / 0.02"                          ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDmass1b         = new TH1F("ak8jetSDmass1b"         ,";Soft Drop mass of t-jet candidate (GeV);Events / 5 GeV"          ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1b_bb      = new TH1F("ak8jetSDmass1b_bb"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1b_b       = new TH1F("ak8jetSDmass1b_b"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1b_cc      = new TH1F("ak8jetSDmass1b_cc"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1b_c       = new TH1F("ak8jetSDmass1b_c"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1b_l       = new TH1F("ak8jetSDmass1b_l"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDsubjet0pt1b    = new TH1F("ak8jetSDsubjet0pt1b"    ,";p_{T} of leading Soft Drop subjet (GeV);Events / 10 GeV"         ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet0mass1b  = new TH1F("ak8jetSDsubjet0mass1b"  ,";Mass of leading Soft Drop subjet (GeV);Events / 5 GeV"           ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet0CSV1b   = new TH1F("ak8jetSDsubjet0CSV1b"   ,";CSVv2 of leading Soft Drop subjet;Events / 0.02"                 ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjet1pt1b    = new TH1F("ak8jetSDsubjet1pt1b"    ,";p_{T} of subleading Soft Drop subjet (GeV);Events / 10 GeV"      ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet1mass1b  = new TH1F("ak8jetSDsubjet1mass1b"  ,";Mass of subleading Soft Drop subjet (GeV);Events / 5 GeV"        ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet1CSV1b   = new TH1F("ak8jetSDsubjet1CSV1b"   ,";CSVv2 of subleading Soft Drop subjet;Events / 0.02"              ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjetMaxCSV1b = new TH1F("ak8jetSDsubjetMaxCSV1b" ,";Max CSV of Soft Drop subjets;Events / 0.02"                      ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDm011b          = new TH1F("ak8jetSDm011b"          ,";Mass of two hardest Soft Drop subjets (GeV);Events / 5 GeV"      ,50,  0.0, 250.);
  TH1F* h_lepPt1b                = new TH1F("lepPt1b"                ,";Lepton p_{T} (GeV);Events / 10 GeV"                              ,50,  0.0, 500.);
  TH1F* h_lepEta1b               = new TH1F("lepEta1b"               ,";Lepton #eta;Events / 0.1"                                        ,50, -2.5, 2.5);
  TH1F* h_lepEta1b_bb            = new TH1F("lepEta1b_bb"            ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1b_b             = new TH1F("lepEta1b_b"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1b_cc            = new TH1F("lepEta1b_cc"            ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1b_c             = new TH1F("lepEta1b_c"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1b_l             = new TH1F("lepEta1b_l"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepAbsEta1b            = new TH1F("lepAbsEta1b"            ,";Lepton |#eta|;Events / 0.1"                                      ,25,  0.0, 2.5);
  TH1F* h_lepSignEta1b           = new TH1F("lepSignEta1b"           ,";Lepton |#eta| * Q;Events / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepPhi1b               = new TH1F("lepPhi1b"               ,";Lepton #phi;Events / 0.1"                                        ,70, -3.5, 3.5);
  TH1F* h_lepBJetdR1b            = new TH1F("lepBJetdR1b"            ,";#Delta R(lepton, b-jet cand);Events / 0.1"                       ,35,  0.0, 3.5);
  TH1F* h_lepTJetdR1b            = new TH1F("lepTJetdR1b"            ,";#Delta R(lepton, t-jet cand);Events / 0.1"                       ,35,  0.0, 3.5);
  TH1F* h_lepBJetPtRel1b         = new TH1F("lepBJetPtRel1b"         ,";p_{T}^{rel}(lepton, b-jet cand) (GeV);Events / 2 GeV"            ,50,  0.0, 100.);
  TH1F* h_lepMiniIso1b           = new TH1F("lepMiniIso1b"           ,";lepton miniIso;Events / 0.02"                                    ,50,  0.0, 1.0);

  TH1F* h_metPt1t                = new TH1F("metPt1t"                ,";Missing E_{T} (GeV);Events / 5 GeV"                              ,50,  0.0, 250.);
  TH1F* h_ht1t                   = new TH1F("ht1t"                   ,";H_{T} (GeV);Events / 20 GeV"                                     ,60,400.0, 1600.);
  TH1F* h_htLep1t                = new TH1F("htLep1t"                ,";H_{T}^{lep} (GeV);Events / 20 GeV"                               ,50,  0.0, 1000.);
  TH1F* h_nAK4jet1t              = new TH1F("nAK4jet1t"              ,";# AK4 jets;Events"                                               ,10, -0.5, 9.5);
  TH1F* h_nBjet1t                = new TH1F("nBjet1t"                ,";# b-jet cands;Events"                                            ,5,  -0.5, 4.5);
  TH1F* h_ak4jetPt1t             = new TH1F("ak4jetPt1t"             ,";p_{T} of b-jet candidate (GeV);Events / 10 GeV"                  ,60,  0.0, 600.);
  TH1F* h_ak4jetEta1t            = new TH1F("ak4jetEta1t"            ,";#eta of b-jet candidate;Events / 0.1"                            ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t_bb         = new TH1F("ak4jetEta1t_bb"         ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t_b          = new TH1F("ak4jetEta1t_b"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t_cc         = new TH1F("ak4jetEta1t_cc"         ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t_c          = new TH1F("ak4jetEta1t_c"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t_l          = new TH1F("ak4jetEta1t_l"          ,";#eta of AK4 jets;Jets / 0.1"                                     ,50, -2.5, 2.5);
  TH1F* h_ak4jetPhi1t            = new TH1F("ak4jetPhi1t"            ,";#phi of b-jet candidate;Events / 0.1"                            ,70, -3.5, 3.5);
  TH1F* h_ak4jetMass1t           = new TH1F("ak4jetMass1t"           ,";Mass of b-jet candidate (GeV);Events / 2 GeV"                    ,50,  0.0, 100.);
  TH1F* h_ak4jetCSV1t            = new TH1F("ak4jetCSV1t"            ,";CSVv2 of b-jet candidate;Events / 0.02"                          ,50,  0.0, 1.0);
  TH1F* h_ak4jetVtxMass1t        = new TH1F("ak4jetVtxMass1t"        ,";Secondary vertex mass of b-jet candidate (GeV);Events / 0.1 GeV" ,50,  0.0, 5.0);
  TH1F* h_nAK8jet1t              = new TH1F("nAK8jet1t"              ,";# AK8 jets;Events"                                               ,10, -0.5, 9.5);
  TH1F* h_nTjet1t                = new TH1F("nTjet1t"                ,";# t-jet cands;Events"                                            ,5,  -0.5, 4.5);
  TH1F* h_ak8jetPt1t             = new TH1F("ak8jetPt1t"             ,";p_{T} of t-jet (GeV);Events / 20 GeV"                            ,60,  0.0, 1200.);
  TH1F* h_ak8jetEta1t            = new TH1F("ak8jetEta1t"            ,";#eta of t-jet;Events / 0.1"                                      ,50, -2.5, 2.5);
  TH1F* h_ak8jetPhi1t            = new TH1F("ak8jetPhi1t"            ,";#phi of t-jet;Events / 0.1"                                      ,70, -3.5, 3.5);
  TH1F* h_ak8jetY1t              = new TH1F("ak8jetY1t"              ,";Rapidity of t-jet;Events / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_ak8jetMass1t           = new TH1F("ak8jetMass1t"           ,";Mass of t-jet (GeV);Events / 5 GeV"                              ,50,  0.0, 250.);
  TH1F* h_ak8jetMassPruned1t     = new TH1F("ak8jetMassPruned1t"     ,";Pruned mass of t-jet (GeV);Events / 5 GeV"                       ,50,  0.0, 250.);
  TH1F* h_ak8jetMassFiltered1t   = new TH1F("ak8jetMassFiltered1t"   ,";Filtered mass of t-jet (GeV);Events / 5 GeV"                     ,50,  0.0, 250.);
  TH1F* h_ak8jetMassTrimmed1t    = new TH1F("ak8jetMassTrimmed1t"    ,";Trimmed mass of t-jet (GeV);Events / 5 Gev"                      ,50,  0.0, 250.);
  TH1F* h_ak8jetTau11t           = new TH1F("ak8jetTau11t"           ,";#tau_{1} of t-jet;Events / 0.01"                                 ,70,  0.0, 0.7);
  TH1F* h_ak8jetTau21t           = new TH1F("ak8jetTau21t"           ,";#tau_{2} of t-jet;Events / 0.01"                                 ,50,  0.0, 0.5);
  TH1F* h_ak8jetTau31t           = new TH1F("ak8jetTau31t"           ,";#tau_{3} of t-jet;Events / 0.01"                                 ,30,  0.0, 0.3);
  TH1F* h_ak8jetTau321t          = new TH1F("ak8jetTau321t"          ,";#tau_{3}/#tau_{2} of t-jet;Events / 0.01"                        ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t_bb       = new TH1F("ak8jetTau321t_bb"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t_b        = new TH1F("ak8jetTau321t_b"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t_cc       = new TH1F("ak8jetTau321t_cc"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t_c        = new TH1F("ak8jetTau321t_c"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t_l        = new TH1F("ak8jetTau321t_l"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t          = new TH1F("ak8jetTau211t"          ,";#tau_{2}/#tau_{1} of t-jet;Events / 0.01"                        ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t_bb       = new TH1F("ak8jetTau211t_bb"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t_b        = new TH1F("ak8jetTau211t_b"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t_cc       = new TH1F("ak8jetTau211t_cc"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t_c        = new TH1F("ak8jetTau211t_c"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t_l        = new TH1F("ak8jetTau211t_l"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                       ,100, 0.0, 1.0);
  TH1F* h_ak8jetCSV1t            = new TH1F("ak8jetCSV1t"            ,";CSVv2 of t-jet;Events / 0.02"                                    ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDmass1t         = new TH1F("ak8jetSDmass1t"         ,";Soft Drop mass of t-jet (GeV);Events / 5 GeV"                    ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t_bb      = new TH1F("ak8jetSDmass1t_bb"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t_b       = new TH1F("ak8jetSDmass1t_b"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t_cc      = new TH1F("ak8jetSDmass1t_cc"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t_c       = new TH1F("ak8jetSDmass1t_c"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t_l       = new TH1F("ak8jetSDmass1t_l"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"                   ,50,  0.0, 250.);
  TH1F* h_ak8jetSDsubjet0pt1t    = new TH1F("ak8jetSDsubjet0pt1t"    ,";p_{T} of leading Soft Drop subjet (GeV);Events / 10 GeV"         ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet0mass1t  = new TH1F("ak8jetSDsubjet0mass1t"  ,";Mass of leading Soft Drop subjet (GeV);Events / 5 GeV"           ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet0CSV1t   = new TH1F("ak8jetSDsubjet0CSV1t"   ,";CSVv2 of leading Soft Drop subjet;Events / 0.02"                 ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjet1pt1t    = new TH1F("ak8jetSDsubjet1pt1t"    ,";p_{T} of subleading Soft Drop subjet (GeV);Events / 10 GeV"      ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet1mass1t  = new TH1F("ak8jetSDsubjet1mass1t"  ,";Mass of subleading Soft Drop subjet (GeV);Events / 5 GeV"        ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet1CSV1t   = new TH1F("ak8jetSDsubjet1CSV1t"   ,";CSVv2 of subleading Soft Drop subjet;Events / 0.02"              ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjetMaxCSV1t = new TH1F("ak8jetSDsubjetMaxCSV1t" ,";Max CSV of Soft Drop subjets;Events / 0.02"                      ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDm011t          = new TH1F("ak8jetSDm011t"          ,";Mass of two hardest Soft Drop subjets (GeV);Events / 5 GeV"      ,50,  0.0, 250.);
  TH1F* h_lepPt1t                = new TH1F("lepPt1t"                ,";Lepton p_{T} (GeV);Events / 10 GeV"                              ,50,  0.0, 500.);
  TH1F* h_lepEta1t               = new TH1F("lepEta1t"               ,";Lepton #eta;Events / 0.1"                                        ,50, -2.5, 2.5);
  TH1F* h_lepEta1t_bb            = new TH1F("lepEta1t_bb"            ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1t_b             = new TH1F("lepEta1t_b"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1t_cc            = new TH1F("lepEta1t_cc"            ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1t_c             = new TH1F("lepEta1t_c"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepEta1t_l             = new TH1F("lepEta1t_l"             ,";Lepton #eta;Leptons / 0.1"                                       ,50, -2.5, 2.5);
  TH1F* h_lepAbsEta1t            = new TH1F("lepAbsEta1t"            ,";Lepton |#eta|;Events / 0.1"                                      ,25,  0.0, 2.5);
  TH1F* h_lepSignEta1t           = new TH1F("lepSignEta1t"           ,";Lepton |#eta| * Q;Events / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepPhi1t               = new TH1F("lepPhi1t"               ,";Lepton #phi;Events / 0.1"                                        ,70, -3.5, 3.5);
  TH1F* h_lepBJetdR1t            = new TH1F("lepBJetdR1t"            ,";#Delta R(lepton, b-jet cand);Events / 0.1"                       ,35,  0.0, 3.5);
  TH1F* h_lepTJetdR1t            = new TH1F("lepTJetdR1t"            ,";#Delta R(lepton, t-jet);Events / 0.1"                            ,35,  0.0, 3.5);
  TH1F* h_lepBJetPtRel1t         = new TH1F("lepBJetPtRel1t"         ,";p_{T}^{rel}(lepton, b-jet cand) (GeV);Events / 2 GeV"            ,50,  0.0, 100.);
  TH1F* h_lepMiniIso1t           = new TH1F("lepMiniIso1t"           ,";lepton miniIso;Events / 0.02"                                    ,50,  0.0, 1.0);

  TH1F* h_metPt1t1b                = new TH1F("metPt1t1b"                ,";Missing E_{T} (GeV);Events / 5 GeV"                         ,50,  0.0, 250.);
  TH1F* h_ht1t1b                   = new TH1F("ht1t1b"                   ,";H_{T} (GeV);Events / 20 GeV"                                ,60,400.0, 1600.);
  TH1F* h_htLep1t1b                = new TH1F("htLep1t1b"                ,";H_{T}^{lep} (GeV);Events / 20 GeV"                          ,50,  0.0, 1000.);
  TH1F* h_nAK4jet1t1b              = new TH1F("nAK4jet1t1b"              ,";# AK4 jets;Events"                                          ,10, -0.5, 9.5);
  TH1F* h_nBjet1t1b                = new TH1F("nBjet1t1b"                ,";# b-jet cands;Events"                                       ,5,  -0.5, 4.5);
  TH1F* h_ak4jetPt1t1b             = new TH1F("ak4jetPt1t1b"             ,";p_{T} of b-jet (GeV);Events / 10 GeV"                       ,60,  0.0, 600.);
  TH1F* h_ak4jetEta1t1b            = new TH1F("ak4jetEta1t1b"            ,";#eta of b-jet;Events / 0.1"                                 ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t1b_bb         = new TH1F("ak4jetEta1t1b_bb"         ,";#eta of AK4 jets;Jets / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t1b_b          = new TH1F("ak4jetEta1t1b_b"          ,";#eta of AK4 jets;Jets / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t1b_cc         = new TH1F("ak4jetEta1t1b_cc"         ,";#eta of AK4 jets;Jets / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t1b_c          = new TH1F("ak4jetEta1t1b_c"          ,";#eta of AK4 jets;Jets / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_ak4jetEta1t1b_l          = new TH1F("ak4jetEta1t1b_l"          ,";#eta of AK4 jets;Jets / 0.1"                                ,50, -2.5, 2.5);
  TH1F* h_ak4jetPhi1t1b            = new TH1F("ak4jetPhi1t1b"            ,";#phi of b-jet;Events / 0.1"                                 ,70, -3.5, 3.5);
  TH1F* h_ak4jetMass1t1b           = new TH1F("ak4jetMass1t1b"           ,";Mass of b-jet (GeV);Events / 2 GeV"                         ,50,  0.0, 100.);
  TH1F* h_ak4jetCSV1t1b            = new TH1F("ak4jetCSV1t1b"            ,";CSVv2 of b-jet;Events / 0.02"                               ,50,  0.0, 1.0);
  TH1F* h_ak4jetVtxMass1t1b        = new TH1F("ak4jetVtxMass1t1b"        ,";Secondary vertex mass of b-jet (GeV);Events / 0.1 GeV"      ,50,  0.0, 5.0);
  TH1F* h_nAK8jet1t1b              = new TH1F("nAK8jet1t1b"              ,";# AK8 jets;Events"                                          ,10, -0.5, 9.5);
  TH1F* h_nTjet1t1b                = new TH1F("nTjet1t1b"                ,";# t-jet cands;Events"                                       ,5,  -0.5, 4.5);
  TH1F* h_ak8jetPt1t1b             = new TH1F("ak8jetPt1t1b"             ,";p_{T} of t-jet (GeV);Events / 20 GeV"                       ,60,  0.0, 1200.);
  TH1F* h_ak8jetEta1t1b            = new TH1F("ak8jetEta1t1b"            ,";#eta of t-jet;Events / 0.1"                                 ,50, -2.5, 2.5);
  TH1F* h_ak8jetPhi1t1b            = new TH1F("ak8jetPhi1t1b"            ,";#phi of t-jet;Events / 0.1"                                 ,70, -3.5, 3.5);
  TH1F* h_ak8jetY1t1b              = new TH1F("ak8jetY1t1b"              ,";Rapidity of t-jet;Events / 0.1"                             ,50, -2.5, 2.5);
  TH1F* h_ak8jetMass1t1b           = new TH1F("ak8jetMass1t1b"           ,";Mass of t-jet (GeV);Events / 5 GeV"                         ,50,  0.0, 250.);
  TH1F* h_ak8jetMassPruned1t1b     = new TH1F("ak8jetMassPruned1t1b"     ,";Pruned mass of t-jet (GeV);Events / 5 GeV"                  ,50,  0.0, 250.);
  TH1F* h_ak8jetMassFiltered1t1b   = new TH1F("ak8jetMassFiltered1t1b"   ,";Filtered mass of t-jet (GeV);Events / 5 GeV"                ,50,  0.0, 250.);
  TH1F* h_ak8jetMassTrimmed1t1b    = new TH1F("ak8jetMassTrimmed1t1b"    ,";Trimmed mass of t-jet (GeV);Events / 5 Gev"                 ,50,  0.0, 250.);
  TH1F* h_ak8jetTau11t1b           = new TH1F("ak8jetTau11t1b"           ,";#tau_{1} of t-jet;Events / 0.01"                            ,70,  0.0, 0.7);
  TH1F* h_ak8jetTau21t1b           = new TH1F("ak8jetTau21t1b"           ,";#tau_{2} of t-jet;Events / 0.01"                            ,50,  0.0, 0.5);
  TH1F* h_ak8jetTau31t1b           = new TH1F("ak8jetTau31t1b"           ,";#tau_{3} of t-jet;Events / 0.01"                            ,30,  0.0, 0.3);
  TH1F* h_ak8jetTau321t1b          = new TH1F("ak8jetTau321t1b"          ,";#tau_{3}/#tau_{2} of t-jet;Events / 0.01"                   ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t1b_bb       = new TH1F("ak8jetTau321t1b_bb"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t1b_b        = new TH1F("ak8jetTau321t1b_b"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t1b_cc       = new TH1F("ak8jetTau321t1b_cc"       ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t1b_c        = new TH1F("ak8jetTau321t1b_c"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau321t1b_l        = new TH1F("ak8jetTau321t1b_l"        ,";#tau_{3}/#tau_{2} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b          = new TH1F("ak8jetTau211t1b"          ,";#tau_{2}/#tau_{1} of t-jet;Events / 0.01"                   ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b_bb       = new TH1F("ak8jetTau211t1b_bb"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b_b        = new TH1F("ak8jetTau211t1b_b"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b_cc       = new TH1F("ak8jetTau211t1b_cc"       ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b_c        = new TH1F("ak8jetTau211t1b_c"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetTau211t1b_l        = new TH1F("ak8jetTau211t1b_l"        ,";#tau_{2}/#tau_{1} of AK8 jets;Jets / 0.01"                  ,100, 0.0, 1.0);
  TH1F* h_ak8jetCSV1t1b            = new TH1F("ak8jetCSV1t1b"            ,";CSVv2 of t-jet;Events / 0.02"                               ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDmass1t1b         = new TH1F("ak8jetSDmass1t1b"         ,";Soft Drop mass of t-jet (GeV);Events / 5 GeV"               ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t1b_bb      = new TH1F("ak8jetSDmass1t1b_bb"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"              ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t1b_b       = new TH1F("ak8jetSDmass1t1b_b"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"              ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t1b_cc      = new TH1F("ak8jetSDmass1t1b_cc"      ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"              ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t1b_c       = new TH1F("ak8jetSDmass1t1b_c"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"              ,50,  0.0, 250.);
  TH1F* h_ak8jetSDmass1t1b_l       = new TH1F("ak8jetSDmass1t1b_l"       ,";Soft Drop mass of AK8 jets (GeV);Jets / 5 GeV"              ,50,  0.0, 250.);
  TH1F* h_ak8jetSDsubjet0pt1t1b    = new TH1F("ak8jetSDsubjet0pt1t1b"    ,";p_{T} of leading Soft Drop subjet (GeV);Events / 10 GeV"    ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet0mass1t1b  = new TH1F("ak8jetSDsubjet0mass1t1b"  ,";Mass of leading Soft Drop subjet (GeV);Events / 5 GeV"      ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet0CSV1t1b   = new TH1F("ak8jetSDsubjet0CSV1t1b"   ,";CSVv2 of leading Soft Drop subjet;Events / 0.02"            ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjet1pt1t1b    = new TH1F("ak8jetSDsubjet1pt1t1b"    ,";p_{T} of subleading Soft Drop subjet (GeV);Events / 10 GeV" ,80,  0.0, 800.);
  TH1F* h_ak8jetSDsubjet1mass1t1b  = new TH1F("ak8jetSDsubjet1mass1t1b"  ,";Mass of subleading Soft Drop subjet (GeV);Events / 5 GeV"   ,40,  0.0, 200.);
  TH1F* h_ak8jetSDsubjet1CSV1t1b   = new TH1F("ak8jetSDsubjet1CSV1t1b"   ,";CSVv2 of subleading Soft Drop subjet;Events / 0.02"         ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDsubjetMaxCSV1t1b = new TH1F("ak8jetSDsubjetMaxCSV1t1b" ,";Max CSV of Soft Drop subjets;Events / 0.02"                 ,50,  0.0, 1.0);
  TH1F* h_ak8jetSDm011t1b          = new TH1F("ak8jetSDm011t1b"          ,";Mass of two hardest Soft Drop subjets (GeV);Events / 5 GeV" ,50,  0.0, 250.);  
  TH1F* h_lepPt1t1b                = new TH1F("lepPt1t1b"                ,";Lepton p_{T} (GeV);Events / 10 GeV"                         ,50,  0.0, 500.);
  TH1F* h_lepEta1t1b               = new TH1F("lepEta1t1b"               ,";Lepton #eta;Events / 0.1"                                   ,50, -2.5, 2.5);
  TH1F* h_lepEta1t1b_bb            = new TH1F("lepEta1t1b_bb"            ,";Lepton #eta;Leptons / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepEta1t1b_b             = new TH1F("lepEta1t1b_b"             ,";Lepton #eta;Leptons / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepEta1t1b_cc            = new TH1F("lepEta1t1b_cc"            ,";Lepton #eta;Leptons / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepEta1t1b_c             = new TH1F("lepEta1t1b_c"             ,";Lepton #eta;Leptons / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepEta1t1b_l             = new TH1F("lepEta1t1b_l"             ,";Lepton #eta;Leptons / 0.1"                                  ,50, -2.5, 2.5);
  TH1F* h_lepAbsEta1t1b            = new TH1F("lepAbsEta1t1b"            ,";Lepton |#eta|;Events / 0.1"                                 ,25,  0.0, 2.5);
  TH1F* h_lepSignEta1t1b           = new TH1F("lepSignEta1t1b"           ,";Lepton |#eta| * Q;Events / 0.1"                             ,50, -2.5, 2.5);
  TH1F* h_lepPhi1t1b               = new TH1F("lepPhi1t1b"               ,";Lepton #phi;Events / 0.1"                                   ,70, -3.5, 3.5);
  TH1F* h_lepBJetdR1t1b            = new TH1F("lepBJetdR1t1b"            ,";#Delta R(lepton, b-jet);Events / 0.1"                       ,35,  0.0, 3.5);
  TH1F* h_lepTJetdR1t1b            = new TH1F("lepTJetdR1t1b"            ,";#Delta R(lepton, t-jet);Events / 0.1"                       ,35,  0.0, 3.5);
  TH1F* h_lepBJetPtRel1t1b         = new TH1F("lepBJetPtRel1t1b"         ,";p_{T}^{rel}(lepton, b-jet) (GeV);Events / 2 GeV"            ,50,  0.0, 100.);
  TH1F* h_lepMiniIso1t1b           = new TH1F("lepMiniIso1t1b"           ,";lepton miniIso;Events / 0.02"                               ,50,  0.0, 1.0);

  TH2F* h_2DisoPt15                = new TH2F("2DisoPt15"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt20                = new TH2F("2DisoPt20"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt25                = new TH2F("2DisoPt25"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt30                = new TH2F("2DisoPt30"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt35                = new TH2F("2DisoPt35"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt40                = new TH2F("2DisoPt40"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);
  TH2F* h_2DisoPt45                = new TH2F("2DisoPt45"                ,";#Delta R(lepton, closest jet);p_{T}^{rel}(lepton, closest jet)",32,0.0,1.6,50,0.0,100.0);

  TH1F* h_2DisoScanPoints          = new TH1F("2DisoScanPoints"          ,";Scan points;Nevents",169,-0.5,168.5);
  TH1F* h_MiniIsoScanPoints        = new TH1F("MiniIsoScanPoints"        ,";Scan points;Nevents",31,-0.5,30.5);

  TH1F* h_elPtRaw               = new TH1F("elPtRaw"               ,";Electron p_{T} (GeV);Events / 10 GeV", 50, 0.0, 500.);
  TH1F* h_elEtaRaw              = new TH1F("elEtaRaw"              ,";Electron #eta;Events / 0.1", 50, -2.5, 2.5);
  TH1F* h_elMiniIsoRaw          = new TH1F("elMiniIsoRaw"          ,";Electron miniIso;Events / 10 GeV", 60, 0.0, 1.5);
  TH1F* h_elMiniIsoRaw_b        = new TH1F("elMiniIsoRaw_b"        ,";Electron miniIso;Events / 10 GeV", 60, 0.0, 1.5);
  TH1F* h_elMiniIsoRaw_e        = new TH1F("elMiniIsoRaw_e"        ,";Electron miniIso;Events / 10 GeV", 60, 0.0, 1.5);
  TH2F* h_elPtVsMiniIsoRaw      = new TH2F("elPtVsMiniIsoRaw"      ,";Electron p_{T} (GeV);Electron miniIso", 50, 0.0, 500., 60, 0.0, 1.5);
  TH2F* h_elEtaVsMiniIsoRaw     = new TH2F("elEtaVsMiniIsoRaw"     ,";Electron #eta;Electron miniIso", 50, -2.5, 2.5, 60, 0.0, 1.5);
  TH1F* h_muPtRaw               = new TH1F("muPtRaw"               ,";Muon p_{T} (GeV);Events / 10 GeV", 50, 0.0, 500.);
  TH1F* h_muEtaRaw              = new TH1F("muEtaRaw"              ,";Muon #eta;Events / 0.1", 50, -2.5, 2.5);
  TH1F* h_muMiniIsoRaw          = new TH1F("muMiniIsoRaw"          ,";Muon miniIso;Events / 10 GeV", 60,  0.0, 1.5);
  TH1F* h_muMiniIsoRaw_b        = new TH1F("muMiniIsoRaw_b"        ,";Muon miniIso;Events / 10 GeV", 60,  0.0, 1.5);
  TH1F* h_muMiniIsoRaw_e        = new TH1F("muMiniIsoRaw_e"        ,";Muon miniIso;Events / 10 GeV", 60,  0.0, 1.5);
  TH2F* h_muPtVsMiniIsoRaw      = new TH2F("muPtVsMiniIsoRaw"      ,";Muon p_{T} (GeV);Muon miniIso", 50, 0.0, 500., 60, 0.0, 1.5);
  TH2F* h_muEtaVsMiniIsoRaw     = new TH2F("muEtaVsMiniIsoRaw"     ,";Muon #eta;Muon miniIso", 50, -2.5, 2.5, 60, 0.0, 1.5);

  TH1F* h_lepMETdPhiRaw         = new TH1F("lepMETdPhiRaw"         ,";#Delta #phi(lepton, MET);Events / 0.1", 35, 0.0, 3.5);
  TH1F* h_ak4METdPhiRaw         = new TH1F("ak4METdPhiRaw"         ,";#Delta #phi(AK4 jet, MET);Events / 0.1", 35, 0.0, 3.5);
  TH1F* h_ak8METdPhiRaw         = new TH1F("ak8METdPhiRaw"         ,";#Delta #phi(AK8 jet, MET);Events / 0.1", 35, 0.0, 3.5);

  TH1F* h_lepPhiRaw        = new TH1F("lepPhiRaw",";lepton #phi;Events",70,-3.5,3.5);
  TH1F* h_metPhiRaw        = new TH1F("metPhiRaw",";MET #phi;Events",70,-3.5,3.5);
  TH1F* h_jetPhiRaw        = new TH1F("jetPhiRaw",";Leading AK4 jet #phi;Events",70,-3.5,3.5);
  TH1F* h_jetPtRaw         = new TH1F("jetPtRaw",";Leading AK4 jet p_{T};Events",100,50.,1050.);
  TH1F* h_metPtRaw         = new TH1F("metPtRaw",";MET p_{T};Events",50,0.,500.);

  TH2F* h_bTagSFvsPt            = new TH2F("bTagSFvsPt"            ,";b-tag SF;b-tagged jet p_{T}",50,0.9,1.0,60,0.0,600.0);
  TH2F* h_bTagSFvsCSV           = new TH2F("bTagSFvsCSV"           ,";b-tag SF;b-tagged jet CSV",50,0.9,1.0,50,0.8,1.0);

  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;
  
  int nPassSemiLep = 0;
  int nPassParton = 0;
  int nPassParticle = 0;
  
  int passStep1a = 0;
  int passStep1b = 0;
  int passStep2 = 0;
  int passStep3 = 0;
  int passStep4 = 0;
  int passStep5 = 0;
  
  int nPassCSV = 0;
  int nPassBtag = 0;
  int nPassMassCut = 0;
  int nPassTau32Cut = 0;
  int nPassTopTag = 0;
  
  int n1b = 0;
  int n1t = 0;
  int n1t1b = 0;

  int nTrueClose = 0;
  int nTrueLead = 0;
  int nLeadClose = 0;
  int nTrueLeadClose = 0;
  int nCloseTag = 0;
  int nCloseTrueTag = 0;
  int nLeadTag = 0;
  int nLeadTrueTag = 0;
  int nLeadCloseTag = 0;
  int nLeadCloseTrueTag = 0;
  
  float noIsoCount = 0.;
  float MiniIsoCounts[30] = {0.};
  float Count2DIso[7][6][4] = {0.};
  
  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {
    
    bool passParton = false;
    bool passReco = false;
    
    tree->GetEntry(i,0);
    
    if (i%10000 == 0) cout << "Event " << i << endl;
    
    if (oddOrEven != 0){
      if (i % 2 == oddOrEven - 1) continue;
    }
    // ----------------------------------------------------------------------------------------------------------------
    
    float weight = eventWeight_nom->at(0);
    if (!isData && systematic == "puUp") weight = eventWeight_puUp->at(0);
    if (!isData && systematic == "puDown") weight = eventWeight_puDown->at(0);
    if (sample.Contains("PowhegPythia8") && systematic == "Q2Up") weight *= eventWeight_Q2Up->at(0);
    if (sample.Contains("PowhegPythia8") && systematic == "Q2Down") weight *= eventWeight_Q2Down->at(0);
    if (sample.Contains("PowhegPythia8") && systematic == "PDFUp") weight *= (1.0 + eventWeight_PDF->at(0));
    if (sample.Contains("PowhegPythia8") && systematic == "PDFDown") weight *= (1.0 - eventWeight_PDF->at(0));
    if (sample.Contains("PowhegPythia8") && systematic == "ASUp") weight *= eventWeight_alphaUp->at(0);
    if (sample.Contains("PowhegPythia8") && systematic == "ASDown") weight *= eventWeight_alphaDown->at(0);

    float unfold_w_ptup = 1.0;
    float unfold_w_ptdn = 1.0;

    // Look at truth information (TTbar signal only)
    if (isSignal && sample.Contains("PowhegPythia8")) {
      
      // Do channel selection at parton level
      if ((int)truthChannel->size() == 0){
	cout << "Error: no truthChannel information!" << endl;
	continue;
      }
      if (truthChannel->at(0) == 2) {
	cout << "Error: signal event is not semileptonic at parton level!" << endl;
	continue;
      }
      if (truthChannel->at(0) == 0 && channel == "el") continue;
      if (truthChannel->at(0) == 1 && channel == "mu") continue;
      nPassSemiLep += 1;

      // Get parton-level info
      if ((int)genTopPt->size() != 0){
	nPassParton += 1;
	passParton = true;
	
	h_genTopPt->Fill(genTopPt->at(0),weight);
	h_genTopEta->Fill(genTopEta->at(0),weight);
	h_genTopPhi->Fill(genTopPhi->at(0),weight);
	h_genTTbarMass->Fill(genTTbarMass->at(0),weight);
	
	h_ptGenTop->Fill(genTopPt->at(0),weight);
	
	h_ptGenTop2->Fill(genTopPt->at(0),weight);
	h_ptGenTop3->Fill(genTopPt->at(0),weight);
	h_ptGenTop4->Fill(genTopPt->at(0),weight);
	h_ptGenTop5->Fill(genTopPt->at(0),weight);
	h_ptGenTop6->Fill(genTopPt->at(0),weight);

	float oldpt = genTopPt->at(0);
	float w_ptup = (1.0 + 0.0008*oldpt);
	float w_ptdn = (1.0 - 0.0008*oldpt);

	unfold_w_ptup = w_ptup;
	unfold_w_ptdn = w_ptdn;

	h_ptGenTopMod->Fill(oldpt,weight*w_ptup);
	h_ptGenTopMod2->Fill(oldpt,weight*w_ptup);
	h_ptGenTopMod3->Fill(oldpt,weight*w_ptup);
	h_ptGenTopMod4->Fill(oldpt,weight*w_ptup);
	h_ptGenTopMod5->Fill(oldpt,weight*w_ptup);
	h_ptGenTopMod6->Fill(oldpt,weight*w_ptup);

	h_ptGenTopModDown->Fill(oldpt,weight*w_ptdn);
	h_ptGenTopModDown2->Fill(oldpt,weight*w_ptdn);
	h_ptGenTopModDown3->Fill(oldpt,weight*w_ptdn);
	h_ptGenTopModDown4->Fill(oldpt,weight*w_ptdn);
	h_ptGenTopModDown5->Fill(oldpt,weight*w_ptdn);
	h_ptGenTopModDown6->Fill(oldpt,weight*w_ptdn);

	if (channel == "mu"){
	  if ((int)genMuPt->size() == 0) {
	    cout << "Error: no muon in mu channel!" << endl;
	    continue;
	  }
	  h_genLepPt->Fill(genMuPt->at(0),weight);
	  h_genLepEta->Fill(genMuEta->at(0),weight);
	  h_genLepPhi->Fill(genMuPhi->at(0),weight);
	}
	
	if (channel == "el"){
	  if ((int)genElPt->size() == 0){
	    cout << "Error: no electron in el channel!" << endl;
	    continue;
	  }
	  h_genLepPt->Fill(genElPt->at(0),weight);
	  h_genLepEta->Fill(genElEta->at(0),weight);
	  h_genLepPhi->Fill(genElPhi->at(0),weight);
	}
      }
      
      /*
      // TODO: implement further particle-level selection here
      //h_ngenAK4jet->Fill((int)genAK4jetPt->size(),weight);
      if ((int)genAK4jetPt->size() != 0){
        for (int it=0; it<(int)genAK4jetPt->size(); it++) {
	  //h_genAK4jetPt->Fill(genAK4jetPt->at(it),weight);
	  //h_genAK4jetEta->Fill(genAK4jetEta->at(it),weight);
	  //h_genAK4jetPhi->Fill(genAK4jetPhi->at(it),weight);
	  //h_genAK4jetMass->Fill(genAK4jetMass->at(it),weight);
	}
      }
      
      //h_ngenAK8jet->Fill((int)genAK8jetPt->size(),weight);
      if ((int)genAK8jetPt->size() != 0){
	for (int it=0; it<(int)genAK8jetPt->size(); it++) {
	  //h_genAK8jetPt->Fill(genAK8jetPt->at(it),weight);
	  //h_genAK8jetEta->Fill(genAK8jetEta->at(it),weight);
	  //h_genAK8jetPhi->Fill(genAK8jetPhi->at(it),weight);
	  //h_genAK8jetMass->Fill(genAK8jetMass->at(it),weight);
	}
      }
      
      if (channel == "mu") {
	if ((int)partMuPt->size() != 0) {
	  //h_partLepPt->Fill(partMuPt->at(0),weight);
	  //h_partLepEta->Fill(partMuEta->at(0),weight);
	  //h_partLepPhi->Fill(partMuPhi->at(0),weight);
	}
      }
      
      if (channel == "el") {
	if ((int)partElPt->size() != 0){
	  //h_partLepPt->Fill(partElPt->at(0),weight);
	  //h_partLepEta->Fill(partElEta->at(0),weight);
	  //h_partLepPhi->Fill(partElPhi->at(0),weight);
	}
      }
      */
    } // end if(isSignal)   

    // ------------------------
    // R E C O   L E V E L
    // ------------------------
    // Recall here that events are only stored in trees if they pass tree-level reco preselection, aka
    // ==1 Medium lepton with pt > 50
    // >=1 AK8 jet with pt > 350
    // >=1 AK4 jet with pt > 50

    //Construct AK4 and AK8 jet collections, with proper JEC / JER
    std::vector<TLorentzVector> ak4Jets;
    if ((int)ak4jetPt->size() == 0) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }

    // Count number of heavy flavor jets
    int nB = 0;
    int nC = 0;
    
    for (int iak4 = 0; iak4 < (int)ak4jetPt->size(); iak4++){
      TLorentzVector jetP4;
      if (systematic == "JECUp" && !isData){
	jetP4.SetPtEtaPhiM(ak4jetPt->at(iak4),ak4jetEta->at(iak4),ak4jetPhi->at(iak4),ak4jetMass->at(iak4));
	jetP4 *= (1.0 + ak4jetJECunc->at(iak4));
      }
      else if (systematic == "JECDown" && !isData){
	jetP4.SetPtEtaPhiM(ak4jetPt->at(iak4),ak4jetEta->at(iak4),ak4jetPhi->at(iak4),ak4jetMass->at(iak4));
	jetP4 *= (1.0 - ak4jetJECunc->at(iak4));
      }
      else if (systematic == "JERUp" && !isData){
      	jetP4.SetPtEtaPhiM(ak4jetPtJERup->at(iak4),ak4jetEtaJERup->at(iak4),ak4jetPhiJERup->at(iak4),ak4jetMassJERup->at(iak4));
      }
      else if (systematic == "JERDown" && !isData){
      	jetP4.SetPtEtaPhiM(ak4jetPtJERdown->at(iak4),ak4jetEtaJERdown->at(iak4),ak4jetPhiJERdown->at(iak4),ak4jetMassJERdown->at(iak4));
      }
      else {
	jetP4.SetPtEtaPhiM(ak4jetPt->at(iak4),ak4jetEta->at(iak4),ak4jetPhi->at(iak4),ak4jetMass->at(iak4));
      }
      ak4Jets.push_back(jetP4);

      if (sample.Contains("WJets")){
	if (ak4jetHadronFlavour->at(iak4) == 5) nB += 1;
	if (ak4jetHadronFlavour->at(iak4) == 4) nC += 1;
      }
    }

    std::vector<TLorentzVector> ak8Jets;
    if ((int)ak8jetPt->size() == 0) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }

    for (int iak8 = 0; iak8 < (int)ak8jetPt->size(); iak8++){
      TLorentzVector jetP4;
      if (systematic == "JECUp" && !isData){
	jetP4.SetPtEtaPhiM(ak8jetPt->at(iak8),ak8jetEta->at(iak8),ak8jetPhi->at(iak8),ak8jetMass->at(iak8));
	jetP4 *= (1.0 + ak8jetJECunc->at(iak8));
      }
      else if (systematic == "JECDown" && !isData){
	jetP4.SetPtEtaPhiM(ak8jetPt->at(iak8),ak8jetEta->at(iak8),ak8jetPhi->at(iak8),ak8jetMass->at(iak8));
	jetP4 *= (1.0 - ak8jetJECunc->at(iak8));
      }
      else if (systematic == "JERUp" && !isData){
	jetP4.SetPtEtaPhiM(ak8jetPtJERup->at(iak8),ak8jetEtaJERup->at(iak8),ak8jetPhiJERup->at(iak8),ak8jetMassJERup->at(iak8));
      }
      else if (systematic == "JERDown" && !isData){
	jetP4.SetPtEtaPhiM(ak8jetPtJERdown->at(iak8),ak8jetEtaJERdown->at(iak8),ak8jetPhiJERdown->at(iak8),ak8jetMassJERdown->at(iak8));
      }
      else {
	jetP4.SetPtEtaPhiM(ak8jetPt->at(iak8),ak8jetEta->at(iak8),ak8jetPhi->at(iak8),ak8jetMass->at(iak8));
      }
      ak8Jets.push_back(jetP4);
    }

    // Get (and count / categorize) muons
    TLorentzVector goodMu;
    float goodMuMiniIso = 0.;
    int nGoodMu = 0;
    int nMuForVeto = 0;
    
    // Loop over muons
    // We only expect 0 or 1 muons, but leaving this as a loop for now in case we change something later
    if ((int)muPt->size() != 0) {
      for (int imu = 0; imu < (int)muPt->size(); imu++){
	float dRclosest = 99.;
	float ptRelClosest = -1.0;
	if ((int)ak4Jets.size() != 0){
	  TLorentzVector muP4;
	  muP4.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),0.105);
	  for (int ij = 0; ij < (int)ak4Jets.size(); ij++){
	    float dRtemp = muP4.DeltaR(ak4Jets.at(ij));
	    if (dRtemp < dRclosest) {
	      dRclosest = dRtemp;
	      ptRelClosest = muP4.Perp(ak4Jets.at(ij).Vect());
	    }
	  }
	}
	nMuForVeto += 1; //Veto on medium muons
	
	h_muPtRaw->Fill(muPt->at(imu),weight);
	h_muEtaRaw->Fill(muEta->at(imu),weight);
	h_muMiniIsoRaw->Fill(muMiniIso->at(imu),weight);
	if (std::abs(muEta->at(imu)) < 1.479) h_muMiniIsoRaw_b->Fill(muMiniIso->at(imu),weight);
	else h_muMiniIsoRaw_e->Fill(muMiniIso->at(imu),weight);
	h_muPtVsMiniIsoRaw->Fill(muPt->at(imu),muMiniIso->at(imu),weight);
	h_muEtaVsMiniIsoRaw->Fill(muEta->at(imu),muMiniIso->at(imu),weight);
	
	// Now define 'good' muons.
	if (iso == "Loose" ||
	    (!isQCD && iso == "MiniIso10" && muMiniIso->at(imu) < 0.10 ) ||
	    (!isQCD && iso == "MiniIso20" && muMiniIso->at(imu) < 0.20 ) ||
	    (!isQCD && iso == "2DisoPt25" && (muPtRelPt25->at(imu) > 45.0 || mudRPt25->at(imu) > 0.4)) ||
	    (!isQCD && iso == "2DisoPt45" && (muPtRelPt45->at(imu) > 35.0 || mudRPt45->at(imu) > 0.4)) ||
	    (!isQCD && iso == "2DisoB2G"  && (muPtRelPt15->at(imu) > 20.0 || mudRPt15->at(imu) > 0.4)) ||
	    (!isQCD && iso == "2DisoIHNY" && (muPtRelPt25->at(imu) > 25.0 || mudRPt25->at(imu) > 0.5)) ||
	    (isQCD && iso == "MiniIso10" && muMiniIso->at(imu) > 0.10 && muMiniIso->at(imu) < 0.20) ||
	    (isQCD && iso == "MiniIso20" && muMiniIso->at(imu) > 0.20 && muMiniIso->at(imu) < 0.30)){
	  if (lepID == "Medium" || (lepID == "Tight" && muTight->at(imu))){
	    goodMu.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),0.105);
	    goodMuMiniIso = muMiniIso->at(imu);
	    nGoodMu += 1;
	  }
	}
      } //End muon loop
    }

    // Get (and count / categorize) electrons
    TLorentzVector goodEl;
    int iGoodEl = -1;
    int nGoodEl = 0; 
    int nElForVeto = 0;
    
    // Loop over electrons
    // We only expect 0 or 1 electrons, but leaving this as a loop for now in case we change something later
    if ((int)elPt->size() != 0) {
      for (int iel = 0; iel < (int)elPt->size(); iel++){
	float dRclosest = 99.;
	float ptRelClosest = -1.0;
	if ((int)ak4Jets.size() != 0){
	  TLorentzVector elP4;
	  elP4.SetPtEtaPhiM(elPt->at(iel),elEta->at(iel),elPhi->at(iel),0.0);
	  for (int ij = 0; ij < (int)ak4Jets.size(); ij++){
	    float dRtemp = elP4.DeltaR(ak4Jets.at(ij));
	    if (dRtemp < dRclosest) {
	      dRclosest = dRtemp;
	      ptRelClosest = elP4.Perp(ak4Jets.at(ij).Vect());
	    }
	  }
	}
	nElForVeto += 1;

	if (elTight->at(iel)){
	  h_elPtRaw->Fill(elPt->at(iel),weight);
	  h_elEtaRaw->Fill(elEta->at(iel),weight);
	  h_elMiniIsoRaw->Fill(elMiniIso->at(iel),weight);
	  if (std::abs(elEta->at(iel)) < 1.479) h_elMiniIsoRaw_b->Fill(elMiniIso->at(iel),weight);
	  else h_elMiniIsoRaw_e->Fill(elMiniIso->at(iel),weight);
	  h_elPtVsMiniIsoRaw->Fill(elPt->at(iel),elMiniIso->at(iel),weight);
	  h_elEtaVsMiniIsoRaw->Fill(elEta->at(iel),elMiniIso->at(iel),weight);
	}
	
	// Now define 'good' electrons
	if (iso == "Loose" ||
	    (!isQCD && iso == "MiniIso10" && elMiniIso->at(iel) < 0.10 ) ||
	    (!isQCD && iso == "MiniIso20" && elMiniIso->at(iel) < 0.20 ) ||
	    (!isQCD && iso == "2DisoPt25" && (elPtRelPt25->at(iel) > 25.0 || eldRPt25->at(iel) > 0.4)) ||
	    (!isQCD && iso == "2DisoPt45" && (elPtRelPt45->at(iel) > 25.0 || eldRPt45->at(iel) > 0.4)) ||
	    (!isQCD && iso == "2DisoB2G"  && (elPtRelPt15->at(iel) > 20.0 || eldRPt15->at(iel) > 0.4)) ||
	    (!isQCD && iso == "2DisoIHNY" && (elPtRelPt25->at(iel) > 25.0 || eldRPt25->at(iel) > 0.5)) ||
	    (isQCD && iso == "MiniIso10" && elMiniIso->at(iel) > 0.10 && elMiniIso->at(iel) < 0.20) ||
	    (isQCD && iso == "MiniIso20" && elMiniIso->at(iel) > 0.20 && elMiniIso->at(iel) < 0.30) ||
	    (isQCD && iso == "2DisoIHNY" && elPtRelPt25->at(iel) < 25.0 && eldRPt25->at(iel) < 0.5)){
	  if (lepID == "Medium" || (lepID == "Tight" && elTight->at(iel))) {
	    nGoodEl += 1;
	    iGoodEl = iel;
	  }
	}
      }
    }
    
    // Isolation studies
    if (!isData && (int)ak4Jets.size() > 1 && (int)ak8Jets.size() != 0 && ak8Jets.at(0).Perp() > MINTOPPT){//Loose jet 'preselection' -- same as regular preselection but no hemisphere requirement
      float MiniIsoCuts[30] = {0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29};
      float PtRelCuts[6] = {20.,25.,30.,35.,40.,45.};
      float dRCuts[4] = {0.4,0.45,0.4,0.55};
      
      // Make 2D plots of isolation quantities for each jet collection, for leading lepton
      if (channel == "mu" && (int)muPt->size() != 0){
	if (lepID == "Medium" || (lepID == "Tight" && muTight->at(0))){
	  h_2DisoPt15->Fill(mudRPt15->at(0),muPtRelPt15->at(0),weight);
	  h_2DisoPt20->Fill(mudRPt20->at(0),muPtRelPt20->at(0),weight);
	  h_2DisoPt25->Fill(mudRPt25->at(0),muPtRelPt25->at(0),weight);
	  h_2DisoPt30->Fill(mudRPt30->at(0),muPtRelPt30->at(0),weight);
	  h_2DisoPt35->Fill(mudRPt35->at(0),muPtRelPt35->at(0),weight);
	  h_2DisoPt40->Fill(mudRPt40->at(0),muPtRelPt40->at(0),weight);
	  h_2DisoPt45->Fill(mudRPt45->at(0),muPtRelPt45->at(0),weight);
	}
	
	noIsoCount += weight;
	
	for (int ii = 0; ii < 30; ii ++){
	  for (int jj = 0; jj < (int)muPt->size(); jj++){
	    if (lepID == "Medium" || (lepID == "Tight" && muTight->at(jj))){
	      if (muMiniIso->at(jj) < MiniIsoCuts[ii]) MiniIsoCounts[ii] += weight;
	    }
	  }
	}
	
	for (int ii = 0; ii < 6; ii++){
	  for (int jj = 0; jj < 4; jj++){
	    for (int kk = 0; kk < (int)muPt->size(); kk++){
	      if (lepID == "Medium" || (lepID == "Tight" && muTight->at(kk))){
		if (muPtRelPt15->at(kk) > PtRelCuts[ii] || mudRPt15->at(kk) > dRCuts[jj]) Count2DIso[0][ii][jj] += weight;
		if (muPtRelPt20->at(kk) > PtRelCuts[ii] || mudRPt20->at(kk) > dRCuts[jj]) Count2DIso[1][ii][jj] += weight;
		if (muPtRelPt25->at(kk) > PtRelCuts[ii] || mudRPt25->at(kk) > dRCuts[jj]) Count2DIso[2][ii][jj] += weight;
		if (muPtRelPt30->at(kk) > PtRelCuts[ii] || mudRPt30->at(kk) > dRCuts[jj]) Count2DIso[3][ii][jj] += weight;
		if (muPtRelPt35->at(kk) > PtRelCuts[ii] || mudRPt35->at(kk) > dRCuts[jj]) Count2DIso[4][ii][jj] += weight;
		if (muPtRelPt40->at(kk) > PtRelCuts[ii] || mudRPt40->at(kk) > dRCuts[jj]) Count2DIso[5][ii][jj] += weight;
		if (muPtRelPt45->at(kk) > PtRelCuts[ii] || mudRPt45->at(kk) > dRCuts[jj]) Count2DIso[6][ii][jj] += weight;
	      }
	    }
	  }
	}
      }
      
      if (channel == "el" && (int)elPtRelPt15->size() != 0){
	if (lepID == "Medium" || (lepID == "Tight" && elTight->at(0))){
	  h_2DisoPt15->Fill(eldRPt15->at(0),elPtRelPt15->at(0),weight);
	  h_2DisoPt20->Fill(eldRPt20->at(0),elPtRelPt20->at(0),weight);
	  h_2DisoPt25->Fill(eldRPt25->at(0),elPtRelPt25->at(0),weight);
	  h_2DisoPt30->Fill(eldRPt30->at(0),elPtRelPt30->at(0),weight);
	  h_2DisoPt35->Fill(eldRPt35->at(0),elPtRelPt35->at(0),weight);
	  h_2DisoPt40->Fill(eldRPt40->at(0),elPtRelPt40->at(0),weight);
	  h_2DisoPt45->Fill(eldRPt45->at(0),elPtRelPt45->at(0),weight);
	}
	
	noIsoCount += weight;
	
	for (int ii = 0; ii < 30; ii ++){
	  for (int jj = 0; jj < (int)elPt->size(); jj++){
	    if (lepID == "Medium" || (lepID == "Tight" && elTight->at(jj))){
	      if (elMiniIso->at(jj) < MiniIsoCuts[ii]) MiniIsoCounts[ii] += weight;
	    }
	  }
	}
	
	for (int ii = 0; ii < 6; ii++){
	  for (int jj = 0; jj < 4; jj++){
	    for (int kk = 0; kk < (int)elPtRelPt15->size(); kk++){
	      if (lepID == "Medium" || (lepID == "Tight" && elTight->at(kk))){
		if (elPtRelPt15->at(kk) > PtRelCuts[ii] || eldRPt15->at(kk) > dRCuts[jj]) Count2DIso[0][ii][jj] += weight;
		if (elPtRelPt20->at(kk) > PtRelCuts[ii] || eldRPt20->at(kk) > dRCuts[jj]) Count2DIso[1][ii][jj] += weight;
		if (elPtRelPt25->at(kk) > PtRelCuts[ii] || eldRPt25->at(kk) > dRCuts[jj]) Count2DIso[2][ii][jj] += weight;
		if (elPtRelPt30->at(kk) > PtRelCuts[ii] || eldRPt30->at(kk) > dRCuts[jj]) Count2DIso[3][ii][jj] += weight;
		if (elPtRelPt35->at(kk) > PtRelCuts[ii] || eldRPt35->at(kk) > dRCuts[jj]) Count2DIso[4][ii][jj] += weight;
		if (elPtRelPt40->at(kk) > PtRelCuts[ii] || eldRPt40->at(kk) > dRCuts[jj]) Count2DIso[5][ii][jj] += weight;
		if (elPtRelPt45->at(kk) > PtRelCuts[ii] || eldRPt45->at(kk) > dRCuts[jj]) Count2DIso[6][ii][jj] += weight;
	      }
	    }
	  }
	}
      }
    }

    // ------------------------
    // Preselection
    // ------------------------
    
    // Require exactly one good mu/el and veto additional el/mu
    if (iso == "Loose"){
      if (channel == "mu" && !(nMuForVeto == 1 && nElForVeto == 0)) continue; 
      if (channel == "el" && !(nElForVeto == 1 && nMuForVeto == 0)) continue;
    }    
    else {
      if (channel == "mu" && !(nGoodMu == 1 && nMuForVeto == 1 && nElForVeto == 0)) {
	if (passParton && isSignal) {
	  response.Miss(genTopPt->at(0),weight*weight_response);
	  response2.Miss(genTopPt->at(0),weight*weight_response);
	  response3.Miss(genTopPt->at(0),weight*weight_response);
	  response4.Miss(genTopPt->at(0),weight*weight_response);
	  response5.Miss(genTopPt->at(0),weight*weight_response);
	  response5fine.Miss(genTopPt->at(0),weight*weight_response);
	  response6.Miss(genTopPt->at(0),weight*weight_response);
	  response6fine.Miss(genTopPt->at(0),weight*weight_response);
	}
	continue;
      }
      if (channel == "el" && !(nGoodEl == 1 && nElForVeto == 1 && nMuForVeto == 0)) {
	if (passParton && isSignal) {
	  response.Miss(genTopPt->at(0),weight*weight_response);
	  response2.Miss(genTopPt->at(0),weight*weight_response);
	  response3.Miss(genTopPt->at(0),weight*weight_response);
	  response4.Miss(genTopPt->at(0),weight*weight_response);
	  response5.Miss(genTopPt->at(0),weight*weight_response);
	  response5fine.Miss(genTopPt->at(0),weight*weight_response);
	  response6.Miss(genTopPt->at(0),weight*weight_response);
	  response6fine.Miss(genTopPt->at(0),weight*weight_response);
	}
	continue;
      }
    }

    TLorentzVector refLep;
    float refLepMiniIso;
    float lepQ;
    if (channel == "mu") {
      refLep = goodMu;
      refLepMiniIso = goodMuMiniIso;
      lepQ = muCharge->at(0);
      if (!isData){
	if (systematic == "lepUp") weight *= getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"up",usePost);
	else if (systematic == "lepDown") weight *= getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"down",usePost);
	else weight *= getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"nom",usePost);
      }
    }
    if (channel == "el") {
      refLep.SetPtEtaPhiM(elPt->at(iGoodEl),elEta->at(iGoodEl),elPhi->at(iGoodEl),0.0);
      refLepMiniIso = elMiniIso->at(iGoodEl);
      lepQ = elCharge->at(iGoodEl);
      if (!isData) {
      	if (systematic == "lepUp") weight *= getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"up",usePost);
      	else if (systematic == "lepDown") weight *= getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"down",usePost);
	else weight *= getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"nom",usePost);
      }
    }

    h_lepMETdPhiRaw->Fill(std::abs(metPhi->at(0) - refLep.Phi()));
    h_ak4METdPhiRaw->Fill(std::abs(metPhi->at(0) - ak4Jets.at(0).Phi()));
    h_ak8METdPhiRaw->Fill(std::abs(metPhi->at(0) - ak8Jets.at(0).Phi()));

    passStep1a += 1;

    if (channel == "el"){
      h_lepPhiRaw->Fill(refLep.Phi());
      h_metPhiRaw->Fill(metPhi->at(0));
      h_jetPhiRaw->Fill(ak4Jets.at(0).Phi());
      h_jetPtRaw->Fill(ak4Jets.at(0).Perp());
      h_metPtRaw->Fill(metPt->at(0));
    }
    
    // Triangular cut, if using
    if (channel == "el" && doTriangular){
      float dphi_emet = std::abs(metPhi->at(0) - refLep.Phi());
      if (dphi_emet > 3.14159) dphi_emet = std::abs(2*3.14159 - dphi_emet);
      float dphi_jetmet = std::abs(metPhi->at(0) - ak4Jets.at(0).Phi());
      if (dphi_jetmet > 3.14159) dphi_jetmet = std::abs(2*3.14159 - dphi_jetmet);
      
      if (std::abs(dphi_emet-1.5) > 1.5 * metPt->at(0) / 75.0 || std::abs(dphi_jetmet-1.5) > 1.5 * metPt->at(0) / 75.0){
	if (passParton && isSignal) {
	  response.Miss(genTopPt->at(0),weight*weight_response);
	  response2.Miss(genTopPt->at(0),weight*weight_response);
	  response3.Miss(genTopPt->at(0),weight*weight_response);
	  response4.Miss(genTopPt->at(0),weight*weight_response);
	  response5.Miss(genTopPt->at(0),weight*weight_response);
	  response5fine.Miss(genTopPt->at(0),weight*weight_response);
	  response6.Miss(genTopPt->at(0),weight*weight_response);
	  response6fine.Miss(genTopPt->at(0),weight*weight_response);
	}
	continue;
      }
    } // End triangular cut
    
    passStep1b += 1;
    
    // Require at least two jets
    if ((int)ak4Jets.size() < 2) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }
    passStep2 += 1;
    
    // Define leptonic jets (AK4 jets with dR(jet,lep) < pi/2) and b jet candidate (leptonic jet closest to lepton)
    float bJetCandDR = 99.;
    float bJetCandPt = 0.;
    int ibJetClose = -1;
    int ibJetPt = -1;
    int nLepJet = 0;
    int nTrueBJet = 0;
    for (int it = 0; it < (int)ak4Jets.size(); it++){
      float dRtemp = refLep.DeltaR(ak4Jets.at(it));
      if (doHemiCuts && dRtemp > 3.1415 / 2.) continue;
      if (doHemiCuts && dRtemp < 0.1) continue;
      nLepJet += 1;
      if (dRtemp < bJetCandDR) {
	bJetCandDR = dRtemp;
	ibJetClose = it;
      }
      if (ak4Jets.at(it).Perp() > bJetCandPt){
	bJetCandPt = ak4Jets.at(it).Perp();
	ibJetPt = it;
      }
      if (isSignal && ak4jetHadronFlavour->at(it) == 5){
	nTrueBJet += 1;
      }
    }

    // Require at least one leptonic jet
    if (nLepJet < 1) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }
    passStep3 += 1;
      
    // Define hadronic jets (AK8 jets with dR(jet,lep) > pi/2) and top jet candidate (hadronic jet with highest pt)
    float topJetCandPt = 0.;
    int itopJetCand = -1;
    int nHadJet = 0;
    for (int it = 0; it < (int)ak8Jets.size(); it++){
      if (ak8Jets.at(it).Perp() < MINTOPPT) continue;
      float dRtemp = refLep.DeltaR(ak8Jets.at(it));
      if (doHemiCuts && dRtemp < 3.1415 / 2.) continue;
      nHadJet += 1;
      if (ak8Jets.at(it).Perp() > topJetCandPt) {
	topJetCandPt = ak8Jets.at(it).Perp();
	itopJetCand = it;
      }
    }
    
    // Require at least one hadronic jet
    if (nHadJet < 1) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }
    passStep4 += 1;
    
    // Require MET > 35 GeV
    if (metPt->at(0) < metCut) {
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
      continue;
    }
    passStep5 += 1;

    if (isSignal){
      h_nTrueLepBJet->Fill(nTrueBJet);
      if (ak4jetHadronFlavour->at(ibJetClose) == 5) nTrueClose += 1;
      if (ak4jetHadronFlavour->at(ibJetPt) == 5) nTrueLead += 1;
      if (ibJetClose == ibJetPt) nLeadClose += 1;
      if ((ibJetClose == ibJetPt) && ak4jetHadronFlavour->at(ibJetClose) == 5) nTrueLeadClose += 1;
    }

    // ----------------------------
    // Fill preselection histograms
    // ----------------------------
    if (!isData){
      if (channel == "mu"){
	if (systematic == "lepUp") h_lepSF->Fill(getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"up",usePost));
	else if (systematic == "lepDown") h_lepSF->Fill(getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"down",usePost));
	else h_lepSF->Fill(getMuonSF(refLep.Eta(),refLep.Perp(),h_muID,h_muTrig,"nom",usePost));
      }
      if (channel == "el"){
        if (systematic == "lepUp") h_lepSF->Fill(getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"up",usePost));
        else if (systematic == "lepDown") h_lepSF->Fill(getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"down",usePost));
        else h_lepSF->Fill(getElectronSF(refLep.Eta(),refLep.Perp(),h_elReco,h_elID,h_elIso,h_elTrig,"nom",usePost));
      }
      h_puWeight->Fill(eventWeight_nom->at(0));
    }
    h_metPtPre->Fill(metPt->at(0),weight);
    h_htPre->Fill(ht->at(0),weight);
    h_htLepPre->Fill(metPt->at(0)+refLep.Perp(),weight);
    h_nAK4jetPre->Fill((int)ak4Jets.size(),weight);
    h_nBjetPre->Fill(nLepJet,weight);
    h_ak4jetPtPre->Fill(ak4Jets.at(ibJetPt).Perp(),weight);
    h_ak4jetEtaPre->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
    h_ak4jetPhiPre->Fill(ak4Jets.at(ibJetPt).Phi(),weight);
    h_ak4jetMassPre->Fill(ak4Jets.at(ibJetPt).M(),weight);
    h_ak4jetCSVPre->Fill(ak4jetCSV->at(ibJetPt),weight);
    h_ak4jetVtxMassPre->Fill(ak4jetVtxMass->at(ibJetPt),weight);
    h_nAK8jetPre->Fill((int)ak8Jets.size(),weight);
    h_nTjetPre->Fill(nHadJet,weight);
    h_ak8jetPtPre->Fill(ak8Jets.at(itopJetCand).Perp(),weight);
    h_ak8jetEtaPre->Fill(ak8Jets.at(itopJetCand).Eta(),weight);
    h_ak8jetPhiPre->Fill(ak8Jets.at(itopJetCand).Phi(),weight);
    h_ak8jetYPre->Fill(ak8Jets.at(itopJetCand).Rapidity(),weight);
    h_ak8jetMassPre->Fill(ak8Jets.at(itopJetCand).M(),weight);
    h_ak8jetMassPrunedPre->Fill(ak8jetMassPruned->at(itopJetCand),weight);
    h_ak8jetMassFilteredPre->Fill(ak8jetMassFiltered->at(itopJetCand),weight);
    h_ak8jetMassTrimmedPre->Fill(ak8jetMassTrimmed->at(itopJetCand),weight);
    h_ak8jetTau1Pre->Fill(ak8jetTau1->at(itopJetCand),weight);
    h_ak8jetTau2Pre->Fill(ak8jetTau2->at(itopJetCand),weight);
    h_ak8jetTau3Pre->Fill(ak8jetTau3->at(itopJetCand),weight);
    h_ak8jetTau32Pre->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
    h_ak8jetTau21Pre->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
    h_ak8jetCSVPre->Fill(ak8jetCSV->at(itopJetCand),weight);
    h_ak8jetSDmassPre->Fill(ak8jetSDmass->at(itopJetCand),weight);
    h_ak8jetSDsubjet0ptPre->Fill(ak8jetSDsubjet0pt->at(itopJetCand),weight);
    h_ak8jetSDsubjet0massPre->Fill(ak8jetSDsubjet0mass->at(itopJetCand),weight);
    h_ak8jetSDsubjet0CSVPre->Fill(ak8jetSDsubjet0CSV->at(itopJetCand),weight);
    h_ak8jetSDsubjet1ptPre->Fill(ak8jetSDsubjet1pt->at(itopJetCand),weight);
    h_ak8jetSDsubjet1massPre->Fill(ak8jetSDsubjet1mass->at(itopJetCand),weight);
    h_ak8jetSDsubjet1CSVPre->Fill(ak8jetSDsubjet1CSV->at(itopJetCand),weight);
    if (ak8jetSDsubjet0pt->at(itopJetCand) >= 0.0 && ak8jetSDsubjet1pt->at(itopJetCand) >= 0.0){ //two valid subjets
      TLorentzVector P40;
      TLorentzVector P41;
      P40.SetPtEtaPhiM(ak8jetSDsubjet0pt->at(itopJetCand),ak8jetSDsubjet0eta->at(itopJetCand),ak8jetSDsubjet0phi->at(itopJetCand),ak8jetSDsubjet0mass->at(itopJetCand));
      P41.SetPtEtaPhiM(ak8jetSDsubjet1pt->at(itopJetCand),ak8jetSDsubjet1eta->at(itopJetCand),ak8jetSDsubjet1phi->at(itopJetCand),ak8jetSDsubjet1mass->at(itopJetCand));
      h_ak8jetSDm01Pre->Fill((P40+P41).M(),weight);
      h_ak8jetSDsubjetMaxCSVPre->Fill(max(ak8jetSDsubjet0CSV->at(itopJetCand),ak8jetSDsubjet1CSV->at(itopJetCand)),weight);
    }
    h_lepPtPre->Fill(refLep.Perp(),weight);
    h_lepEtaPre->Fill(refLep.Eta(),weight);
    h_lepAbsEtaPre->Fill(abs(refLep.Eta()),weight);
    h_lepSignEtaPre->Fill(lepQ*abs(refLep.Eta()),weight);
    h_lepPhiPre->Fill(refLep.Phi(),weight);
    h_lepBJetdRPre->Fill(refLep.DeltaR(ak4Jets.at(ibJetPt)),weight);
    h_lepTJetdRPre->Fill(refLep.DeltaR(ak8Jets.at(itopJetCand)),weight);
    h_lepBJetPtRelPre->Fill(refLep.Perp(ak4Jets.at(ibJetPt).Vect()),weight);
    h_lepMiniIsoPre->Fill(refLepMiniIso,weight);
    if (channel == "mu" && nGoodMu > 1) h_leadLepPtPre->Fill(muPt->at(0),weight);
    else if (channel == "el" && nGoodEl > 1) h_leadLepPtPre->Fill(elPt->at(0),weight);
    else h_leadLepPtPre->Fill(refLep.Perp(),weight);

    if (sample.Contains("WJets")){
      if (nB >= 2) {
	h_lepEtaPre_bb->Fill(refLep.Eta(),weight);
	h_ak4jetEtaPre_bb->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	h_ak8jetTau32Pre_bb->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	h_ak8jetTau21Pre_bb->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	h_ak8jetSDmassPre_bb->Fill(ak8jetSDmass->at(itopJetCand),weight);      
      }
      else if (nB == 1) {
	h_lepEtaPre_b->Fill(refLep.Eta(),weight);
	h_ak4jetEtaPre_b->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	h_ak8jetTau32Pre_b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	h_ak8jetTau21Pre_b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	h_ak8jetSDmassPre_b->Fill(ak8jetSDmass->at(itopJetCand),weight);      
      }
      else if (nC >= 2) {
	h_lepEtaPre_cc->Fill(refLep.Eta(),weight);
	h_ak4jetEtaPre_cc->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	h_ak8jetTau32Pre_cc->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	h_ak8jetTau21Pre_cc->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	h_ak8jetSDmassPre_cc->Fill(ak8jetSDmass->at(itopJetCand),weight);      
      }
      else if (nC == 1) {
	h_lepEtaPre_c->Fill(refLep.Eta(),weight);
	h_ak4jetEtaPre_c->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	h_ak8jetTau32Pre_c->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	h_ak8jetTau21Pre_c->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	h_ak8jetSDmassPre_c->Fill(ak8jetSDmass->at(itopJetCand),weight);      
      }
      else {
	h_lepEtaPre_l->Fill(refLep.Eta(),weight);
	h_ak4jetEtaPre_l->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	h_ak8jetTau32Pre_l->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	h_ak8jetTau21Pre_l->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	h_ak8jetSDmassPre_l->Fill(ak8jetSDmass->at(itopJetCand),weight);      
      }
    }

    // ------------------------------------------------------
    // Determine if b-jet and t-jet candidates pass tagging
    // ------------------------------------------------------

    // Using top tagging 3% WP (Loose) as described in https://indico.cern.ch/event/540674/contributions/2196235/attachments/1287602/1915958/TopTagging_2016-06-08.pdf
    float tau32cut = 0.81;
    float btagcut = 0.46;
    float lowmasscut = 105.;
    float highmasscut = 220.;
    bool passTopTag = false;
    float toptagSF = 1.0;
    float toptagUnc = 0.5;
    if (usePost) {
      if (sample.Contains("PowhegPythia8") || sample.Contains("SingleTop")){
	toptagSF = 0.96; //Posterior top-tag SF from combB fit
	toptagUnc = 0.04; //TODO: update w/ postfit values
      }
      else {
	toptagSF = 0.82; //Posterior top-mis-tag SF from combB fit
	toptagUnc = 0.05;
      }
    }
    if (ak8jetSDmass->at(itopJetCand) > lowmasscut && ak8jetSDmass->at(itopJetCand) < highmasscut){
      nPassMassCut += 1;
      if ((ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand)) < tau32cut) {
	nPassTau32Cut += 1;
	//if (max(ak8jetSDsubjet0CSV->at(itopJetCand),ak8jetSDsubjet1CSV->at(itopJetCand)) > btagcut){
	nPassTopTag += 1;
	passTopTag = true;
	if (!isData){
	  if (systematic == "TopTagUp") toptagSF *= (1.0+toptagUnc);
	  if (systematic == "TopTagDown") toptagSF *= (1.0-toptagUnc);
	}
	//}
      }
    }

    // Using medium b-tagging WP as described in https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
    float minCSV = 0.8484;
    bool passBtag = false;
    double btagSF = 1.0;

    if (isSignal){
      if (ak4jetCSV->at(ibJetClose) > minCSV) {
	//if (ak4jetVtxMass->at(ibJetClose) > 0.0){
	nCloseTag += 1;
	if (ibJetClose == ibJetPt) {
	  nLeadCloseTag += 1;
	  if (ak4jetHadronFlavour->at(ibJetClose) == 5) nLeadCloseTrueTag += 1;
	}
	if (ak4jetHadronFlavour->at(ibJetClose) == 5) nCloseTrueTag += 1;
	//}
      }
    }

    if (ak4jetCSV->at(ibJetPt) > minCSV) {
      nPassCSV += 1;
      //if (ak4jetVtxMass->at(ibJetPt) > 0.0){
      nPassBtag += 1;
      if (isSignal) nLeadTag += 1;
      passBtag = true;
      if (!isData){
	if (ak4jetHadronFlavour->at(ibJetPt) == 5){
	  if (isSignal) nLeadTrueTag += 1;
	  if (systematic == "BTagUp")        btagSF = reader.eval_auto_bounds("up",BTagEntry::FLAV_B, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else if (systematic == "BTagDown") btagSF = reader.eval_auto_bounds("down",BTagEntry::FLAV_B, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else                               btagSF = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp()); 
	}
	else if (ak4jetHadronFlavour->at(ibJetPt) == 4){
	  if (systematic == "BTagUp")        btagSF = reader.eval_auto_bounds("up",BTagEntry::FLAV_C, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else if (systematic == "BTagDown") btagSF = reader.eval_auto_bounds("down",BTagEntry::FLAV_C, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else                               btagSF = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp()); 
	}
	else {
	  if (systematic == "BTagUp")        btagSF = reader.eval_auto_bounds("up",BTagEntry::FLAV_UDSG, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else if (systematic == "BTagDown") btagSF = reader.eval_auto_bounds("down",BTagEntry::FLAV_UDSG, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp());
	  else                               btagSF = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, ak4Jets.at(ibJetPt).Eta(), ak4Jets.at(ibJetPt).Perp()); 
	}
	h_bTagSFvsPt->Fill(btagSF,ak4Jets.at(ibJetPt).Perp(),weight);
	h_bTagSFvsCSV->Fill(btagSF,ak4jetCSV->at(ibJetPt),weight);
      }
      //}
    }

    // ----------------------------
    // Fill response
    // ----------------------------

    if (passTopTag && passBtag){
      h_ptRecoTop->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop2->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop3->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop4->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop5->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop5fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop6->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);
      h_ptRecoTop6fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF);

      h_ptRecoTopMod->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod2->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod3->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod4->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod5->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod5fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod6->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopMod6fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptup);
      h_ptRecoTopModDown->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown2->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown3->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown4->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown5->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown5fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown6->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);
      h_ptRecoTopModDown6fine->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF*btagSF*unfold_w_ptdn);

      if (isSignal){
	if (passParton) {
	  response.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response2.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response3.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response4.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response5.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response5fine.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response6.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	  response6fine.Fill(ak8Jets.at(itopJetCand).Perp(),genTopPt->at(0),weight*btagSF*toptagSF*weight_response);
	}
	else {
	  response.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response2.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response3.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response4.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response5.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response5fine.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response6.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	  response6fine.Fake(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF*weight_response);
	}
      }
    }
    else if (passTopTag){
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response2.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response3.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response4.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response5.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response6.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*toptagSF*weight_response);
      }
    }
    else{
      if (passParton && isSignal) {
	response.Miss(genTopPt->at(0),weight*weight_response);
	response2.Miss(genTopPt->at(0),weight*weight_response);
	response3.Miss(genTopPt->at(0),weight*weight_response);
	response4.Miss(genTopPt->at(0),weight*weight_response);
	response5.Miss(genTopPt->at(0),weight*weight_response);
	response5fine.Miss(genTopPt->at(0),weight*weight_response);
	response6.Miss(genTopPt->at(0),weight*weight_response);
	response6fine.Miss(genTopPt->at(0),weight*weight_response);
      }
    }

    // ----------------------------
    // Fill 1b histograms
    // ----------------------------
    
    if (passBtag){
      h_metPt1b->Fill(metPt->at(0),weight*btagSF);
      h_ht1b->Fill(ht->at(0),weight*btagSF);
      h_htLep1b->Fill(metPt->at(0)+refLep.Perp(),weight*btagSF);
      h_nAK4jet1b->Fill((int)ak4Jets.size(),weight*btagSF);
      h_nBjet1b->Fill(nLepJet,weight*btagSF);
      h_ak4jetPt1b->Fill(ak4Jets.at(ibJetPt).Perp(),weight*btagSF);
      h_ak4jetEta1b->Fill(ak4Jets.at(ibJetPt).Eta(),weight*btagSF);
      h_ak4jetPhi1b->Fill(ak4Jets.at(ibJetPt).Phi(),weight*btagSF);
      h_ak4jetMass1b->Fill(ak4Jets.at(ibJetPt).M(),weight*btagSF);
      h_ak4jetCSV1b->Fill(ak4jetCSV->at(ibJetPt),weight*btagSF);
      h_ak4jetVtxMass1b->Fill(ak4jetVtxMass->at(ibJetPt),weight*btagSF);
      h_nAK8jet1b->Fill((int)ak8Jets.size(),weight*btagSF);
      h_nTjet1b->Fill(nHadJet,weight*btagSF);
      h_ak8jetPt1b->Fill(ak8Jets.at(itopJetCand).Perp(),weight*btagSF);
      h_ak8jetEta1b->Fill(ak8Jets.at(itopJetCand).Eta(),weight*btagSF);
      h_ak8jetPhi1b->Fill(ak8Jets.at(itopJetCand).Phi(),weight*btagSF);
      h_ak8jetY1b->Fill(ak8Jets.at(itopJetCand).Rapidity(),weight*btagSF);
      h_ak8jetMass1b->Fill(ak8Jets.at(itopJetCand).M(),weight*btagSF);
      h_ak8jetMassPruned1b->Fill(ak8jetMassPruned->at(itopJetCand),weight*btagSF);
      h_ak8jetMassFiltered1b->Fill(ak8jetMassFiltered->at(itopJetCand),weight*btagSF);
      h_ak8jetMassTrimmed1b->Fill(ak8jetMassTrimmed->at(itopJetCand),weight*btagSF);
      h_ak8jetTau11b->Fill(ak8jetTau1->at(itopJetCand),weight*btagSF);
      h_ak8jetTau21b->Fill(ak8jetTau2->at(itopJetCand),weight*btagSF);
      h_ak8jetTau31b->Fill(ak8jetTau3->at(itopJetCand),weight*btagSF);
      h_ak8jetTau321b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight*btagSF);
      h_ak8jetTau211b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight*btagSF);
      h_ak8jetCSV1b->Fill(ak8jetCSV->at(itopJetCand),weight*btagSF);
      h_ak8jetSDmass1b->Fill(ak8jetSDmass->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet0pt1b->Fill(ak8jetSDsubjet0pt->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet0mass1b->Fill(ak8jetSDsubjet0mass->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet0CSV1b->Fill(ak8jetSDsubjet0CSV->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet1pt1b->Fill(ak8jetSDsubjet1pt->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet1mass1b->Fill(ak8jetSDsubjet1mass->at(itopJetCand),weight*btagSF);
      h_ak8jetSDsubjet1CSV1b->Fill(ak8jetSDsubjet1CSV->at(itopJetCand),weight*btagSF);
      if (ak8jetSDsubjet0pt->at(itopJetCand) >= 0.0 && ak8jetSDsubjet1pt->at(itopJetCand) >= 0.0){ //two valid subjets
	TLorentzVector P40;
	TLorentzVector P41;
	P40.SetPtEtaPhiM(ak8jetSDsubjet0pt->at(itopJetCand),ak8jetSDsubjet0eta->at(itopJetCand),ak8jetSDsubjet0phi->at(itopJetCand),ak8jetSDsubjet0mass->at(itopJetCand));
	P41.SetPtEtaPhiM(ak8jetSDsubjet1pt->at(itopJetCand),ak8jetSDsubjet1eta->at(itopJetCand),ak8jetSDsubjet1phi->at(itopJetCand),ak8jetSDsubjet1mass->at(itopJetCand));
	h_ak8jetSDm011b->Fill((P40+P41).M(),weight*btagSF);
	h_ak8jetSDsubjetMaxCSV1b->Fill(max(ak8jetSDsubjet0CSV->at(itopJetCand),ak8jetSDsubjet1CSV->at(itopJetCand)),weight*btagSF);
      }
      h_lepPt1b->Fill(refLep.Perp(),weight*btagSF);
      h_lepEta1b->Fill(refLep.Eta(),weight*btagSF);
      h_lepAbsEta1b->Fill(abs(refLep.Eta()),weight*btagSF);
      h_lepSignEta1b->Fill(lepQ*abs(refLep.Eta()),weight*btagSF);
      h_lepPhi1b->Fill(refLep.Phi(),weight*btagSF);
      h_lepBJetdR1b->Fill(refLep.DeltaR(ak4Jets.at(ibJetPt)),weight*btagSF);
      h_lepTJetdR1b->Fill(refLep.DeltaR(ak8Jets.at(itopJetCand)),weight*btagSF);
      h_lepBJetPtRel1b->Fill(refLep.Perp(ak4Jets.at(ibJetPt).Vect()),weight*btagSF);
      h_lepMiniIso1b->Fill(refLepMiniIso,weight*btagSF);

      if (sample.Contains("WJets")){
	if (nB >= 2) {
	  h_lepEta1b_bb->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1b_bb->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321b_bb->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211b_bb->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1b_bb->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nB == 1) {
	  h_lepEta1b_b->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1b_b->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321b_b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211b_b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1b_b->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC >= 2) {
	  h_lepEta1b_cc->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1b_cc->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321b_cc->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211b_cc->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1b_cc->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC == 1) {
	  h_lepEta1b_c->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1b_c->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321b_c->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211b_c->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1b_c->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else {
	  h_lepEta1b_l->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1b_l->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321b_l->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211b_l->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1b_l->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
      }

      n1b += 1;
    }

    // ----------------------------
    // Fill 1t histograms
    // ----------------------------
    
    if (passTopTag){
      h_metPt1t->Fill(metPt->at(0),weight*toptagSF);
      h_ht1t->Fill(ht->at(0),weight*toptagSF);
      h_htLep1t->Fill(metPt->at(0)+refLep.Perp(),weight*toptagSF);
      h_nAK4jet1t->Fill((int)ak4Jets.size(),weight*toptagSF);
      h_nBjet1t->Fill(nLepJet,weight*toptagSF);
      h_ak4jetPt1t->Fill(ak4Jets.at(ibJetPt).Perp(),weight*toptagSF);
      h_ak4jetEta1t->Fill(ak4Jets.at(ibJetPt).Eta(),weight*toptagSF);
      h_ak4jetPhi1t->Fill(ak4Jets.at(ibJetPt).Phi(),weight*toptagSF);
      h_ak4jetMass1t->Fill(ak4Jets.at(ibJetPt).M(),weight*toptagSF);
      h_ak4jetCSV1t->Fill(ak4jetCSV->at(ibJetPt),weight*toptagSF);
      h_ak4jetVtxMass1t->Fill(ak4jetVtxMass->at(ibJetPt),weight*toptagSF);
      h_nAK8jet1t->Fill((int)ak8Jets.size(),weight*toptagSF);
      h_nTjet1t->Fill(nHadJet,weight*toptagSF);
      h_ak8jetPt1t->Fill(ak8Jets.at(itopJetCand).Perp(),weight*toptagSF);
      h_ak8jetEta1t->Fill(ak8Jets.at(itopJetCand).Eta(),weight*toptagSF);
      h_ak8jetPhi1t->Fill(ak8Jets.at(itopJetCand).Phi(),weight*toptagSF);
      h_ak8jetY1t->Fill(ak8Jets.at(itopJetCand).Rapidity(),weight*toptagSF);
      h_ak8jetMass1t->Fill(ak8Jets.at(itopJetCand).M(),weight*toptagSF);
      h_ak8jetMassPruned1t->Fill(ak8jetMassPruned->at(itopJetCand),weight*toptagSF);
      h_ak8jetMassFiltered1t->Fill(ak8jetMassFiltered->at(itopJetCand),weight*toptagSF);
      h_ak8jetMassTrimmed1t->Fill(ak8jetMassTrimmed->at(itopJetCand),weight*toptagSF);
      h_ak8jetTau11t->Fill(ak8jetTau1->at(itopJetCand),weight*toptagSF);
      h_ak8jetTau21t->Fill(ak8jetTau2->at(itopJetCand),weight*toptagSF);
      h_ak8jetTau31t->Fill(ak8jetTau3->at(itopJetCand),weight*toptagSF);
      h_ak8jetTau321t->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight*toptagSF);
      h_ak8jetTau211t->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight*toptagSF);
      h_ak8jetCSV1t->Fill(ak8jetCSV->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDmass1t->Fill(ak8jetSDmass->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet0pt1t->Fill(ak8jetSDsubjet0pt->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet0mass1t->Fill(ak8jetSDsubjet0mass->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet0CSV1t->Fill(ak8jetSDsubjet0CSV->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet1pt1t->Fill(ak8jetSDsubjet1pt->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet1mass1t->Fill(ak8jetSDsubjet1mass->at(itopJetCand),weight*toptagSF);
      h_ak8jetSDsubjet1CSV1t->Fill(ak8jetSDsubjet1CSV->at(itopJetCand),weight*toptagSF);
      if (ak8jetSDsubjet0pt->at(itopJetCand) >= 0.0 && ak8jetSDsubjet1pt->at(itopJetCand) >= 0.0){ //two valid subjets
	TLorentzVector P40;
	TLorentzVector P41;
	P40.SetPtEtaPhiM(ak8jetSDsubjet0pt->at(itopJetCand),ak8jetSDsubjet0eta->at(itopJetCand),ak8jetSDsubjet0phi->at(itopJetCand),ak8jetSDsubjet0mass->at(itopJetCand));
	P41.SetPtEtaPhiM(ak8jetSDsubjet1pt->at(itopJetCand),ak8jetSDsubjet1eta->at(itopJetCand),ak8jetSDsubjet1phi->at(itopJetCand),ak8jetSDsubjet1mass->at(itopJetCand));
	h_ak8jetSDm011t->Fill((P40+P41).M(),weight*toptagSF);
	h_ak8jetSDsubjetMaxCSV1t->Fill(max(ak8jetSDsubjet0CSV->at(itopJetCand),ak8jetSDsubjet1CSV->at(itopJetCand)),weight*toptagSF);
      }
      h_lepPt1t->Fill(refLep.Perp(),weight*toptagSF);
      h_lepEta1t->Fill(refLep.Eta(),weight*toptagSF);
      h_lepAbsEta1t->Fill(abs(refLep.Eta()),weight*toptagSF);
      h_lepSignEta1t->Fill(lepQ*abs(refLep.Eta()),weight*toptagSF);
      h_lepPhi1t->Fill(refLep.Phi(),weight*toptagSF);
      h_lepBJetdR1t->Fill(refLep.DeltaR(ak4Jets.at(ibJetPt)),weight*toptagSF);
      h_lepTJetdR1t->Fill(refLep.DeltaR(ak8Jets.at(itopJetCand)),weight*toptagSF);
      h_lepBJetPtRel1t->Fill(refLep.Perp(ak4Jets.at(ibJetPt).Vect()),weight*toptagSF);
      h_lepMiniIso1t->Fill(refLepMiniIso,weight*toptagSF);

      if (sample.Contains("WJets")){
	if (nB >= 2) {
	  h_lepEta1t_bb->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t_bb->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t_bb->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t_bb->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t_bb->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nB == 1) {
	  h_lepEta1t_b->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t_b->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t_b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t_b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t_b->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC >= 2) {
	  h_lepEta1t_cc->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t_cc->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t_cc->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t_cc->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t_cc->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC == 1) {
	  h_lepEta1t_c->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t_c->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t_c->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t_c->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t_c->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else {
	  h_lepEta1t_l->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t_l->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t_l->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t_l->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t_l->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
      }
      n1t += 1;
    }

    // ----------------------------
    // Fill 1t1b histograms
    // ----------------------------
    
    if (passTopTag && passBtag){
      h_metPt1t1b->Fill(metPt->at(0),weight*btagSF*toptagSF);
      h_ht1t1b->Fill(ht->at(0),weight*btagSF*toptagSF);
      h_htLep1t1b->Fill(metPt->at(0)+refLep.Perp(),weight*btagSF*toptagSF);
      h_nAK4jet1t1b->Fill((int)ak4Jets.size(),weight*btagSF*toptagSF);
      h_nBjet1t1b->Fill(nLepJet,weight*btagSF*toptagSF);
      h_ak4jetPt1t1b->Fill(ak4Jets.at(ibJetPt).Perp(),weight*btagSF*toptagSF);
      h_ak4jetEta1t1b->Fill(ak4Jets.at(ibJetPt).Eta(),weight*btagSF*toptagSF);
      h_ak4jetPhi1t1b->Fill(ak4Jets.at(ibJetPt).Phi(),weight*btagSF*toptagSF);
      h_ak4jetMass1t1b->Fill(ak4Jets.at(ibJetPt).M(),weight*btagSF*toptagSF);
      h_ak4jetCSV1t1b->Fill(ak4jetCSV->at(ibJetPt),weight*btagSF*toptagSF);
      h_ak4jetVtxMass1t1b->Fill(ak4jetVtxMass->at(ibJetPt),weight*btagSF*toptagSF);
      h_nAK8jet1t1b->Fill((int)ak8Jets.size(),weight*btagSF*toptagSF);
      h_nTjet1t1b->Fill(nHadJet,weight*btagSF*toptagSF);
      h_ak8jetPt1t1b->Fill(ak8Jets.at(itopJetCand).Perp(),weight*btagSF*toptagSF);
      h_ak8jetEta1t1b->Fill(ak8Jets.at(itopJetCand).Eta(),weight*btagSF*toptagSF);
      h_ak8jetPhi1t1b->Fill(ak8Jets.at(itopJetCand).Phi(),weight*btagSF*toptagSF);
      h_ak8jetY1t1b->Fill(ak8Jets.at(itopJetCand).Rapidity(),weight*btagSF*toptagSF);
      h_ak8jetMass1t1b->Fill(ak8Jets.at(itopJetCand).M(),weight*btagSF*toptagSF);
      h_ak8jetMassPruned1t1b->Fill(ak8jetMassPruned->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetMassFiltered1t1b->Fill(ak8jetMassFiltered->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetMassTrimmed1t1b->Fill(ak8jetMassTrimmed->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetTau11t1b->Fill(ak8jetTau1->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetTau21t1b->Fill(ak8jetTau2->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetTau31t1b->Fill(ak8jetTau3->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetTau321t1b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetTau211t1b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetCSV1t1b->Fill(ak8jetCSV->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDmass1t1b->Fill(ak8jetSDmass->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet0pt1t1b->Fill(ak8jetSDsubjet0pt->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet0mass1t1b->Fill(ak8jetSDsubjet0mass->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet0CSV1t1b->Fill(ak8jetSDsubjet0CSV->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet1pt1t1b->Fill(ak8jetSDsubjet1pt->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet1mass1t1b->Fill(ak8jetSDsubjet1mass->at(itopJetCand),weight*btagSF*toptagSF);
      h_ak8jetSDsubjet1CSV1t1b->Fill(ak8jetSDsubjet1CSV->at(itopJetCand),weight*btagSF*toptagSF);
      if (ak8jetSDsubjet0pt->at(itopJetCand) >= 0.0 && ak8jetSDsubjet1pt->at(itopJetCand) >= 0.0){ //two valid subjets
	TLorentzVector P40;
	TLorentzVector P41;
	P40.SetPtEtaPhiM(ak8jetSDsubjet0pt->at(itopJetCand),ak8jetSDsubjet0eta->at(itopJetCand),ak8jetSDsubjet0phi->at(itopJetCand),ak8jetSDsubjet0mass->at(itopJetCand));
	P41.SetPtEtaPhiM(ak8jetSDsubjet1pt->at(itopJetCand),ak8jetSDsubjet1eta->at(itopJetCand),ak8jetSDsubjet1phi->at(itopJetCand),ak8jetSDsubjet1mass->at(itopJetCand));
	h_ak8jetSDm011t1b->Fill((P40+P41).M(),weight*btagSF*toptagSF);
	h_ak8jetSDsubjetMaxCSV1t1b->Fill(max(ak8jetSDsubjet0CSV->at(itopJetCand),ak8jetSDsubjet1CSV->at(itopJetCand)),weight*btagSF*toptagSF);
      }
      h_lepPt1t1b->Fill(refLep.Perp(),weight*btagSF*toptagSF);
      h_lepEta1t1b->Fill(refLep.Eta(),weight*btagSF*toptagSF);
      h_lepAbsEta1t1b->Fill(abs(refLep.Eta()),weight*btagSF*toptagSF);
      h_lepSignEta1t1b->Fill(lepQ*abs(refLep.Eta()),weight*btagSF*toptagSF);
      h_lepPhi1t1b->Fill(refLep.Phi(),weight*btagSF*toptagSF);
      h_lepBJetdR1t1b->Fill(refLep.DeltaR(ak4Jets.at(ibJetPt)),weight*btagSF*toptagSF);
      h_lepTJetdR1t1b->Fill(refLep.DeltaR(ak8Jets.at(itopJetCand)),weight*btagSF*toptagSF);
      h_lepBJetPtRel1t1b->Fill(refLep.Perp(ak4Jets.at(ibJetPt).Vect()),weight*btagSF*toptagSF);
      h_lepMiniIso1t1b->Fill(refLepMiniIso,weight*btagSF*toptagSF);

      if (sample.Contains("WJets")){
	if (nB >= 2) {
	  h_lepEta1t1b_bb->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t1b_bb->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t1b_bb->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t1b_bb->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t1b_bb->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nB == 1) {
	  h_lepEta1t1b_b->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t1b_b->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t1b_b->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t1b_b->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t1b_b->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC >= 2) {
	  h_lepEta1t1b_cc->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t1b_cc->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t1b_cc->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t1b_cc->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t1b_cc->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else if (nC == 1) {
	  h_lepEta1t1b_c->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t1b_c->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t1b_c->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t1b_c->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t1b_c->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
	else {
	  h_lepEta1t1b_l->Fill(refLep.Eta(),weight);
	  h_ak4jetEta1t1b_l->Fill(ak4Jets.at(ibJetPt).Eta(),weight);
	  h_ak8jetTau321t1b_l->Fill(ak8jetTau3->at(itopJetCand) / ak8jetTau2->at(itopJetCand),weight);
	  h_ak8jetTau211t1b_l->Fill(ak8jetTau2->at(itopJetCand) / ak8jetTau1->at(itopJetCand),weight);
	  h_ak8jetSDmass1t1b_l->Fill(ak8jetSDmass->at(itopJetCand),weight);      
	}
      }
      n1t1b += 1;
    }

  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------

  tree->Delete();

  // Store isolation information needed to produce ROC curves
  h_2DisoScanPoints->SetBinContent(1,noIsoCount);
  for (int ii = 0; ii < 7; ii++){
    for (int jj = 0; jj < 6; jj++){
      for (int kk = 0; kk < 4; kk++){
	h_2DisoScanPoints->SetBinContent(2+ii*6*4+jj*4+kk,Count2DIso[ii][jj][kk]);
      }
    }
  }
  h_MiniIsoScanPoints->SetBinContent(1,noIsoCount);
  for (int ii = 0; ii < 30; ii++){
    h_MiniIsoScanPoints->SetBinContent(ii+2,MiniIsoCounts[ii]);
  }


  cout << "Finished looping over events!" << endl;
  cout << "For sample " << sample << " : " << endl;
  cout << nevt      << " events in tree" << endl;
  if (isSignal) cout << nPassSemiLep << " events in channel" << endl;
  cout << endl;
  cout << passStep1a << " events with ==1 lepton" << endl;
  cout << passStep1b << " events passing Triangular cut" << endl;
  cout << passStep2 << " events with >=2 jets" << endl;
  cout << passStep3 << " events with >=1 leptonic jet" << endl;
  cout << passStep4 << " events with >=1 hadronic jet" << endl;
  cout << passStep5 << " events pass MET cut" << endl;
  cout << endl;
  cout << n1b       << " events in 1b region" << endl;
  cout << n1t       << " events in 1t region" << endl;
  cout << n1t1b     << " events in 1t1b region" << endl;
  cout << endl;
  cout << "B-tagging stats:" << endl;
  cout << passStep5 << " b-jet candidates" << endl;
  cout << nPassCSV  << " passing CSV cut" << endl;
  //cout << nPassBtag << " passing VtxMass cut" << endl;
  cout << endl;
  cout << "Top-tagging stats:" << endl;
  cout << passStep5 << " top-jet candidates" << endl;
  cout << nPassMassCut << " pass SD mass cut" << endl;
  cout << nPassTau32Cut << " pass Tau32 cut" << endl;
  cout << nPassTopTag << " pass max CSV cut" << endl;
  cout << endl;
  if (isSignal){
    cout << "B-tagging troubleshoot: " << endl;
    cout << "Number of events w/ close jet true b: " << nTrueClose << endl;
    cout << "Number of events w/ lead jet true b: " << nTrueLead << endl;
    cout << "Number of events w/ close jet and lead jet same: " << nLeadClose << endl;
    cout << "Number of events w/ close jet and lead jet same, and true b: " << nTrueLeadClose << endl;
    cout << "Number of b-tagged events w/ close jet: " << nCloseTag << endl;
    cout << "Number of b-tagged events w/ close jet true b: " << nCloseTrueTag << endl;
    cout << "Number of b-tagged events w/ lead jet: " << nLeadTag << endl;
    cout << "Number of b-tagged events w/ lead jet true b: " << nLeadTrueTag << endl;
    cout << "Number of b-tagged events w/ lead jet and close jet same: " << nLeadCloseTag << endl;
    cout << "Number of b-tagged events w/ lead jet and close jet same, true b: " << nLeadCloseTrueTag << endl;
  }
    
  // -------------------------------------------------------------------------------------------
  // output file for histograms
  // -------------------------------------------------------------------------------------------

  TString outname = sample + "_" + channel;
  if (sample.Contains("Data")) outname = sample;
  if (!isData) outname = outname + "_" + systematic;
  if (isQCD) outname = outname + "_qcd";
  if (oddOrEven == 1) outname = outname + "_odd";
  if (oddOrEven == 2) outname = outname + "_even";
  if (!isData && usePost) outname = outname+"_post";
  TFile* fout = new TFile(OUTDIR + "/hists_"+outname+".root","recreate");

  // -------------------------------------------------------------------------------------------
  // * * * * * * D R A W   A N D   S A V E   P L O T S * * * * * * 
  // -------------------------------------------------------------------------------------------

  if (!isData){
    h_lepSF->Write();
    h_puWeight->Write();
  }

  if (isSignal){
    h_nTrueLepBJet->Write();

    h_genTopPt->Write();
    h_genTopEta->Write();
    h_genTopPhi->Write();
    h_genTTbarMass->Write();
    h_genLepPt->Write();
    h_genLepEta->Write();
    h_genLepPhi->Write();
    
    h_ptGenTop->Write();
    response.Write();

    h_ptGenTop2->Write();
    h_ptGenTop3->Write();
    h_ptGenTop4->Write();
    h_ptGenTop5->Write();
    h_ptGenTop6->Write();

    h_ptGenTopMod->Write();
    h_ptGenTopMod2->Write();
    h_ptGenTopMod3->Write();
    h_ptGenTopMod4->Write();
    h_ptGenTopMod5->Write();
    h_ptGenTopMod6->Write();
    h_ptGenTopModDown->Write();
    h_ptGenTopModDown2->Write();
    h_ptGenTopModDown3->Write();
    h_ptGenTopModDown4->Write();
    h_ptGenTopModDown5->Write();
    h_ptGenTopModDown6->Write();

    response2.Write();
    response3.Write();
    response4.Write();
    response5.Write();
    response5fine.Write();
    response6.Write();
    response6fine.Write();


    /*
    h_genAK4jetPt->Write();
    h_genAK4jetEta->Write();
    h_genAK4jetPhi->Write();
    h_genAK4jetMass->Write();
    h_ngenAK4jet->Write();
    h_genAK8jetPt->Write();
    h_genAK8jetEta->Write();
    h_genAK8jetPhi->Write();
    h_genAK8jetMass->Write();
    h_ngenAK8jet->Write();
    h_partLepPt->Write();
    h_partLepEta->Write();
    h_partLepPhi->Write();
    */
  }

  h_ptRecoTop->Write();
  h_ptRecoTop2->Write();
  h_ptRecoTop3->Write();
  h_ptRecoTop4->Write();
  h_ptRecoTop5->Write();
  h_ptRecoTop5fine->Write();
  h_ptRecoTop6->Write();
  h_ptRecoTop6fine->Write();

  h_ptRecoTopMod->Write();
  h_ptRecoTopMod2->Write();
  h_ptRecoTopMod3->Write();
  h_ptRecoTopMod4->Write();
  h_ptRecoTopMod5->Write();
  h_ptRecoTopMod5fine->Write();
  h_ptRecoTopMod6->Write();
  h_ptRecoTopMod6fine->Write();
  h_ptRecoTopModDown->Write();
  h_ptRecoTopModDown2->Write();
  h_ptRecoTopModDown3->Write();
  h_ptRecoTopModDown4->Write();
  h_ptRecoTopModDown5->Write();
  h_ptRecoTopModDown5fine->Write();
  h_ptRecoTopModDown6->Write();
  h_ptRecoTopModDown6fine->Write();

  h_metPtPre->Write();
  h_htPre->Write();
  h_htLepPre->Write();
  h_nAK4jetPre->Write();
  h_nBjetPre->Write();
  h_ak4jetPtPre->Write();
  h_ak4jetEtaPre->Write();
  h_ak4jetPhiPre->Write();
  h_ak4jetMassPre->Write();
  h_ak4jetCSVPre->Write();
  h_ak4jetVtxMassPre->Write();
  h_nAK8jetPre->Write();
  h_nTjetPre->Write();
  h_ak8jetPtPre->Write();
  h_ak8jetEtaPre->Write();
  h_ak8jetPhiPre->Write();
  h_ak8jetYPre->Write();
  h_ak8jetMassPre->Write();
  h_ak8jetMassPrunedPre->Write();
  h_ak8jetMassFilteredPre->Write();
  h_ak8jetMassTrimmedPre->Write();
  h_ak8jetTau1Pre->Write();
  h_ak8jetTau2Pre->Write();
  h_ak8jetTau3Pre->Write();
  h_ak8jetTau32Pre->Write();
  h_ak8jetTau21Pre->Write();
  h_ak8jetCSVPre->Write();
  h_ak8jetSDmassPre->Write();
  h_ak8jetSDsubjet0ptPre->Write();
  h_ak8jetSDsubjet0massPre->Write();
  h_ak8jetSDsubjet0CSVPre->Write();
  h_ak8jetSDsubjet1ptPre->Write();
  h_ak8jetSDsubjet1massPre->Write();
  h_ak8jetSDsubjet1CSVPre->Write();
  h_ak8jetSDsubjetMaxCSVPre->Write();
  h_ak8jetSDm01Pre->Write();
  h_lepPtPre->Write();
  h_lepEtaPre->Write();
  h_lepAbsEtaPre->Write();
  h_lepSignEtaPre->Write();
  h_lepPhiPre->Write();
  h_lepBJetdRPre->Write();
  h_lepTJetdRPre->Write();
  h_lepBJetPtRelPre->Write();
  h_lepMiniIsoPre->Write();
  h_leadLepPtPre->Write();

  if (sample.Contains("WJets")){
    h_ak4jetEtaPre_bb->Write();
    h_ak4jetEtaPre_b->Write();
    h_ak4jetEtaPre_cc->Write();
    h_ak4jetEtaPre_c->Write();
    h_ak4jetEtaPre_l->Write();
    h_ak8jetTau32Pre_bb->Write();
    h_ak8jetTau32Pre_b->Write();
    h_ak8jetTau32Pre_cc->Write();
    h_ak8jetTau32Pre_c->Write();
    h_ak8jetTau32Pre_l->Write();
    h_ak8jetTau21Pre_bb->Write();
    h_ak8jetTau21Pre_b->Write();
    h_ak8jetTau21Pre_cc->Write();
    h_ak8jetTau21Pre_c->Write();
    h_ak8jetTau21Pre_l->Write();
    h_ak8jetSDmassPre_bb->Write();
    h_ak8jetSDmassPre_b->Write();
    h_ak8jetSDmassPre_cc->Write();
    h_ak8jetSDmassPre_c->Write();
    h_ak8jetSDmassPre_l->Write();
    h_lepEtaPre_bb->Write();
    h_lepEtaPre_b->Write();
    h_lepEtaPre_cc->Write();
    h_lepEtaPre_c->Write();
    h_lepEtaPre_l->Write();
  }


  // ------------
  // 1b
  // ------------

  h_metPt1b->Write();
  h_ht1b->Write();
  h_htLep1b->Write();
  h_nAK4jet1b->Write();
  h_nBjet1b->Write();
  h_ak4jetPt1b->Write();
  h_ak4jetEta1b->Write();
  h_ak4jetPhi1b->Write();
  h_ak4jetMass1b->Write();
  h_ak4jetCSV1b->Write();
  h_ak4jetVtxMass1b->Write();
  h_nAK8jet1b->Write();
  h_nTjet1b->Write();
  h_ak8jetPt1b->Write();
  h_ak8jetEta1b->Write();
  h_ak8jetPhi1b->Write();
  h_ak8jetY1b->Write();
  h_ak8jetMass1b->Write();
  h_ak8jetMassPruned1b->Write();
  h_ak8jetMassFiltered1b->Write();
  h_ak8jetMassTrimmed1b->Write();
  h_ak8jetTau11b->Write();
  h_ak8jetTau21b->Write();
  h_ak8jetTau31b->Write();
  h_ak8jetTau321b->Write();
  h_ak8jetTau211b->Write();
  h_ak8jetCSV1b->Write();
  h_ak8jetSDmass1b->Write();
  h_ak8jetSDsubjet0pt1b->Write();
  h_ak8jetSDsubjet0mass1b->Write();
  h_ak8jetSDsubjet0CSV1b->Write();
  h_ak8jetSDsubjet1pt1b->Write();
  h_ak8jetSDsubjet1mass1b->Write();
  h_ak8jetSDsubjet1CSV1b->Write();
  h_ak8jetSDsubjetMaxCSV1b->Write();
  h_ak8jetSDm011b->Write();
  h_lepPt1b->Write();
  h_lepEta1b->Write();
  h_lepAbsEta1b->Write();
  h_lepSignEta1b->Write();
  h_lepPhi1b->Write();
  h_lepBJetdR1b->Write();
  h_lepTJetdR1b->Write();
  h_lepBJetPtRel1b->Write();
  h_lepMiniIso1b->Write();

  if (sample.Contains("WJets")){
    h_ak4jetEta1b_bb->Write();
    h_ak4jetEta1b_b->Write();
    h_ak4jetEta1b_cc->Write();
    h_ak4jetEta1b_c->Write();
    h_ak4jetEta1b_l->Write();
    h_ak8jetTau321b_bb->Write();
    h_ak8jetTau321b_b->Write();
    h_ak8jetTau321b_cc->Write();
    h_ak8jetTau321b_c->Write();
    h_ak8jetTau321b_l->Write();
    h_ak8jetTau211b_bb->Write();
    h_ak8jetTau211b_b->Write();
    h_ak8jetTau211b_cc->Write();
    h_ak8jetTau211b_c->Write();
    h_ak8jetTau211b_l->Write();
    h_ak8jetSDmass1b_bb->Write();
    h_ak8jetSDmass1b_b->Write();
    h_ak8jetSDmass1b_cc->Write();
    h_ak8jetSDmass1b_c->Write();
    h_ak8jetSDmass1b_l->Write();
    h_lepEta1b_bb->Write();
    h_lepEta1b_b->Write();
    h_lepEta1b_cc->Write();
    h_lepEta1b_c->Write();
    h_lepEta1b_l->Write();
  }
    
  // ------------
  // 1t
  // ------------
  
  h_metPt1t->Write();
  h_ht1t->Write();
  h_htLep1t->Write();
  h_nAK4jet1t->Write();
  h_nBjet1t->Write();
  h_ak4jetPt1t->Write();
  h_ak4jetEta1t->Write();
  h_ak4jetPhi1t->Write();
  h_ak4jetMass1t->Write();
  h_ak4jetCSV1t->Write();
  h_ak4jetVtxMass1t->Write();
  h_nAK8jet1t->Write();
  h_nTjet1t->Write();
  h_ak8jetPt1t->Write();
  h_ak8jetEta1t->Write();
  h_ak8jetPhi1t->Write();
  h_ak8jetY1t->Write();
  h_ak8jetMass1t->Write();
  h_ak8jetMassPruned1t->Write();
  h_ak8jetMassFiltered1t->Write();
  h_ak8jetMassTrimmed1t->Write();
  h_ak8jetTau11t->Write();
  h_ak8jetTau21t->Write();
  h_ak8jetTau31t->Write();
  h_ak8jetTau321t->Write();
  h_ak8jetTau211t->Write();
  h_ak8jetCSV1t->Write();
  h_ak8jetSDmass1t->Write();
  h_ak8jetSDsubjet0pt1t->Write();
  h_ak8jetSDsubjet0mass1t->Write();
  h_ak8jetSDsubjet0CSV1t->Write();
  h_ak8jetSDsubjet1pt1t->Write();
  h_ak8jetSDsubjet1mass1t->Write();
  h_ak8jetSDsubjet1CSV1t->Write();
  h_ak8jetSDsubjetMaxCSV1t->Write();
  h_ak8jetSDm011t->Write();
  h_lepPt1t->Write();
  h_lepEta1t->Write();
  h_lepAbsEta1t->Write();
  h_lepSignEta1t->Write();
  h_lepPhi1t->Write();
  h_lepBJetdR1t->Write();
  h_lepTJetdR1t->Write();
  h_lepBJetPtRel1t->Write();
  h_lepMiniIso1t->Write();

  if (sample.Contains("WJets")){
    h_ak4jetEta1t_bb->Write();
    h_ak4jetEta1t_b->Write();
    h_ak4jetEta1t_cc->Write();
    h_ak4jetEta1t_c->Write();
    h_ak4jetEta1t_l->Write();
    h_ak8jetTau321t_bb->Write();
    h_ak8jetTau321t_b->Write();
    h_ak8jetTau321t_cc->Write();
    h_ak8jetTau321t_c->Write();
    h_ak8jetTau321t_l->Write();
    h_ak8jetTau211t_bb->Write();
    h_ak8jetTau211t_b->Write();
    h_ak8jetTau211t_cc->Write();
    h_ak8jetTau211t_c->Write();
    h_ak8jetTau211t_l->Write();
    h_ak8jetSDmass1t_bb->Write();
    h_ak8jetSDmass1t_b->Write();
    h_ak8jetSDmass1t_cc->Write();
    h_ak8jetSDmass1t_c->Write();
    h_ak8jetSDmass1t_l->Write();
    h_lepEta1t_bb->Write();
    h_lepEta1t_b->Write();
    h_lepEta1t_cc->Write();
    h_lepEta1t_c->Write();
    h_lepEta1t_l->Write();
  }

  /////////////
  // 1t1b
  /////////////
  
  h_metPt1t1b->Write();
  h_ht1t1b->Write();
  h_htLep1t1b->Write();
  h_nAK4jet1t1b->Write();
  h_nBjet1t1b->Write();
  h_ak4jetPt1t1b->Write();
  h_ak4jetEta1t1b->Write();
  h_ak4jetPhi1t1b->Write();
  h_ak4jetMass1t1b->Write();
  h_ak4jetCSV1t1b->Write();
  h_ak4jetVtxMass1t1b->Write();
  h_nAK8jet1t1b->Write();
  h_nTjet1t1b->Write();
  h_ak8jetPt1t1b->Write();
  h_ak8jetEta1t1b->Write();
  h_ak8jetPhi1t1b->Write();
  h_ak8jetY1t1b->Write();
  h_ak8jetMass1t1b->Write();
  h_ak8jetMassPruned1t1b->Write();
  h_ak8jetMassFiltered1t1b->Write();
  h_ak8jetMassTrimmed1t1b->Write();
  h_ak8jetTau11t1b->Write();
  h_ak8jetTau21t1b->Write();
  h_ak8jetTau31t1b->Write();
  h_ak8jetTau321t1b->Write();
  h_ak8jetTau211t1b->Write();
  h_ak8jetCSV1t1b->Write();
  h_ak8jetSDmass1t1b->Write();
  h_ak8jetSDsubjet0pt1t1b->Write();
  h_ak8jetSDsubjet0mass1t1b->Write();
  h_ak8jetSDsubjet0CSV1t1b->Write();
  h_ak8jetSDsubjet1pt1t1b->Write();
  h_ak8jetSDsubjet1mass1t1b->Write();
  h_ak8jetSDsubjet1CSV1t1b->Write();
  h_ak8jetSDsubjetMaxCSV1t1b->Write();
  h_ak8jetSDm011t1b->Write();
  h_lepPt1t1b->Write();
  h_lepEta1t1b->Write();
  h_lepAbsEta1t1b->Write();
  h_lepSignEta1t1b->Write();
  h_lepPhi1t1b->Write();
  h_lepBJetdR1t1b->Write();
  h_lepTJetdR1t1b->Write();
  h_lepBJetPtRel1t1b->Write();
  h_lepMiniIso1t1b->Write();

  if (sample.Contains("WJets")){
    h_ak4jetEta1t1b_bb->Write();
    h_ak4jetEta1t1b_b->Write();
    h_ak4jetEta1t1b_cc->Write();
    h_ak4jetEta1t1b_c->Write();
    h_ak4jetEta1t1b_l->Write();
    h_ak8jetTau321t1b_bb->Write();
    h_ak8jetTau321t1b_b->Write();
    h_ak8jetTau321t1b_cc->Write();
    h_ak8jetTau321t1b_c->Write();
    h_ak8jetTau321t1b_l->Write();
    h_ak8jetTau211t1b_bb->Write();
    h_ak8jetTau211t1b_b->Write();
    h_ak8jetTau211t1b_cc->Write();
    h_ak8jetTau211t1b_c->Write();
    h_ak8jetTau211t1b_l->Write();
    h_ak8jetSDmass1t1b_bb->Write();
    h_ak8jetSDmass1t1b_b->Write();
    h_ak8jetSDmass1t1b_cc->Write();
    h_ak8jetSDmass1t1b_c->Write();
    h_ak8jetSDmass1t1b_l->Write();
    h_lepEta1t1b_bb->Write();
    h_lepEta1t1b_b->Write();
    h_lepEta1t1b_cc->Write();
    h_lepEta1t1b_c->Write();
    h_lepEta1t1b_l->Write();
  }

  h_2DisoPt15->Write();
  h_2DisoPt20->Write();
  h_2DisoPt25->Write();
  h_2DisoPt30->Write();
  h_2DisoPt35->Write();
  h_2DisoPt40->Write();
  h_2DisoPt45->Write();
  
  h_2DisoScanPoints->Write();
  h_MiniIsoScanPoints->Write();

  h_elPtRaw->Write();
  h_elEtaRaw->Write();
  h_elMiniIsoRaw->Write();
  h_elMiniIsoRaw_b->Write();
  h_elMiniIsoRaw_e->Write();
  h_elPtVsMiniIsoRaw->Write();
  h_elEtaVsMiniIsoRaw->Write();
  h_muPtRaw->Write();
  h_muEtaRaw->Write();
  h_muMiniIsoRaw->Write();
  h_muMiniIsoRaw_b->Write();
  h_muMiniIsoRaw_e->Write();
  h_muPtVsMiniIsoRaw->Write();
  h_muEtaVsMiniIsoRaw->Write();

  h_lepMETdPhiRaw->Write();
  h_ak4METdPhiRaw->Write();
  h_ak8METdPhiRaw->Write();

  if (!isData) {
    h_bTagSFvsPt->Write();
    h_bTagSFvsCSV->Write();
  }

  h_lepPhiRaw->Write();
  h_metPhiRaw->Write();
  h_jetPhiRaw->Write();
  h_jetPtRaw->Write();
  h_metPtRaw->Write();
  
  fout->Close();
  delete fout;

  // Do cleanup (messy, but can't see another option currently)
  h_nTrueLepBJet->Delete();
  
  h_lepSF->Delete();
  h_puWeight->Delete();

  h_genTopPt->Delete();
  h_genTopEta->Delete();
  h_genTopPhi->Delete();
  h_genLepPt->Delete();
  h_genLepEta->Delete();
  h_genLepPhi->Delete();
  h_genTTbarMass->Delete();

  h_ptGenTop->Delete();
  h_ptRecoTop->Delete();
  response.Delete();

  h_ptGenTop2->Delete();
  h_ptGenTop3->Delete();
  h_ptGenTop4->Delete();
  h_ptGenTop5->Delete();
  h_ptGenTop6->Delete();

  h_ptRecoTop2->Delete();
  h_ptRecoTop3->Delete();
  h_ptRecoTop4->Delete();
  h_ptRecoTop5->Delete();
  h_ptRecoTop5fine->Delete();
  h_ptRecoTop6->Delete();
  h_ptRecoTop6fine->Delete();

  h_ptGenTopMod->Delete();
  h_ptGenTopMod2->Delete();
  h_ptGenTopMod3->Delete();
  h_ptGenTopMod4->Delete();
  h_ptGenTopMod5->Delete();
  h_ptGenTopMod6->Delete();
  h_ptGenTopModDown->Delete();
  h_ptGenTopModDown2->Delete();
  h_ptGenTopModDown3->Delete();
  h_ptGenTopModDown4->Delete();
  h_ptGenTopModDown5->Delete();
  h_ptGenTopModDown6->Delete();

  h_ptRecoTopMod->Delete();
  h_ptRecoTopMod2->Delete();
  h_ptRecoTopMod3->Delete();
  h_ptRecoTopMod4->Delete();
  h_ptRecoTopMod5->Delete();
  h_ptRecoTopMod5fine->Delete();
  h_ptRecoTopMod6->Delete();
  h_ptRecoTopMod6fine->Delete();
  h_ptRecoTopModDown->Delete();
  h_ptRecoTopModDown2->Delete();
  h_ptRecoTopModDown3->Delete();
  h_ptRecoTopModDown4->Delete();
  h_ptRecoTopModDown5->Delete();
  h_ptRecoTopModDown5fine->Delete();
  h_ptRecoTopModDown6->Delete();
  h_ptRecoTopModDown6fine->Delete();

  response2.Delete();
  response3.Delete();
  response4.Delete();
  response5.Delete();
  response5fine.Delete();
  response6.Delete();
  response6fine.Delete();


  /*
  h_ngenAK4jet->Delete();
  h_genAK4jetPt->Delete();
  h_genAK4jetEta->Delete();
  h_genAK4jetPhi->Delete();
  h_genAK4jetMass->Delete();
  h_ngenAK8jet->Delete();
  h_genAK8jetPt->Delete();
  h_genAK8jetEta->Delete();
  h_genAK8jetPhi->Delete();
  h_genAK8jetMass->Delete();
  h_partLepPt->Delete();
  h_partLepEta->Delete();
  h_partLepPhi->Delete();
  */

  h_metPtPre                ->Delete();
  h_htPre                   ->Delete();
  h_htLepPre                ->Delete();
  h_nAK4jetPre              ->Delete();
  h_nBjetPre                ->Delete();
  h_ak4jetPtPre             ->Delete();
  h_ak4jetEtaPre            ->Delete();
  h_ak4jetPhiPre            ->Delete();
  h_ak4jetMassPre           ->Delete();
  h_ak4jetCSVPre            ->Delete();
  h_ak4jetVtxMassPre        ->Delete();
  h_nAK8jetPre              ->Delete();
  h_nTjetPre                ->Delete();
  h_ak8jetPtPre             ->Delete();
  h_ak8jetEtaPre            ->Delete();
  h_ak8jetPhiPre            ->Delete();
  h_ak8jetYPre              ->Delete();
  h_ak8jetMassPre           ->Delete();
  h_ak8jetMassPrunedPre     ->Delete();
  h_ak8jetMassFilteredPre   ->Delete();
  h_ak8jetMassTrimmedPre    ->Delete();
  h_ak8jetTau1Pre           ->Delete();
  h_ak8jetTau2Pre           ->Delete();
  h_ak8jetTau3Pre           ->Delete();
  h_ak8jetTau32Pre          ->Delete();
  h_ak8jetTau21Pre          ->Delete();
  h_ak8jetCSVPre            ->Delete();
  h_ak8jetSDmassPre         ->Delete();
  h_ak8jetSDsubjet0ptPre    ->Delete();
  h_ak8jetSDsubjet0massPre  ->Delete();
  h_ak8jetSDsubjet0CSVPre   ->Delete();
  h_ak8jetSDsubjet1ptPre    ->Delete();
  h_ak8jetSDsubjet1massPre  ->Delete();
  h_ak8jetSDsubjet1CSVPre   ->Delete();
  h_ak8jetSDsubjetMaxCSVPre ->Delete();
  h_ak8jetSDm01Pre          ->Delete();
  h_lepPtPre                ->Delete();
  h_lepEtaPre               ->Delete();
  h_lepAbsEtaPre            ->Delete();
  h_lepSignEtaPre           ->Delete();
  h_lepPhiPre               ->Delete();
  h_lepBJetdRPre            ->Delete();
  h_lepTJetdRPre            ->Delete();
  h_lepBJetPtRelPre         ->Delete();
  h_lepMiniIsoPre           ->Delete();
  h_leadLepPtPre            ->Delete();

  h_ak4jetEtaPre_bb->Delete();
  h_ak4jetEtaPre_b->Delete();
  h_ak4jetEtaPre_cc->Delete();
  h_ak4jetEtaPre_c->Delete();
  h_ak4jetEtaPre_l->Delete();
  h_ak8jetTau32Pre_bb->Delete();
  h_ak8jetTau32Pre_b->Delete();
  h_ak8jetTau32Pre_cc->Delete();
  h_ak8jetTau32Pre_c->Delete();
  h_ak8jetTau32Pre_l->Delete();
  h_ak8jetTau21Pre_bb->Delete();
  h_ak8jetTau21Pre_b->Delete();
  h_ak8jetTau21Pre_cc->Delete();
  h_ak8jetTau21Pre_c->Delete();
  h_ak8jetTau21Pre_l->Delete();
  h_ak8jetSDmassPre_bb->Delete();
  h_ak8jetSDmassPre_b->Delete();
  h_ak8jetSDmassPre_cc->Delete();
  h_ak8jetSDmassPre_c->Delete();
  h_ak8jetSDmassPre_l->Delete();
  h_lepEtaPre_bb->Delete();
  h_lepEtaPre_b->Delete();
  h_lepEtaPre_cc->Delete();
  h_lepEtaPre_c->Delete();
  h_lepEtaPre_l->Delete();

  h_metPt1b                ->Delete();
  h_ht1b                   ->Delete();
  h_htLep1b                ->Delete();
  h_nAK4jet1b              ->Delete();
  h_nBjet1b                ->Delete();
  h_ak4jetPt1b             ->Delete();
  h_ak4jetEta1b            ->Delete();
  h_ak4jetPhi1b            ->Delete();
  h_ak4jetMass1b           ->Delete();
  h_ak4jetCSV1b            ->Delete();
  h_ak4jetVtxMass1b        ->Delete();
  h_nAK8jet1b              ->Delete();
  h_nTjet1b                ->Delete();
  h_ak8jetPt1b             ->Delete();
  h_ak8jetEta1b            ->Delete();
  h_ak8jetPhi1b            ->Delete();
  h_ak8jetY1b              ->Delete();
  h_ak8jetMass1b           ->Delete();
  h_ak8jetMassPruned1b     ->Delete();
  h_ak8jetMassFiltered1b   ->Delete();
  h_ak8jetMassTrimmed1b    ->Delete();
  h_ak8jetTau11b           ->Delete();
  h_ak8jetTau21b           ->Delete();
  h_ak8jetTau31b           ->Delete();
  h_ak8jetTau321b          ->Delete();
  h_ak8jetTau211b          ->Delete();
  h_ak8jetCSV1b            ->Delete();
  h_ak8jetSDmass1b         ->Delete();
  h_ak8jetSDsubjet0pt1b    ->Delete();
  h_ak8jetSDsubjet0mass1b  ->Delete();
  h_ak8jetSDsubjet0CSV1b   ->Delete();
  h_ak8jetSDsubjet1pt1b    ->Delete();
  h_ak8jetSDsubjet1mass1b  ->Delete();
  h_ak8jetSDsubjet1CSV1b   ->Delete();
  h_ak8jetSDsubjetMaxCSV1b ->Delete();
  h_ak8jetSDm011b          ->Delete();
  h_lepPt1b                ->Delete();
  h_lepEta1b               ->Delete();
  h_lepAbsEta1b            ->Delete();
  h_lepSignEta1b           ->Delete();
  h_lepPhi1b               ->Delete();
  h_lepBJetdR1b            ->Delete();
  h_lepTJetdR1b            ->Delete();
  h_lepBJetPtRel1b         ->Delete();
  h_lepMiniIso1b           ->Delete();

  h_ak4jetEta1b_bb->Delete();
  h_ak4jetEta1b_b->Delete();
  h_ak4jetEta1b_cc->Delete();
  h_ak4jetEta1b_c->Delete();
  h_ak4jetEta1b_l->Delete();
  h_ak8jetTau321b_bb->Delete();
  h_ak8jetTau321b_b->Delete();
  h_ak8jetTau321b_cc->Delete();
  h_ak8jetTau321b_c->Delete();
  h_ak8jetTau321b_l->Delete();
  h_ak8jetTau211b_bb->Delete();
  h_ak8jetTau211b_b->Delete();
  h_ak8jetTau211b_cc->Delete();
  h_ak8jetTau211b_c->Delete();
  h_ak8jetTau211b_l->Delete();
  h_ak8jetSDmass1b_bb->Delete();
  h_ak8jetSDmass1b_b->Delete();
  h_ak8jetSDmass1b_cc->Delete();
  h_ak8jetSDmass1b_c->Delete();
  h_ak8jetSDmass1b_l->Delete();
  h_lepEta1b_bb->Delete();
  h_lepEta1b_b->Delete();
  h_lepEta1b_cc->Delete();
  h_lepEta1b_c->Delete();
  h_lepEta1b_l->Delete();

  h_metPt1t                ->Delete();
  h_ht1t                   ->Delete();
  h_htLep1t                ->Delete();
  h_nAK4jet1t              ->Delete();
  h_nBjet1t                ->Delete();
  h_ak4jetPt1t             ->Delete();
  h_ak4jetEta1t            ->Delete();
  h_ak4jetPhi1t            ->Delete();
  h_ak4jetMass1t           ->Delete();
  h_ak4jetCSV1t            ->Delete();
  h_ak4jetVtxMass1t        ->Delete();
  h_nAK8jet1t              ->Delete();
  h_nTjet1t                ->Delete();
  h_ak8jetPt1t             ->Delete();
  h_ak8jetEta1t            ->Delete();
  h_ak8jetPhi1t            ->Delete();
  h_ak8jetY1t              ->Delete();
  h_ak8jetMass1t           ->Delete();
  h_ak8jetMassPruned1t     ->Delete();
  h_ak8jetMassFiltered1t   ->Delete();
  h_ak8jetMassTrimmed1t    ->Delete();
  h_ak8jetTau11t           ->Delete();
  h_ak8jetTau21t           ->Delete();
  h_ak8jetTau31t           ->Delete();
  h_ak8jetTau321t          ->Delete();
  h_ak8jetTau211t          ->Delete();
  h_ak8jetCSV1t            ->Delete();
  h_ak8jetSDmass1t         ->Delete();
  h_ak8jetSDsubjet0pt1t    ->Delete();
  h_ak8jetSDsubjet0mass1t  ->Delete();
  h_ak8jetSDsubjet0CSV1t   ->Delete();
  h_ak8jetSDsubjet1pt1t    ->Delete();
  h_ak8jetSDsubjet1mass1t  ->Delete();
  h_ak8jetSDsubjet1CSV1t   ->Delete();
  h_ak8jetSDsubjetMaxCSV1t ->Delete();
  h_ak8jetSDm011t          ->Delete();
  h_lepPt1t                ->Delete();
  h_lepEta1t               ->Delete();
  h_lepAbsEta1t            ->Delete();
  h_lepSignEta1t           ->Delete();
  h_lepPhi1t               ->Delete();
  h_lepBJetdR1t            ->Delete();
  h_lepTJetdR1t            ->Delete();
  h_lepBJetPtRel1t         ->Delete();
  h_lepMiniIso1t           ->Delete();

  h_ak4jetEta1t_bb->Delete();
  h_ak4jetEta1t_b->Delete();
  h_ak4jetEta1t_cc->Delete();
  h_ak4jetEta1t_c->Delete();
  h_ak4jetEta1t_l->Delete();
  h_ak8jetTau321t_bb->Delete();
  h_ak8jetTau321t_b->Delete();
  h_ak8jetTau321t_cc->Delete();
  h_ak8jetTau321t_c->Delete();
  h_ak8jetTau321t_l->Delete();
  h_ak8jetTau211t_bb->Delete();
  h_ak8jetTau211t_b->Delete();
  h_ak8jetTau211t_cc->Delete();
  h_ak8jetTau211t_c->Delete();
  h_ak8jetTau211t_l->Delete();
  h_ak8jetSDmass1t_bb->Delete();
  h_ak8jetSDmass1t_b->Delete();
  h_ak8jetSDmass1t_cc->Delete();
  h_ak8jetSDmass1t_c->Delete();
  h_ak8jetSDmass1t_l->Delete();
  h_lepEta1t_bb->Delete();
  h_lepEta1t_b->Delete();
  h_lepEta1t_cc->Delete();
  h_lepEta1t_c->Delete();
  h_lepEta1t_l->Delete();

  h_metPt1t1b                ->Delete();
  h_ht1t1b                   ->Delete();
  h_htLep1t1b                ->Delete();
  h_nAK4jet1t1b              ->Delete();
  h_nBjet1t1b                ->Delete();
  h_ak4jetPt1t1b             ->Delete();
  h_ak4jetEta1t1b            ->Delete();
  h_ak4jetPhi1t1b            ->Delete();
  h_ak4jetMass1t1b           ->Delete();
  h_ak4jetCSV1t1b            ->Delete();
  h_ak4jetVtxMass1t1b        ->Delete();
  h_nAK8jet1t1b              ->Delete();
  h_nTjet1t1b                ->Delete();
  h_ak8jetPt1t1b             ->Delete();
  h_ak8jetEta1t1b            ->Delete();
  h_ak8jetPhi1t1b            ->Delete();
  h_ak8jetY1t1b              ->Delete();
  h_ak8jetMass1t1b           ->Delete();
  h_ak8jetMassPruned1t1b     ->Delete();
  h_ak8jetMassFiltered1t1b   ->Delete();
  h_ak8jetMassTrimmed1t1b    ->Delete();
  h_ak8jetTau11t1b           ->Delete();
  h_ak8jetTau21t1b           ->Delete();
  h_ak8jetTau31t1b           ->Delete();
  h_ak8jetTau321t1b          ->Delete();
  h_ak8jetTau211t1b          ->Delete();
  h_ak8jetCSV1t1b            ->Delete();
  h_ak8jetSDmass1t1b         ->Delete();
  h_ak8jetSDsubjet0pt1t1b    ->Delete();
  h_ak8jetSDsubjet0mass1t1b  ->Delete();
  h_ak8jetSDsubjet0CSV1t1b   ->Delete();
  h_ak8jetSDsubjet1pt1t1b    ->Delete();
  h_ak8jetSDsubjet1mass1t1b  ->Delete();
  h_ak8jetSDsubjet1CSV1t1b   ->Delete();
  h_ak8jetSDsubjetMaxCSV1t1b ->Delete();
  h_ak8jetSDm011t1b          ->Delete();
  h_lepPt1t1b                ->Delete();
  h_lepEta1t1b               ->Delete();
  h_lepAbsEta1t1b            ->Delete();
  h_lepSignEta1t1b           ->Delete();
  h_lepPhi1t1b               ->Delete();
  h_lepBJetdR1t1b            ->Delete();
  h_lepTJetdR1t1b            ->Delete();
  h_lepBJetPtRel1t1b         ->Delete();
  h_lepMiniIso1t1b           ->Delete();

  h_ak4jetEta1t1b_bb->Delete();
  h_ak4jetEta1t1b_b->Delete();
  h_ak4jetEta1t1b_cc->Delete();
  h_ak4jetEta1t1b_c->Delete();
  h_ak4jetEta1t1b_l->Delete();
  h_ak8jetTau321t1b_bb->Delete();
  h_ak8jetTau321t1b_b->Delete();
  h_ak8jetTau321t1b_cc->Delete();
  h_ak8jetTau321t1b_c->Delete();
  h_ak8jetTau321t1b_l->Delete();
  h_ak8jetTau211t1b_bb->Delete();
  h_ak8jetTau211t1b_b->Delete();
  h_ak8jetTau211t1b_cc->Delete();
  h_ak8jetTau211t1b_c->Delete();
  h_ak8jetTau211t1b_l->Delete();
  h_ak8jetSDmass1t1b_bb->Delete();
  h_ak8jetSDmass1t1b_b->Delete();
  h_ak8jetSDmass1t1b_cc->Delete();
  h_ak8jetSDmass1t1b_c->Delete();
  h_ak8jetSDmass1t1b_l->Delete();
  h_lepEta1t1b_bb->Delete();
  h_lepEta1t1b_b->Delete();
  h_lepEta1t1b_cc->Delete();
  h_lepEta1t1b_c->Delete();
  h_lepEta1t1b_l->Delete();

  h_2DisoPt15               ->Delete();
  h_2DisoPt20               ->Delete();
  h_2DisoPt25               ->Delete();
  h_2DisoPt30               ->Delete();
  h_2DisoPt35               ->Delete();
  h_2DisoPt40               ->Delete();
  h_2DisoPt45               ->Delete();

  h_2DisoScanPoints->Delete();
  h_MiniIsoScanPoints->Delete();

  h_elPtRaw->Delete();
  h_elEtaRaw->Delete();
  h_elMiniIsoRaw->Delete();
  h_elMiniIsoRaw_b->Delete();
  h_elMiniIsoRaw_e->Delete();
  h_elPtVsMiniIsoRaw->Delete();
  h_elEtaVsMiniIsoRaw->Delete();
  h_muPtRaw->Delete();
  h_muEtaRaw->Delete();
  h_muMiniIsoRaw->Delete();
  h_muMiniIsoRaw_b->Delete();
  h_muMiniIsoRaw_e->Delete();
  h_muPtVsMiniIsoRaw->Delete();
  h_muEtaVsMiniIsoRaw->Delete();

  h_lepMETdPhiRaw->Delete();
  h_ak4METdPhiRaw->Delete();
  h_ak8METdPhiRaw->Delete();

  if (!isData) {
    h_bTagSFvsPt->Delete();
    h_bTagSFvsCSV->Delete();
  }

  h_lepPhiRaw->Delete();
  h_metPhiRaw->Delete();
  h_jetPhiRaw->Delete();
  h_jetPtRaw->Delete();
  h_metPtRaw->Delete();
  
  h_muID->Delete();
  h_muTrig->Delete();
  h_elReco->Delete();
  h_elID->Delete();
  h_elIso->Delete();
  h_elTrig->Delete();
}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


