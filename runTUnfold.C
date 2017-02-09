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
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TUnfold.h"
#include "TUnfoldDensity.h"
#include "TSpline.h"
#include "TSystem.h"
#include "TPaletteAxis.h"
#include "TExec.h"
#include "RooUnfold/src/RooUnfold.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldDagostini.h"
#include "RooUnfold/src/RooUnfoldErrors.h"
#include "RooUnfold/src/RooUnfoldInvert.h"
#include "RooUnfold/src/RooUnfoldParms.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
#include "RooUnfold/src/RooUnfoldTUnfold.h"
#include "TUnfold/TUnfold.h"
#include "TUnfold/TUnfoldBinning.h"
//#include "TUnfold/TUnfoldBinningXML.h"
#include "TUnfold/TUnfoldDensity.h"
#include "TUnfold/TUnfoldSys.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;


void doUnfold(TString channel, TString toUnfold, bool isClosure, bool doNormalize, TString whichSyst, TString regMode);

void runTUnfold() {
  gSystem->Load("RooUnfold/libRooUnfold.so");
  gSystem->Load("TUnfold/libunfold.a");
  gStyle->SetOptStat(0);

  cout<<"TUnfold version is "<<TUnfold::GetTUnfoldVersion()<< endl;
  
  //Closure test
  doUnfold("mu","pt",true,false,"none","LCurve");
  doUnfold("el","pt",true,false,"none","LCurve");
  //Full unfolding

  doUnfold("mu","pt",false,false,"Up","LCurve");
  doUnfold("el","pt",false,false,"Up","LCurve");
  //doUnfold("mu","pt",false,false,"Down","LCurve");
  //doUnfold("el","pt",false,false,"Down","LCurve");

}

void doUnfold(TString channel, TString toUnfold, bool isClosure, bool doNormalize, TString whichSyst, TString regMode){

  // switch on histogram errors
  TH1::AddDirectory(kFALSE);
  TH1::SetDefaultSumw2();

  float eff_closure = 1.0;
  if (isClosure) eff_closure = 2.0;
  
  float bkg_xsecs[13] = {831.76 / 182123200.,
			 831.76 / 182075200.,
			 136.02 * 0.322 / 3279200.,
			 80.95 * 0.322 / 1682400.,
			 35.9 / 998400.,
			 35.9 / 985000.,
			 1345.0 * 1.21 / 27529599.,
			 359.7 * 1.21 / 4963240.,
			 48.91 * 1.21 / 1963464.,
			 12.05 * 1.21 / 3722395.,
			 5.501 * 1.21 / 6314257.,
			 1.329 * 1.21 / 6768156.,
			 0.03216 * 1.21 / 253561.
  };
  
  float LUMI = 12337.98;
  if (channel == "el") LUMI = 12267.67;

  const int nSYST = 9;
  TString syst_names[9] = {"lep", "pu", "JER", "JEC", "BTag", "TopTag", "PDF", "Q2", "AS"};
  
  // --------------------------------------------------------------------------------------
  // Get nominal response matrix
  // --------------------------------------------------------------------------------------

  TFile* f_nom;
  if (isClosure) f_nom = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_odd_post.root","read");
  else f_nom = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_post.root","read");

  RooUnfoldResponse * response = (RooUnfoldResponse*) f_nom->Get("response_"+toUnfold+"5fine")->Clone();
  response->SetName("response_"+toUnfold);

  TH2F* h_response = (TH2F*) response->HresponseNoOverflow();
  TVectorD tru = response->Vtruth();
  for (int j = 1; j < response->GetNbinsTruth(); j++){
    double ntru = 0.0;
    for (int i = 1; i < response->GetNbinsMeasured(); i++){
      ntru += h_response->GetBinContent(i,j);
    }
    h_response->SetBinContent(response->GetNbinsMeasured()+1, j, tru[j-1]-ntru);
  }
  
  cout << "Number of x (meas) bins, response: " << h_response->GetXaxis()->GetNbins() << endl;
  cout << "Number of y (true) bins, response: " << h_response->GetYaxis()->GetNbins() << endl;	
  f_nom->Close();
  delete f_nom;

  // --------------------------------------------------------------------------------------
  // Get 'data'
  // --------------------------------------------------------------------------------------

  TFile* f_data;
  if (!isClosure) f_data = TFile::Open("histfiles_80X/hists_Data_"+channel+".root","read");
  else f_data = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_even_post.root","read");
  TH1F* h_data = (TH1F*) f_data->Get(toUnfold+"RecoTop5fine")->Clone();
  f_data->Close();
  delete f_data;
  cout << "Number of data bins: " << h_data->GetXaxis()->GetNbins() << endl;

  // Correct for fake events
  if (response->FakeEntries()){
    TVectorD fakes = response->Vfakes();
    double fac = response->Vmeasured().Sum();
    if (fac != 0.0) {
      TVectorD* measVec = RooUnfoldResponse::H2V(h_data,response->GetNbinsMeasured(),0);
      fac = measVec->Sum() / fac;
    }
    for (int ii = 1; ii < response->GetNbinsMeasured(); ii++){
      h_data->SetBinContent(ii,h_data->GetBinContent(ii) - fac * fakes[ii-1]);
    }
  }

  if (isClosure) h_data->Scale(bkg_xsecs[0] * LUMI * eff_closure);

  TH1F* hMeas = (TH1F*) h_data->Clone();

  // --------------------------------------------------------------------------------------
  // Get 'truth'
  // --------------------------------------------------------------------------------------  

  TFile* f_true;
  if (!isClosure) f_true = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_post.root","read");
  else f_true = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_even_post.root","read");
  TH1F* hTrue = (TH1F*) f_true->Get(toUnfold+"GenTop5")->Clone();
  f_true->Close();
  delete f_true;
  
  hTrue->Scale(bkg_xsecs[0] * LUMI * eff_closure);
  hTrue->Scale(1.0/(LUMI*0.438/3.));

  for (int ib = 1; ib < hTrue->GetNbinsX()+1; ib++){
    hTrue->SetBinContent(ib,hTrue->GetBinContent(ib) / hTrue->GetBinWidth(ib));
    hTrue->SetBinError(ib,hTrue->GetBinError(ib) / hTrue->GetBinWidth(ib));
  }

  // --------------------------------------------------------------------------------------
  // Get backgrounds if not closure test
  // --------------------------------------------------------------------------------------

  TH1F* h_ttbar;
  TH1F* h_singletop;
  TH1F* h_wjets;
  TH1F* h_qcd;
  
  if (!isClosure) {
    
    TFile* f_ttbar = TFile::Open("histfiles_80X/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_post.root","read");
    h_ttbar = (TH1F*) f_ttbar->Get(toUnfold+"RecoTop5fine")->Clone();
    h_ttbar->Sumw2();
    h_ttbar->Scale(bkg_xsecs[1] * LUMI);
    f_ttbar->Close();
    delete f_ttbar;
    
    TString singletop_filenames[4] = {
      "histfiles_80X/hists_SingleTop_t_t_"+channel+"_nom_post.root",
      "histfiles_80X/hists_SingleTop_tbar_t_"+channel+"_nom_post.root",
      "histfiles_80X/hists_SingleTop_t_tW_"+channel+"_nom_post.root",
      "histfiles_80X/hists_SingleTop_tbar_tW_"+channel+"_nom_post.root"
    };
    
    for (int jj = 0; jj < 4; jj++){
      TFile* f_singletop = TFile::Open(singletop_filenames[jj],"read");
      TH1F* h_tmp = (TH1F*) f_singletop->Get(toUnfold+"RecoTop5fine")->Clone();
      h_tmp->Sumw2();
      h_tmp->Scale(bkg_xsecs[jj+2] * LUMI);
      if (jj == 0) h_singletop = (TH1F*) h_tmp->Clone();
      else h_singletop->Add(h_tmp);
      f_singletop->Close();
      delete f_singletop;
    }
    
    TString wjets_filenames[7] = {
      "histfiles_80X/hists_WJets_HT100to200_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT200to400_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT400to600_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT600to800_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT800to1200_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT1200to2500_"+channel+"_nom_post.root",
      "histfiles_80X/hists_WJets_HT2500toInf_"+channel+"_nom_post.root"
    };
    
    for (int jj = 0; jj < 7; jj++){
      TFile* f_wjets = TFile::Open(wjets_filenames[jj],"read");
      TH1F* h_tmp = (TH1F*) f_wjets->Get(toUnfold+"RecoTop5fine")->Clone();
      h_tmp->Sumw2();
      h_tmp->Scale(bkg_xsecs[jj+6] * LUMI);
      if (jj == 0) h_wjets = (TH1F*) h_tmp->Clone();
      else h_wjets->Add(h_tmp);
      f_wjets->Close();
      delete f_wjets;
    }
    
    TString qcd_filenames[14] = {
      "histfiles_80X/hists_Data_"+channel+"_qcd.root",
      "histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_PowhegPythia8_nonsemilep_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_SingleTop_t_t_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_SingleTop_tbar_t_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_SingleTop_t_tW_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_SingleTop_tbar_tW_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT100to200_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT200to400_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT400to600_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT600to800_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT800to1200_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT1200to2500_"+channel+"_nom_qcd_post.root",
      "histfiles_80X/hists_WJets_HT2500toInf_"+channel+"_nom_qcd_post.root"
    };
    
    for (int jj = 0; jj < 14; jj++){
      TFile* f_qcd = TFile::Open(qcd_filenames[jj],"read");
      TH1F* h_tmp = (TH1F*) f_qcd->Get(toUnfold+"RecoTop5fine")->Clone();
      h_tmp->Sumw2();
      if (jj == 0) h_qcd = (TH1F*) h_tmp->Clone();
      else {
	h_tmp->Scale(bkg_xsecs[jj-1] * LUMI);
	h_qcd->Add(h_tmp,-1.0);
      }
      f_qcd->Close();
      delete f_qcd;
    }

    // -------------------------------
    // Normalize QCD to MC prediction
    // -------------------------------
    
    float QCD_norm = 0.0;
    TString qcd_MC_filenames[6] = {"histfiles_80X/hists_QCD_HT300to500_"+channel+"_nom_post.root",
				   "histfiles_80X/hists_QCD_HT500to700_"+channel+"_nom_post.root",
				   "histfiles_80X/hists_QCD_HT700to1000_"+channel+"_nom_post.root",
				   "histfiles_80X/hists_QCD_HT1000to1500_"+channel+"_nom_post.root",
				   "histfiles_80X/hists_QCD_HT1500to2000_"+channel+"_nom_post.root",
				   "histfiles_80X/hists_QCD_HT2000toInf_"+channel+"_nom_post.root"};
    float QCDnorms[6] = {347700. / 37828442., 32100. / 44058594., 6831. / 29832311., 1207. / 4963881., 119.9 / 7803965., 25.24 / 4047532.};
    for (int iq = 0; iq < 6; iq++){
      TFile* f_tmp = TFile::Open(qcd_MC_filenames[iq],"read");
      TH1F* h_tmp = (TH1F*) f_tmp->Get(toUnfold+"RecoTop5fine");
      h_tmp->Scale(QCDnorms[iq] * LUMI);
      QCD_norm += h_tmp->Integral();
    }
    h_qcd->Scale(QCD_norm / h_qcd->Integral());
  }

  if (!isClosure){
    if (channel == "mu"){
      hMeas->Add(h_ttbar,-1.0*0.77);
      hMeas->Add(h_singletop,-1.0*1.00);
      hMeas->Add(h_wjets,-1.0*1.02);
      hMeas->Add(h_qcd,-1.0*0.62);
    }
    else {
      hMeas->Add(h_ttbar,-1.0*0.77*1.79);
      hMeas->Add(h_singletop,-1.0*1.00*1.79);
      hMeas->Add(h_wjets,-1.0*1.02*1.79);
      hMeas->Add(h_qcd,-1.0*0.67*1.79);
    }
  }
  
  hMeas->Scale(1.0/(LUMI*0.438/3.));

  for (int ib = 1; ib < hMeas->GetNbinsX()+1; ib++){
    hMeas->SetBinContent(ib,hMeas->GetBinContent(ib) / hMeas->GetBinWidth(ib));
    hMeas->SetBinError(ib,hMeas->GetBinError(ib) / hMeas->GetBinWidth(ib));
  }      

  // --------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------
  // Do full unfolding!!
  // --------------------------------------------------------------------------------------
  // --------------------------------------------------------------------------------------

  TUnfoldDensity::ERegMode regModes[2] = {TUnfoldDensity::kRegModeDerivative, TUnfoldDensity::kRegModeCurvature};
  const int nSCAN = (regMode == "LCurve") ? 1 : 3;
  TUnfoldDensity::EScanTauMode scan_modes[3] = {TUnfoldDensity::kEScanTauRhoAvg,TUnfoldDensity::kEScanTauRhoAvgSys,TUnfoldDensity::kEScanTauRhoSquareAvgSys};

  for (int ii = 0; ii < 2; ii++){
    TUnfoldDensity unfold(h_response,TUnfoldDensity::kHistMapOutputVert,regModes[ii],TUnfoldDensity::kEConstraintNone,TUnfoldDensity::kDensityModeBinWidth);
    cout << "Made unfold object!" << endl;
    unfold.SetInput(h_data);
    cout << "Added data input!" << endl;

    // --------------------------------------------------------------------------------------
    // Subtract backgrounds if not closure test
    // --------------------------------------------------------------------------------------
    
    if (!isClosure){
      if (channel == "mu"){
	unfold.SubtractBackground(h_ttbar,"ttbar",0.77,0.06); //hist, bkg name, fit normalization, uncertainty
	unfold.SubtractBackground(h_singletop,"singletop",1.00,0.05);
	unfold.SubtractBackground(h_wjets,"wjets",1.02,0.06);
	unfold.SubtractBackground(h_qcd,"qcd",0.62,0.17);
      }
      else {
	unfold.SubtractBackground(h_ttbar,"ttbar",0.77*1.79,0.06*1.79); //same as above, w/ extra fudge factor from fit
	unfold.SubtractBackground(h_singletop,"singletop",1.00*1.79,0.05*1.79);
	unfold.SubtractBackground(h_wjets,"wjets",1.02*1.79,0.06*1.79);
	unfold.SubtractBackground(h_qcd,"qcd",0.67*1.79,0.35*1.79);
      }
    }

    // --------------------------------------------------------------------------------------
    // Get systematics
    // --------------------------------------------------------------------------------------
    
    if (whichSyst != "none"){
      for (int jj = 0; jj < nSYST; jj++){
	
	TFile* f_syst = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_"+syst_names[jj]+whichSyst+"_post.root","read");
	
	RooUnfoldResponse* response_syst = (RooUnfoldResponse*) f_syst->Get("response_"+toUnfold+"5fine")->Clone();
	response_syst->SetName("response_"+toUnfold+"_"+syst_names[jj]);
	
	TH2F* h_response_syst = (TH2F*) response_syst->HresponseNoOverflow();
	TVectorD tru_syst = response_syst->Vtruth();
	for (int j = 1; j < response_syst->GetNbinsTruth(); j++){
	  double ntru_syst = 0.0;
	  for (int i = 1; i < response_syst->GetNbinsMeasured(); i++){
	    ntru_syst += h_response_syst->GetBinContent(i,j);
	  }
	  h_response_syst->SetBinContent(response_syst->GetNbinsMeasured()+1, j, tru_syst[j-1]-ntru_syst);
	}
	
	unfold.AddSysError(h_response_syst, syst_names[jj], TUnfoldDensity::kHistMapOutputVert, TUnfoldDensity::kSysErrModeMatrix);
	f_syst->Close();
	delete f_syst;
      }
    }
    
    // --------------------------------------------------------------------------------------
    // Find tau, do unfolding
    // --------------------------------------------------------------------------------------

    for (int jj = 0; jj < nSCAN; jj++){

      TCanvas* c1 = new TCanvas();
      c1->cd();

      if (regMode == "LCurve"){
	TSpline* logTauX = 0;
	TSpline* logTauY = 0;
	TGraph* lCurve = 0;
	TGraph* bestLCurve = new TGraph(1);
	int i_tau = unfold.ScanLcurve(100,0.0001,0.01,&lCurve,&logTauX,&logTauY);
	Double_t tau, x ,y, rho;
	logTauX->GetKnot(i_tau,tau,x);
	logTauY->GetKnot(i_tau,tau,y);
	bestLCurve->SetPoint(1,x,y);
	bestLCurve->SetMarkerColor(2);
	lCurve->Draw();
	bestLCurve->Draw("*");
      }
      
      else {
	TSpline* scanResult = 0;
	int i_tau = unfold.ScanTau(100,0.0001,0.01,&scanResult,scan_modes[jj]);
	
	Double_t t[1],rho[1],x[1],y[1];
	scanResult->GetKnot(i_tau,t[0],rho[0]);
	TGraph *bestRhoLogTau=new TGraph(1,t,rho);
	Double_t *tAll=new Double_t[100];
	Double_t *rhoAll=new Double_t[100];
	for(Int_t i=0;i<100;i++) {
	  scanResult->GetKnot(i,tAll[i],rhoAll[i]);
	}
	TGraph *knots=new TGraph(100,tAll,rhoAll);
	knots->Draw();
	bestRhoLogTau->Draw("*");
      }
      
      c1->Print(TString::Format("UnfoldingPlots/"+regMode+"_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      
      cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
	  <<" / "<<unfold.GetNdf()<<"\n";
      
      // --------------------------------------------------------------------------------------
      // Get central result + combined uncertainties
      // --------------------------------------------------------------------------------------
      
      TH1F* hReco_stat = (TH1F*) unfold.GetOutput("Unfolded");
      TH1F* hReco_total = (TH1F*) hReco_stat->Clone();
      
      TH2* histEmatStat = unfold.GetEmatrixInput("unfolding stat error matrix"); // For stat unc
      TH2* histEmatTotal = unfold.GetEmatrixTotal("unfolding total error matrix"); // For total unc
      
      for(int ii = 0; ii < histEmatStat->GetNbinsX()+2; ii++) {
	hReco_stat->SetBinError(ii,TMath::Sqrt(histEmatStat->GetBinContent(ii,ii)));
	hReco_total->SetBinError(ii,TMath::Sqrt(histEmatTotal->GetBinContent(ii,ii)));
      }
      
      // --------------------------------------------------------------------------------------
      // Get uncertainty breakdown
      // --------------------------------------------------------------------------------------
      
      TH1F* h_sysUnc[nSYST];
      TH1F* h_bkgUnc;
      
      if (whichSyst != "none"){
	for (int kk = 0; kk < nSYST; kk++){
	  h_sysUnc[kk] = (TH1F*) unfold.GetDeltaSysSource(syst_names[kk],"h_"+syst_names[kk]);
	  for (int ib = 1; ib < h_sysUnc[kk]->GetNbinsX()+1; ib++){
	    h_sysUnc[kk]->SetBinContent(ib,fabs(h_sysUnc[kk]->GetBinContent(ib) / hReco_stat->GetBinContent(ib) * 100.));
	    h_sysUnc[kk]->SetBinError(ib,0.0);
	  }
	}
	
	if (!isClosure){
	  TH1F* h_bkgUnc_tmp[4];
	  TString bkgs[4] = {"ttbar","singletop","wjets","qcd"};
	  for (int kk = 0; kk < 4; kk++){
	    h_bkgUnc_tmp[kk] = (TH1F*) unfold.GetDeltaSysBackgroundScale(bkgs[kk],"h_bkg_"+bkgs[kk]);
	    if (kk == 0) {
	      h_bkgUnc = (TH1F*) h_bkgUnc_tmp[kk]->Clone();
	      h_bkgUnc->Reset();
	    }
	    for (int ib = 1; ib < h_bkgUnc_tmp[kk]->GetNbinsX()+1; ib++){
	      h_bkgUnc_tmp[kk]->SetBinContent(ib,fabs((h_bkgUnc_tmp[kk]->GetBinContent(ib) - hReco_stat->GetBinContent(ib))/hReco_stat->GetBinContent(ib) * 100.));
	    }
	  }
	  
	  for (int ib = 1; ib < h_bkgUnc->GetNbinsX()+1; ib++){
	    h_bkgUnc->SetBinContent(ib,sqrt(pow(2,h_bkgUnc_tmp[0]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[1]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[2]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[3]->GetBinContent(ib))));
	    
	  }
	}
      }
      
      // --------------------------------------------------------------------------------------
      // Start normalizing
      // --------------------------------------------------------------------------------------
      
      hReco_stat->Sumw2();
      hReco_stat->Scale(1.0/(LUMI*0.438/3.));
      hReco_total->Sumw2();
      hReco_total->Scale(1.0/(LUMI*0.438/3.));
      
      for (int ib = 1; ib < hReco_stat->GetNbinsX()+1; ib++){
	hReco_stat->SetBinContent(ib,hReco_stat->GetBinContent(ib) / hReco_stat->GetBinWidth(ib));
	hReco_stat->SetBinError(ib,hReco_stat->GetBinError(ib) / hReco_stat->GetBinWidth(ib));
	hReco_total->SetBinContent(ib,hReco_total->GetBinContent(ib) / hReco_total->GetBinWidth(ib));
	hReco_total->SetBinError(ib,hReco_total->GetBinError(ib) / hReco_total->GetBinWidth(ib));
      }
      
      TH1F* hFrac = (TH1F*) hReco_stat->Clone();
      hFrac->Sumw2();
      hFrac->Divide(hTrue);
      
      TH1F* hStatUnc = (TH1F*) hFrac->Clone();
      TH1F* hTotalUnc = (TH1F*) hFrac->Clone();
      
      for (int ib = 1; ib < hFrac->GetNbinsX()+1; ib++){
	hStatUnc->SetBinContent(ib,1.0);
	hStatUnc->SetBinError(ib,hReco_stat->GetBinError(ib) / hReco_stat->GetBinContent(ib));
	hTotalUnc->SetBinContent(ib,1.0);
	hTotalUnc->SetBinError(ib,hReco_total->GetBinError(ib) / hReco_total->GetBinContent(ib));
      }
            
      // --------------------------------------------------------------------------------------
      // Plot unfolding
      // --------------------------------------------------------------------------------------
      
      Float_t f_textSize=0.057;
      gStyle->SetPadLeftMargin(0.12);
      gStyle->SetPadRightMargin(0.1);
      gStyle->SetTextSize(f_textSize);
      gStyle->SetLabelSize(f_textSize,"x");
      gStyle->SetTitleSize(f_textSize,"x");
      gStyle->SetLabelSize(f_textSize,"y");
      gStyle->SetTitleSize(f_textSize,"y");
      gStyle->SetOptStat(0);
      gStyle->SetEndErrorSize(5);
      
      TCanvas *c2 = new TCanvas("c2", "", 700, 625);
      
      float ratio_size = 0.33;
      
      TPad* p1 = new TPad("p1","p1",0,ratio_size,1,1);
      p1->SetBottomMargin(0.03);
      p1->SetTopMargin(0.10);
      p1->SetLeftMargin(0.16);
      p1->SetRightMargin(0.075);
      
      TPad* p2 = new TPad("p2","p2",0,0,1,ratio_size);
      p2->SetTopMargin(0.01);
      p2->SetBottomMargin(0.36);
      p2->SetLeftMargin(0.16);
      p2->SetRightMargin(0.075);
      
      p2->Draw();
      p1->Draw();
      
      p1->cd();
      
      TString xsec_title = ";;d#sigma/dp_{T} [fb/GeV]";
      if (doNormalize) xsec_title = ";;1/#sigma d#sigma/dp_{T} [1/GeV]";
      
      hTrue->SetTitle(xsec_title);
      hTrue->GetYaxis()->SetTitleSize(0.075);    
      hTrue->GetYaxis()->SetTitleOffset(1.1);
      hTrue->GetYaxis()->SetLabelSize(0.065);
      hTrue->GetXaxis()->SetTitleSize(0);
      hTrue->GetXaxis()->SetLabelSize(0);
      
      hReco_total->SetMarkerStyle(21);
      hReco_stat->SetMarkerStyle(21);
      hMeas->SetMarkerStyle(25);
      
      hReco_total->GetXaxis()->SetRangeUser(400.,1200.);
      hReco_stat->GetXaxis()->SetRangeUser(400.,1200.);
      hTrue->GetXaxis()->SetRangeUser(400.,1200.);
      hMeas->GetXaxis()->SetRangeUser(400.,1200.);
      
      hTrue->Draw("hist");
      TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
      setex1->Draw();
      hReco_total->Draw("hist,ep,same");
      hReco_stat->Draw("hist,E1p,same");
      TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
      setex2->Draw();
      hMeas->Draw("same");
      
      TLegend* leg = new TLegend(0.5, 0.55, 0.9, 0.75);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->SetBorderSize(0);
      
      if (!isClosure){ 
	leg->AddEntry( hReco_total, "Unfolded data", "p");
	leg->AddEntry( hTrue, "Generated (Powheg)", "l");
	leg->AddEntry( hMeas, "Raw data", "p");
	leg->AddEntry( hStatUnc,"Stat. uncertainty","f");
	leg->AddEntry( hTotalUnc,"Stat. #oplus syst. uncertainties","f");
      }
      else {
	leg->AddEntry( hReco_total, "Unfolded MC (Powheg)", "p");
	leg->AddEntry( hTrue, "Generated (Powheg)", "l");
	leg->AddEntry( hMeas, "Reco-level (Powheg)", "p");
	leg->AddEntry( hStatUnc,"Stat. uncertainty","f");
	leg->AddEntry( hTotalUnc,"Stat. #oplus syst. uncertainties","f");
      }
      
      leg->Draw();
      
      TLatex* text1 = new TLatex();
      text1->SetNDC();
      text1->SetTextFont(42);
      text1->DrawLatex(0.55,0.8, "#scale[1.0]{L = 12.4 fb^{-1}, #sqrt{s} = 13 TeV}");
      
      // ----------------------------------------------------------------------------------------------------------------
      
      p2->cd();
      
      hFrac->GetXaxis()->SetTitleSize(0.15);
      hFrac->GetYaxis()->SetTitleSize(0.13);
      hFrac->GetXaxis()->SetLabelSize(0.12);
      hFrac->GetYaxis()->SetLabelSize(0.11);
      
      hFrac->GetYaxis()->SetTitle("Theory / Data");
      hFrac->GetXaxis()->SetTitle("Top quark p_{T} (GeV)");
      hFrac->GetYaxis()->SetTitleOffset(0.8);
      hFrac->GetYaxis()->SetNdivisions(505);
      hFrac->GetYaxis()->SetTitleOffset(0.38);
      hFrac->GetXaxis()->SetTitleOffset(1.0);
      
      hFrac->SetAxisRange(0.0,2.0,"Y");
      
      hTotalUnc->SetFillColor(kAzure-9);
      hTotalUnc->SetFillStyle(1001);
      hStatUnc->SetFillColor(18);
      hStatUnc->SetFillStyle(1001);
      
      hFrac->Draw("hist");
      hTotalUnc->Draw("same,e2");
      hStatUnc->Draw("same,e2");
      hFrac->Draw("same,hist");
      hFrac->GetXaxis()->SetRangeUser(400., 1200.);
      
      c2->Update();
      c2->Print(TString::Format("UnfoldingPlots/unfolded_ttbar_xs_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      
      // --------------------------------------------------------------------------------------
      // Plot uncertainty breakdown
      // --------------------------------------------------------------------------------------
      
      if (whichSyst != "none"){
	TCanvas *c3 = new TCanvas("c3", "", 800, 600);
	c3->SetTopMargin(0.08);
	c3->SetRightMargin(0.05);
	c3->SetBottomMargin(0.14);
	c3->SetLeftMargin(0.16);
	
	TH1F* h_sysTotal = (TH1F*) hTotalUnc->Clone();
	TH1F* h_sysStat = (TH1F*) hStatUnc->Clone();
	for (int ib = 1; ib < hTotalUnc->GetNbinsX()+1; ib++){
	  h_sysTotal->SetBinContent(ib, hTotalUnc->GetBinError(ib) * 100.);
	  h_sysStat->SetBinContent(ib, hStatUnc->GetBinError(ib) * 100.);
	}
	
	if (toUnfold == "pt") h_sysTotal->GetXaxis()->SetTitle("Top quark p_{T} (GeV)");
	else if (toUnfold == "y") h_sysTotal->GetXaxis()->SetTitle("Top quark y");
	h_sysTotal->GetYaxis()->SetTitle("Uncertainty [%]");
	if (toUnfold == "pt") h_sysTotal->SetAxisRange(400,1150,"X");
	h_sysTotal->SetMaximum(100.0);
	h_sysTotal->SetMinimum(0.0);
	
	int colors[9] = {632,600,617,417,432,801,864,906,419};
	int markers[9] = {20,21,22,23,33,24,25,26,32};
	for (int is = 0; is < nSYST; is++){
	  h_sysUnc[is]->SetLineColor(colors[is]);
	  h_sysUnc[is]->SetLineWidth(2);
	  h_sysUnc[is]->SetMarkerColor(colors[is]);
	  h_sysUnc[is]->SetMarkerStyle(markers[is]);
	}
	
	if (!isClosure){
	  h_bkgUnc->SetLineColor(426);
	  h_bkgUnc->SetLineWidth(2);
	  h_bkgUnc->SetMarkerColor(426);
	  h_bkgUnc->SetMarkerStyle(27);
	}

	h_sysStat->SetLineColor(1);
	h_sysStat->SetLineWidth(2);
	h_sysStat->SetMarkerColor(1);
	h_sysStat->SetMarkerStyle(34);
	
	TString longnames[9] = {"Lepton ID","Pileup","Jet energy scale","Jet energy resolution",
				"b tagging efficiency","t tagging efficiency",
				"PDF Uncertainty","#mu_{R}, #mu_{F} scales","#alpha_{s} scale"};
	
	TLegend* leg2;
	if (toUnfold == "pt" && doNormalize) leg2 = new TLegend(0.2,0.65,0.45,0.8);
	else leg2 = new TLegend(0.2,0.65,0.45,0.88);
	const int split = (int) nSYST / 2;
	for (int is = 0; is < split+1; is++){
	  leg2->AddEntry(h_sysUnc[is],longnames[is],"lp");
	}
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.04);
	leg2->SetTextFont(42);
	
	TLegend* leg22;
	if (toUnfold == "pt" && doNormalize) leg22 = new TLegend(0.55,0.65,0.8,0.8);
	else leg22 = new TLegend(0.55,0.65,0.8,0.88);
	for (int is = split+1; is < nSYST; is++){
	  leg22->AddEntry(h_sysUnc[is],longnames[is],"lp");
	}
	if (!isClosure) leg22->AddEntry(h_bkgUnc,"Background normalization","lp");
	leg22->AddEntry(h_sysStat,"Statistical uncertainty","lp");
	leg22->AddEntry(h_sysTotal,"Total uncertainty","f");
	
	leg22->SetFillStyle(0);
	leg22->SetBorderSize(0);
	leg22->SetTextSize(0.04);
	leg22->SetTextFont(42);
	
	h_sysTotal->SetFillStyle(3344);

	for (int is = 0; is < nSYST; is++){
	  h_sysUnc[is]->Draw("hist");
	}
	h_sysTotal->Draw("hist");
	for (int is = 0; is < nSYST; is++){
	  h_sysUnc[is]->Draw("ep,same");
	}
	if (!isClosure) h_bkgUnc->Draw("ep,same");
	h_sysStat->Draw("ep,same");
	h_sysTotal->Draw("axis,same");
	leg2->Draw(); 
	leg22->Draw(); 
	
	c3->SaveAs(TString::Format("UnfoldingPlots/unfold_relative_uncertainties_"+toUnfold+"_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      }
      
      // --------------------------------------------------------------------------------------
      // Plot covariance matrix
      // --------------------------------------------------------------------------------------
      
      // Normalize diagonal
      float diags[10] = {0.};
      for (int id = 0; id < histEmatTotal->GetXaxis()->GetNbins(); id++){
	diags[id] = sqrt(histEmatTotal->GetBinContent(id+1,id+1));
      }
      for (int ix = 1; ix < histEmatTotal->GetXaxis()->GetNbins()+1; ix++){
	for (int iy = 1; iy < histEmatTotal->GetYaxis()->GetNbins()+1; iy++){
	  if (diags[ix-1] > 0.0 && diags[iy-1] >= 0.0){
	    float rawval = histEmatTotal->GetBinContent(ix,iy);
	    rawval *= 100.0 / diags[ix-1] / diags[iy-1];
	    int roundval = rawval * 1000;
	    float finalval = (float) roundval / 1000.;
	    histEmatTotal->SetBinContent(ix,iy,finalval);
	  }
	  else histEmatTotal->SetBinContent(ix,iy,0.0);
	}
      }
      
      /*
	TString binnames[5] = {"Bin 1","Bin 2","Bin 3","Bin 4","Bin 5"};
	for (int ii = 0; ii < histEmatTotal->GetXaxis()->GetNbins() - 2; ii++){
	histEmatTotal->GetXaxis()->SetBinLabel(ii+2,binnames[ii]);
	histEmatTotal->GetYaxis()->SetBinLabel(ii+2,binnames[ii]);
	}
      */
      
      TCanvas* c4 = new TCanvas("c4","c4",900,700);
      c4->SetLeftMargin(0.18);
      c4->SetRightMargin(0.12);
      c4->SetBottomMargin(0.1);
      //histEmatTotal->SetAxisRange(1.0,histEmatTotal->GetXaxis()->GetNbins()-1.001,"X");
      //histEmatTotal->SetAxisRange(1.0,histEmatTotal->GetXaxis()->GetNbins()-1.001,"Y");
      histEmatTotal->GetZaxis()->SetLabelSize(0.04);
      histEmatTotal->Draw("COLZ");
      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)histEmatTotal->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.88);
      palette->SetX2NDC(0.93);
      palette->SetY1NDC(0.1);
      palette->SetY2NDC(0.9);
      gPad->Modified();
      gPad->Update();
      
      c4->SaveAs(TString::Format("UnfoldingPlots/total_covariance_"+toUnfold+"_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      
    }// End mode loop
  }// End regularization loop
}
