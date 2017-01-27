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


void doUnfold(TString channel, TString toUnfold, bool isClosure, bool doNormalize, TString whichSyst);

void runTUnfold() {
  gSystem->Load("RooUnfold/libRooUnfold.so");
  gSystem->Load("TUnfold/libunfold.a");
  gStyle->SetOptStat(0);

  cout<<"TUnfold version is "<<TUnfold::GetTUnfoldVersion()<< endl;
  
  //Closure test
  doUnfold("mu","pt",true,false,"none");
  doUnfold("el","pt",true,false,"none");
  //Full unfolding
  /*
  doUnfold("mu","pt",false,false,"none");
  doUnfold("el","pt",false,false,"none");
  doUnfold("mu","pt",false,false,"Up");
  doUnfold("el","pt",false,false,"Up");
  doUnfold("mu","pt",false,false,"Down");
  doUnfold("el","pt",false,false,"Down");
  */
}

void doUnfold(TString channel, TString toUnfold, bool isClosure, bool doNormalize, TString whichSyst){

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
  if (channel == "el") LUMI = 12267.67; //Fudge factor from fit

  TString syst_names[9] = {"lep", "pu", "JER", "JEC", "BTag", "TopTag", "PDF", "Q2", "AS"};
  
  // --------------------------------------------------------------------------------------
  // Get nominal response matrix
  // --------------------------------------------------------------------------------------

  TFile* f_nom;
  if (isClosure) f_nom = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_even_post.root","read");
  else f_nom = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_post.root","read");

  RooUnfoldResponse * response = (RooUnfoldResponse*) f_nom->Get("response_"+toUnfold+"5fine")->Clone();
  response->SetName("response_"+toUnfold);
  TH2F* h_response = (TH2F*) response->Hresponse();
  cout << "Number of x (meas) bins, response: " << h_response->GetXaxis()->GetNbins() << endl;
  cout << "Number of y (true) bins, response: " << h_response->GetYaxis()->GetNbins() << endl;	
  f_nom->Close();
  delete f_nom;

  TCanvas output;
  output.Divide(2);
  output.cd(1);
  h_response->SetLineColor(kBlue);
  h_response->Draw("BOX");
  
  // --------------------------------------------------------------------------------------
  // Get 'data'
  // --------------------------------------------------------------------------------------

  TFile* f_data;
  if (!isClosure) f_data = TFile::Open("histfiles_80X/hists_Data_"+channel+".root","read");
  else f_data = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_odd_post.root","read");
  TH1F* h_data = (TH1F*) f_data->Get(toUnfold+"RecoTop5fine")->Clone();
  f_data->Close();
  delete f_data;
  cout << "Number of data bins: " << h_data->GetXaxis()->GetNbins() << endl;

  if (isClosure) h_data->Scale(bkg_xsecs[0] * LUMI * eff_closure);

  output.cd(2);
  h_data->Draw();
  output.SaveAs("test.pdf");

  // --------------------------------------------------------------------------------------
  // Get 'truth'
  // --------------------------------------------------------------------------------------  

  TFile* f_true;
  if (!isClosure) f_true = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_post.root","read");
  else f_true = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_nom_odd_post.root","read");
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

  // --------------------------------------------------------------------------------------
  // Set up initial TUnfold object
  // --------------------------------------------------------------------------------------

  TUnfoldDensity::ERegMode regModes[2] = {TUnfoldDensity::kRegModeDerivative, TUnfoldDensity::kRegModeCurvature};
  const int nSCAN = 2;
  TUnfoldDensity::EScanTauMode scan_modes[nSCAN] = {TUnfoldDensity::kEScanTauRhoAvg, TUnfoldDensity::kEScanTauRhoMax};

  for (int ii = 0; ii < 2; ii++){
    TUnfoldDensity unfold(h_response,TUnfoldDensity::kHistMapOutputVert,regModes[ii],TUnfoldDensity::kEConstraintArea,TUnfoldDensity::kDensityModeBinWidth);
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
    
      // --------------------------------------------------------------------------------------
      // Get systematics
      // --------------------------------------------------------------------------------------
            
      if (whichSyst != "none"){
	for (int jj = 0; jj < 9; jj++){
	  
	  TFile* f_syst = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_"+channel+"_"+syst_names[jj]+whichSyst+"_post.root","read");
	  
	  RooUnfoldResponse* response_syst = (RooUnfoldResponse*) f_syst->Get("response_"+toUnfold+"5fine")->Clone();
	  response_syst->SetName("response_"+toUnfold+"_"+syst_names[jj]);
	  TH2F* h_response_syst = (TH2F*) response_syst->Hresponse();
	  
	  unfold.AddSysError(h_response_syst, syst_names[jj], TUnfoldDensity::kHistMapOutputVert, TUnfoldDensity::kSysErrModeMatrix);
	  f_syst->Close();
	  delete f_syst;
	}
      }
    }
      
    // --------------------------------------------------------------------------------------
    // Find tau, do unfolding
    // --------------------------------------------------------------------------------------

    for (int jj = 0; jj < nSCAN; jj++){
      TSpline* scanResult = 0;
      int i_tau = unfold.ScanTau(100,0.0,0.0,&scanResult,scan_modes[jj]);

      Double_t t[1],rho[1],x[1],y[1];
      scanResult->GetKnot(i_tau,t[0],rho[0]);
      TGraph *bestRhoLogTau=new TGraph(1,t,rho);
      Double_t *tAll=new Double_t[100];
      Double_t *rhoAll=new Double_t[100];
      for(Int_t i=0;i<100;i++) {
	scanResult->GetKnot(i,tAll[i],rhoAll[i]);
      }
      TGraph *knots=new TGraph(100,tAll,rhoAll);

      TCanvas* c1 = new TCanvas();
      c1->cd();
      scanResult->Draw();
      knots->Draw("*");
      bestRhoLogTau->SetMarkerColor(kRed);
      bestRhoLogTau->Draw("*");
      c1->Print(TString::Format("UnfoldingPlots/scan_tau_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      
      cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
	  <<" / "<<unfold.GetNdf()<<"\n";

      // --------------------------------------------------------------------------------------
      // Get central result + combined uncertainties
      // --------------------------------------------------------------------------------------

      //unfold.DoUnfold(0.01); HACK -- tried to override best tau
      TH1F* hReco = (TH1F*) unfold.GetOutput("Unfolded");
      
      hReco->Scale(1.0/(LUMI*0.438/3.));

      for (int ib = 1; ib < hReco->GetNbinsX()+1; ib++){
	hReco->SetBinContent(ib,hReco->GetBinContent(ib) / hReco->GetBinWidth(ib));
	hReco->SetBinError(ib,hReco->GetBinError(ib) / hReco->GetBinWidth(ib));
      }

      TH1F* hFrac = (TH1F*) hReco->Clone();
      hFrac->Sumw2();
      hFrac->Divide(hTrue);
      
      // --------------------------------------------------------------------------------------
      // Plot unfolding
      // --------------------------------------------------------------------------------------
      
      TCanvas *c2 = new TCanvas("c", "c", 700, 700);
      TPad *pad1 =  new TPad("pad1","pad1",0,0.3,1,1);
      pad1->SetBottomMargin(0.05);
      pad1->Draw();
      pad1->cd();

      hReco->SetMarkerStyle(21);
      //hMeas->SetMarkerStyle(25);

      hReco->GetXaxis()->SetRangeUser(400.,1200.);
      hTrue->GetXaxis()->SetRangeUser(400.,1200.);
      //hMeas->GetXaxis()->SetRangeUser(400.,1200.);

      TString xsec_title = ";;d#sigma/dp_{T} [fb/GeV]";
      if (doNormalize) xsec_title = ";;1/#sigma d#sigma/dp_{T} [1/GeV]";

      hReco->SetTitle(xsec_title);
      hReco->GetYaxis()->SetTitleOffset(1.2);
      hReco->SetMinimum(0.0);
      float max = hTrue->GetMaximum();
      float max2 = hReco->GetMaximum();
      if (max2 > max) max = max2;
      hReco->SetAxisRange(0,max*1.15,"Y");
      hReco->Draw();
      hTrue->Draw("hist same");
      //hMeas->Draw('same');
      hTrue->UseCurrentStyle();
      hTrue->SetLineColor(4);
      hTrue->GetYaxis()->SetTitleSize(25);
      hTrue->GetXaxis()->SetLabelSize(0);

      TLegend* leg = new TLegend(0.5, 0.55, 0.9, 0.75);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->SetBorderSize(0);

      if (isClosure){ 
	leg->AddEntry( hReco, "Unfolded data", "p");
	leg->AddEntry( hTrue, "Generated (Powheg)", "l");
	//leg->AddEntry( hMeas, "Raw data", "p");
      }
      else {
	leg->AddEntry( hReco, "Unfolded MC (Powheg)", "p");
	leg->AddEntry( hTrue, "Generated (Powheg)", "l");
	//leg->AddEntry( hMeas, "Reco-level (Powheg)", "p");
      }
    
      leg->Draw();

      TLatex* text1 = new TLatex();
      text1->SetNDC();
      text1->SetTextFont(42);
      text1->DrawLatex(0.55,0.8, "#scale[1.0]{L = 12.4 fb^{-1}, #sqrt{s} = 13 TeV}");

      c2->cd();
      TPad* pad2 = new TPad("pad2","pad2",0,0.0,1,0.28);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.4);
      pad2->Draw();
      pad2->cd();
      pad2->SetGridy();
      
      hFrac->SetMaximum(1.2 * hFrac->GetMaximum());
      hFrac->SetMinimum(0.8 * hFrac->GetMinimum());
      hFrac->UseCurrentStyle();
      hFrac->GetYaxis()->SetTitleSize(25);
      hFrac->GetYaxis()->SetTitleOffset(2.0);
      hFrac->GetXaxis()->SetTitleOffset(4.0);
      hFrac->GetXaxis()->SetLabelSize(25);
      hFrac->GetYaxis()->SetNdivisions(4,4,0,false);

      hFrac->Draw("e");
      hFrac->GetXaxis()->SetRangeUser(400., 1200.);

      c2->Update();
      c2->Print(TString::Format("UnfoldingPlots/unfolded_ttbar_xs_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      
      // --------------------------------------------------------------------------------------
      // Get uncertainty breakdown
      // --------------------------------------------------------------------------------------
      
      if (whichSyst != "none"){
	TH1F* h_sysUnc[9];
	for (int kk = 0; kk < 9; kk++){
	  h_sysUnc[kk] = (TH1F*) unfold.GetDeltaSysSource(syst_names[kk],"h_"+syst_names[kk]);
	  for (int ib = 1; ib < h_sysUnc[kk]->GetNbinsX()+1; ib++){
	    h_sysUnc[kk]->SetBinContent(ib,fabs((h_sysUnc[kk]->GetBinContent(ib) - hReco->GetBinContent(ib))/hReco->GetBinContent(ib) * 100.));
	  }
	}
	
	TH1F* h_bkgUnc;
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
	      h_bkgUnc_tmp[kk]->SetBinContent(ib,fabs((h_bkgUnc_tmp[kk]->GetBinContent(ib) - hReco->GetBinContent(ib))/hReco->GetBinContent(ib) * 100.));
	    }
	  }
	  
	  for (int ib = 1; ib < h_bkgUnc->GetNbinsX()+1; ib++){
	    h_bkgUnc->SetBinContent(ib,sqrt(pow(2,h_bkgUnc_tmp[0]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[1]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[2]->GetBinContent(ib))+
					    pow(2,h_bkgUnc_tmp[3]->GetBinContent(ib))));
	    
	  }
	}
	
	
	//TH1F* h_rel_unc_tau = (TH1F*) unfold.GetDeltaSysTau("h_tau");
	
	// --------------------------------------------------------------------------------------
	// Plot uncertainty breakdown
	// --------------------------------------------------------------------------------------
	
	TCanvas *c3 = new TCanvas("c3", "", 800, 600);
	c3->SetTopMargin(0.08);
	c3->SetRightMargin(0.05);
	c3->SetBottomMargin(0.14);
	c3->SetLeftMargin(0.16);
	
	if (toUnfold == "pt") h_sysUnc[0]->GetXaxis()->SetTitle("Top quark p_{T} (GeV)");
	else if (toUnfold == "y") h_sysUnc[0]->GetXaxis()->SetTitle("Top quark y");
	h_sysUnc[0]->GetYaxis()->SetTitle("Uncertainty [%]");
	if (toUnfold == "pt") h_sysUnc[0]->SetAxisRange(400,1150,"X");
	h_sysUnc[0]->SetAxisRange(0,100,"Y");
	
	int colors[9] = {632,600,617,417,432,801,864,906,419};
	int markers[9] = {20,21,22,23,33,24,25,26,32};
	for (int is = 0; is < 9; is++){
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
	
	TString longnames[9] = {"Lepton ID","Pileup","Jet energy scale","Jet energy resolution",
				"b tagging efficiency","t tagging efficiency",
				"PDF Uncertainty","#mu_{R}, #mu_{F} scales","#alpha_{s} scale"};
	
	TLegend* leg2;
	if (toUnfold == "pt" && doNormalize) leg2 = new TLegend(0.2,0.65,0.45,0.8);
	else leg2 = new TLegend(0.2,0.65,0.45,0.88);
	for (int is = 0; is < 5; is++){
	  leg2->AddEntry(h_sysUnc[is],longnames[is],"lp");
	}
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.04);
	leg2->SetTextFont(42);
	
	TLegend* leg22;
	if (toUnfold == "pt" && doNormalize) leg22 = new TLegend(0.55,0.65,0.8,0.8);
	else leg22 = new TLegend(0.55,0.65,0.8,0.88);
	for (int is = 5; is < 9; is++){
	  leg22->AddEntry(h_sysUnc[is],longnames[is],"lp");
	}
	leg22->AddEntry(h_bkgUnc,"Background normalization","lp");
	
	leg22->SetFillStyle(0);
	leg22->SetBorderSize(0);
	leg22->SetTextSize(0.04);
	leg22->SetTextFont(42);
	
	h_sysUnc[0]->Draw("hist");
	for (int is = 1; is < 9; is++){
	  h_sysUnc[is]->Draw("ep,same");
	}
	h_bkgUnc->Draw("ep,same");
	h_sysUnc[0]->Draw("hist,axis,same");
	leg2->Draw(); 
	leg22->Draw(); 

	c3->SaveAs(TString::Format("UnfoldingPlots/unfold_relative_uncertainties_"+toUnfold+"_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      }
      
	
      // --------------------------------------------------------------------------------------
      // Get covariance matrix
      // --------------------------------------------------------------------------------------
      /*
      TH2F* h_cov = (TH2F*) unfold.GetEmatrixTotal("h_covariance");
      
      // Normalize diagonal
      float diags[10] = {0.};
      for (int id = 0; id < h_cov->GetXaxis()->GetNbins(); id++){
	int idbin = h_cov->GetBin(id+1,id+1);
	diags[id] = sqrt(h_cov->GetBinContent(idbin));
      }
      for (int ix = 1; ix < h_cov->GetXaxis()->GetNbins()-1; ix++){
	for (int iy = 1; iy < h_cov->GetYaxis()->GetNbins()-1; iy++){
	  int ibin = h_cov->GetBin(ix+1,iy+1);
	  float rawval = h_cov->GetBinContent(ibin);
	  rawval *= 100.0 / diags[ix] / diags[iy];
	  int roundval = rawval * 1000;
	  float finalval = (float) roundval / 1000.;
	  h_cov->SetBinContent(ibin,finalval);
	}
      }
      
      TString binnames[5] = {"Bin 1","Bin 2","Bin 3","Bin 4","Bin 5"};
      for (int ii = 0; ii < h_cov->GetXaxis()->GetNbins() - 2; ii++){
	h_cov->GetXaxis()->SetBinLabel(ii+2,binnames[ii]);
	h_cov->GetYaxis()->SetBinLabel(ii+2,binnames[ii]);
      }
      
      TCanvas* c4 = new TCanvas("c4","c4",900,700);
      c4->SetLeftMargin(0.18);
      c4->SetRightMargin(0.12);
      c4->SetBottomMargin(0.1);
      h_cov->SetAxisRange(1.0,h_cov->GetXaxis()->GetNbins()-1.001,"X");
      h_cov->SetAxisRange(1.0,h_cov->GetXaxis()->GetNbins()-1.001,"Y");
      h_cov->GetZaxis()->SetLabelSize(0.04);
      h_cov->Draw("COLZ");
      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)h_cov->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.88);
      palette->SetX2NDC(0.93);
      palette->SetY1NDC(0.1);
      palette->SetY2NDC(0.9);
      gPad->Modified();
      gPad->Update();

      c4->SaveAs(TString::Format("UnfoldingPlots/total_covariance_"+toUnfold+"_"+channel+"_reg%d_scan%d.pdf",ii,jj));
      */

    }// End mode loop
  }// End regularization loop
}
