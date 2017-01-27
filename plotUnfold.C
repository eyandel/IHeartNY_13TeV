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

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize,char *text); 
void drawCMS(Double_t x1,Double_t y1, Double_t x2,Double_t y2, bool pad, bool prel, bool forPublic);

void symmetrize(TH1F* h_input);

void plot(TString channel, TString toUnfold, bool doNormalized, bool doLogscale, bool doAverageErr);

void plotUnfold(TString channel, TString toUnfold="pt") {
  
  bool doNormalized = false;  // normalized differential?
  bool doLogscale   = true;   // plot distributions on log scale?
  if (toUnfold == "y") doLogscale = false;
  bool doAverageErr = true;   // average up/down systematic uncertainties? 
  
  plot(channel,toUnfold,doNormalized,doLogscale,doAverageErr);
  
}

void plot(TString channel, TString toUnfold, bool doNormalized, bool doLogscale, bool doAverageErr) {
  
  SetPlotStyle();
  
  TString normflag = "";
  if (doNormalized) normflag = "_norm";
  
  // ---------------------------------------------------------------------------------------------------------------
  // get files & histograms
  // ---------------------------------------------------------------------------------------------------------------
  
  cout << endl << "Getting files and hists..." << endl;
  
  const int nSYSTEXP = 7;
  const int nSYSTTH = 3;
  const int nSYST = nSYSTEXP + nSYSTTH;
  TString name_syst[2*nSYST+1] = {"nom", "lepUp", "lepDown", "puUp", "puDown", "JECUp", "JECDown",
				  "JERUp", "JERDown", "BTagUp", "BTagDown", "TopTagUp", "TopTagDown",
				  "BkgUp", "BkgDown", "PDFUp", "PDFDown", "Q2Up", "Q2Down", "ASUp", "ASDown"};
				  
  TString sysnames[nSYST] = {"lep","pu","JEC","JER","BTag","TopTag","Bkg","PDF","Q2","AS"};
				  
  if (channel == "comb") const int nCHANNEL = 2;
  else const int nCHANNEL = 1;
  
  TString channel_name[nCHANNEL];
  if (nCHANNEL == 1) channel_name[0] = channel;
  else if (nCHANNEL == 2) {
    channel_name[0] = "muon";
    channel_name[1] = "ele";
  }
  else return; //shouldn't happen
  
  TFile* f_syst[nCHANNEL][2*nSYST+1];
  TH1F* h_true_tmp[nCHANNEL];
  TH1F* h_unfolded_tmp[nCHANNEL][2*nSYST+1];
  
  for (int ic=0; ic<nCHANNEL; ic++) {
    for (int is=0; is<2*nSYST+1; is++) {
      
      f_syst[ic][is] = new TFile("UnfoldingPlots/unfold_"+toUnfold+"_PowhegPythia8_"+channel_name[ic]+"_"+name_syst[is]+normflag+".root");
      
      h_unfolded_tmp[ic][is] = (TH1F*) f_syst[ic][is]->Get("UnfoldedData")->Clone();
      h_unfolded_tmp[ic][is]->SetName("UnfoldedData"+name_syst[is]+"_"+channel_name[ic]);
      if (toUnfold == "pt") h_unfolded_tmp[ic][is]->SetAxisRange(400,1150,"X");
      
      if (is==0) {
	h_true_tmp[ic] = (TH1F*) f_syst[ic][is]->Get(toUnfold+"GenTop_true")->Clone();
	h_true_tmp[ic]->SetName(toUnfold+"_genTop_"+channel_name[ic]);
	if (toUnfold == "pt") h_true_tmp[ic]->SetAxisRange(400,1150,"X");
      }          
    }//end systematic loop
  }//end channel loop
  
  // ----------------------------------------------------------------------------------------------------------------
  // possibly perform combination of electron+muon channels
  // ----------------------------------------------------------------------------------------------------------------
  
  // histograms that need to be combined:
  // h_true_tmp[ic] --> top pt @ parton-level
  // h_unfolded_tmp[ic][is] --> data unfolded to parton-level
  
  TH1F* h_true;
  TH1F* h_unfolded[2*nSYST+1];
  
  // ----------------------------------------------------------------------------------------------------------------
  // no combination -- just copy the electron/muon histograms
  if (channel != "comb") {
    
    h_true = (TH1F*) h_true_tmp[0]->Clone();
    h_true->SetName(toUnfold+"_genTop");
    
    for (int is=0; is<2*nSYST+1; is++) {
      h_unfolded[is] = (TH1F*) h_unfolded_tmp[0][is]->Clone();
      h_unfolded[is]->SetName("UnfoldedData"+name_syst[is]);
    }
    
  }
  // ----------------------------------------------------------------------------------------------------------------
  // combination of the channels
  else {
    
    int nBIN = h_true_tmp[0]->GetNbinsX();
    
    // first create empty histograms
    h_true = (TH1F*) h_true_tmp[0]->Clone();
    h_true->SetName(toUnfold+"_genTop");
    h_true->Reset();
    
    for (int is=0; is<2*nSYST+1; is++) {
      h_unfolded[is] = (TH1F*) h_unfolded_tmp[0][is]->Clone();
      h_unfolded[is]->SetName("UnfoldedData"+name_syst[is]);
      h_unfolded[is]->Reset();
    }
    
    //then do the combination
    float nel = 0;
    float nmu = 0;
    float ncomb = 0;
    float snel = 0;
    float snmu = 0;
    float sncomb = 0;
    
    for (int ib=0; ib<nBIN; ib++) {
      
      // h_true
      nmu = h_true_tmp[0]->GetBinContent(ib+1);
      nel = h_true_tmp[1]->GetBinContent(ib+1);
      snmu = h_true_tmp[0]->GetBinError(ib+1);
      snel = h_true_tmp[1]->GetBinError(ib+1);
      if (snel==0 || snmu==0) {
	ncomb = 0;
	sncomb = 0;
      }
      else {
	ncomb = (nel/(snel*snel) + nmu/(snmu*snmu)) / (1.0/(snel*snel) + 1.0/(snmu*snmu));
	sncomb = 1.0 / (1.0/(snel*snel) + 1.0/(snmu*snmu));
	sncomb = sqrt(sncomb);
      }
      h_true->SetBinContent(ib+1,ncomb);
      h_true->SetBinError(ib+1,sncomb);
      
      for (int is=0; is<2*nSYST+1; is++) {
	// h_unfolded
	nmu = h_unfolded_tmp[0][is]->GetBinContent(ib+1);
	nel = h_unfolded_tmp[1][is]->GetBinContent(ib+1);
	snmu = h_unfolded_tmp[0][is]->GetBinError(ib+1);
	snel = h_unfolded_tmp[1][is]->GetBinError(ib+1);
	if (snel==0 || snmu==0) {
	  ncomb = 0;
	  sncomb = 0;
	}
	else {
	  ncomb = (nel/(snel*snel) + nmu/(snmu*snmu)) / (1.0/(snel*snel) + 1.0/(snmu*snmu));
	  sncomb = 1.0 / (1.0/(snel*snel) + 1.0/(snmu*snmu));
	  sncomb = sqrt(sncomb);
	}
	h_unfolded[is]->SetBinContent(ib+1,ncomb);
	h_unfolded[is]->SetBinError(ib+1,sncomb);
	
      }//end loop systematics
    }//end bin loop
  }//end do combination
  
  // Delete tmp hists once no longer needed
  for (int ic=0; ic<nCHANNEL; ic++) {
    h_true_tmp[ic]->Delete();
    for (int is=0; is<2*nSYST+1; is++) {
      h_unfolded_tmp[ic][is]->Delete();
    }
  }
  
  // --------------------------------------------------------------------------------------
  
  TH1F* h_dummy = (TH1F*) h_unfolded[0]->Clone("dummy");
  h_dummy->Reset();
  
  TH1F* h_dummy_r = (TH1F*) h_unfolded[0]->Clone("dummy_r");
  h_dummy_r->Reset();
  
  if (doNormalized) {
    if (toUnfold == "pt"){
      h_dummy->GetYaxis()->SetTitle("1/#sigma d#sigma/dp_{T} (1/GeV)");
      //h_dummy->SetAxisRange(0,0.008,"Y");
    }
    else if (toUnfold == "y"){
      h_dummy->GetYaxis()->SetTitle("1/#sigma d#sigma/dy");
      //h_dummy->SetAxisRange(0,0.5,"Y");
    }
  }
  else {
    if (toUnfold == "pt"){
      h_dummy->GetYaxis()->SetTitle("d#sigma/dp_{T} (fb/GeV)");
      //h_dummy->SetAxisRange(0,12,"Y");
    }
    else if (toUnfold == "y"){
      h_dummy->GetYaxis()->SetTitle("d#sigma/dy (fb)");
      //h_dummy->SetAxisRange(0,900,"Y");
    }
  }
  
  // ----------------------------------------------------------------------------------------------------------------
  // colors and stuff
  // ----------------------------------------------------------------------------------------------------------------
  
  h_true->SetLineColor(2);
  h_true->SetLineWidth(3);
  
  h_unfolded[0]->SetLineColor(1);
  h_unfolded[0]->SetLineWidth(2);
  h_unfolded[0]->SetMarkerColor(1);
  h_unfolded[0]->SetMarkerStyle(8);
  
  float tmp_max = 0;
  
  // ----------------------------------------------------------------------------------------------------------------
  // systematics for ratio plot
  // ----------------------------------------------------------------------------------------------------------------
  
  TH1F* h_stat_up = (TH1F*) h_true->Clone("stat");
  TH1F* h_stat_dn = (TH1F*) h_true->Clone("stat");
  h_stat_up->Reset();
  h_stat_dn->Reset();
  
  TH1F* h_systEXP_up = (TH1F*) h_true->Clone("syst_up_exp");
  TH1F* h_systEXP_dn = (TH1F*) h_true->Clone("syst_dn_exp");
  h_systEXP_up->Reset();
  h_systEXP_dn->Reset();
  
  TH1F* h_systTH_up = (TH1F*) h_true->Clone("syst_up_th");
  TH1F* h_systTH_dn = (TH1F*) h_true->Clone("syst_dn_th");
  h_systTH_up->Reset();
  h_systTH_dn->Reset();
  
  TH1F* h_systTOT_up = (TH1F*) h_true->Clone("syst_up_tot");
  TH1F* h_systTOT_dn = (TH1F*) h_true->Clone("syst_dn_tot");
  h_systTOT_up->Reset();
  h_systTOT_dn->Reset();
  
  // experimental+theory
  TH1F* h_syst[nSYST];
  TH2F* h_cov[nSYST];
  const int nBINS = h_true->GetNbinsX();
  for (int ii = 0; ii < nSYST; ii++){
    h_syst[ii] = (TH1F*) h_true->Clone("syst_"+sysnames[ii]);
    h_syst[ii]->Reset();
    h_cov[ii] = new TH2F("cov_"+sysnames[ii],"",nBINS,0,nBINS,nBINS,0,nBINS);
  }
  
  // stat, lumi
  TH1F* h_syst_stat   = (TH1F*) h_true->Clone("syst_stat");
  h_syst_stat->Reset();
  TH1F* h_lumi        = (TH1F*) h_true->Clone("lumi");
  h_lumi->Reset();
  
  TH1F* h_syst_tot    = (TH1F*) h_true->Clone("syst_tot");
  h_syst_tot->Reset();
  
  float count[2*nSYST+1] = {0};
  float sig[2*nSYST+1] = {0};
  
  float xsec_meas = 0;
  float xsec_true = 0;
  
  // ----------------------------------------------------------------------------------------------------------------
  // loop over bins
  // ----------------------------------------------------------------------------------------------------------------
  
  cout << endl << "*** unfolding result (parton-level) for " << channel << " channel ***" << endl;
  
  for (int i=0; i<h_unfolded[0]->GetNbinsX(); i++) {

    for (int is=0; is<2*nSYST+1; is++) {
      count[is] = h_unfolded[is]->GetBinContent(i+1);
      sig[is]   = h_unfolded[is]->GetBinError(i+1);
    }
    float truth = h_true->GetBinContent(i+1);
    float this_lumi = 0.027 * count[0];
    
    // experimental
    float this_systEXP_up = 0;
    float this_systEXP_dn = 0;
    for (int ie = 0; ie < nSYSTEXP; ie++){
      this_systEXP_up += (count[2*ie+1]-count[0])*(count[2*ie+1]-count[0]); //jec
      this_systEXP_dn += (count[2*ie+2]-count[0])*(count[2*ie+2]-count[0]);
    }
    this_systEXP_up = sqrt(this_systEXP_up);
    this_systEXP_dn = sqrt(this_systEXP_dn);
    
    // theory
    float this_systTH_up = 0;
    float this_systTH_dn = 0;
    for (int it = 0; it < nSYSTTH; it++){
      this_systTH_up += (count[2*it+15]-count[0])*(count[2*it+15]-count[0]);
      this_systTH_dn += (count[2*it+16]-count[0])*(count[2*it+16]-count[0]);
    }
    this_systTH_up = sqrt(this_systTH_up);
    this_systTH_dn = sqrt(this_systTH_dn);

    // theory + experimental + statistical (+ lumi)
    float this_systTOT_up = this_systEXP_up*this_systEXP_up + this_systTH_up*this_systTH_up + sig[0]*sig[0];
    if (doNormalized == false) this_systTOT_up += this_lumi*this_lumi;
    this_systTOT_up = sqrt(this_systTOT_up);
    float this_systTOT_dn = this_systEXP_dn*this_systEXP_dn + this_systTH_dn*this_systTH_dn + sig[0]*sig[0];
    if (doNormalized == false) this_systTOT_dn += this_lumi*this_lumi;
    this_systTOT_dn = sqrt(this_systTOT_dn);

    // fill histograms for ratio plot
    h_stat_up->SetBinContent(i+1,count[0]+sig[0]);
    h_stat_dn->SetBinContent(i+1,count[0]-sig[0]);

    if (doAverageErr) {
      h_systEXP_up->SetBinContent(i+1,count[0]+(this_systEXP_up+this_systEXP_dn)/2);
      h_systEXP_dn->SetBinContent(i+1,count[0]-(this_systEXP_up+this_systEXP_dn)/2);   
      h_systTH_up->SetBinContent(i+1,count[0]+(this_systTH_up+this_systTH_dn)/2);
      h_systTH_dn->SetBinContent(i+1,count[0]-(this_systTH_up+this_systTH_dn)/2);  
      h_systTOT_up->SetBinContent(i+1,count[0]+(this_systTOT_up+this_systTOT_dn)/2);
      h_systTOT_dn->SetBinContent(i+1,count[0]-(this_systTOT_up+this_systTOT_dn)/2);  
    }
    else {
      h_systEXP_up->SetBinContent(i+1,count[0]+this_systEXP_up);
      h_systEXP_dn->SetBinContent(i+1,count[0]-this_systEXP_dn);   
      h_systTH_up->SetBinContent(i+1,count[0]+this_systTH_up);
      h_systTH_dn->SetBinContent(i+1,count[0]-this_systTH_dn);  
      h_systTOT_up->SetBinContent(i+1,count[0]+this_systTOT_up);
      h_systTOT_dn->SetBinContent(i+1,count[0]-this_systTOT_dn);  
    }

    // histograms for uncertainty plots
   
    // experimental+systematic
    float systs[2*nSYST];
    for (int is = 0; is < 2*nSYST; is++){
      systs[is] = fabs((count[is+1]-count[0])/count[0])*100;
    }
    
    // statistics    
    float syst_stat   = sig[0]/count[0]*100;

    float max_systs[nSYST] = {0};
    for (int im = 0; im < nSYST; im++){
      if (doAverageErr) max_systs[im] = (systs[2*im]+systs[2*im+1])/2;
      else max_systs[im] = max(systs[2*im],systs[2*im+1]);
    }

    float syst_totalEXP_up = sqrt(systs[0]*systs[0] + systs[2]*systs[2] + systs[4]*systs[4] + systs[6]*systs[6] + systs[8]*systs[8] + systs[10]*systs[10]
				   + systs[12]*systs[12]);
    float syst_totalEXP_dn = sqrt(systs[1]*systs[1] + systs[3]*systs[3] + systs[5]*systs[5] + systs[7]*systs[7] + systs[9]*systs[9] + systs[11]*systs[11]
				   + systs[13]*systs[13]);
    float syst_totalTH_up = sqrt(systs[14]*systs[14] + systs[16]*systs[16]+systs[18]*systs[18]);
    float syst_totalTH_dn = sqrt(systs[15]*systs[15] + systs[17]*systs[17]+systs[19]*systs[19]);
    float syst_totalSYS_up = sqrt(syst_totalEXP_up*syst_totalEXP_up + syst_totalTH_up*syst_totalTH_up);
    float syst_totalSYS_dn = sqrt(syst_totalEXP_dn*syst_totalEXP_dn + syst_totalTH_dn*syst_totalTH_dn);

    float syst_total_up = sqrt(syst_totalSYS_up*syst_totalSYS_up + syst_stat*syst_stat);
    float syst_total_dn = sqrt(syst_totalSYS_dn*syst_totalSYS_dn + syst_stat*syst_stat);

    // Add lumi uncertainty
    if (doNormalized == false) {
      syst_totalSYS_up = sqrt(syst_totalSYS_up*syst_totalSYS_up+2.7*2.7);
      syst_total_up    = sqrt(syst_total_up*syst_total_up+2.7*2.7);
      syst_totalSYS_dn = sqrt(syst_totalSYS_dn*syst_totalSYS_dn+2.7*2.7);
      syst_total_dn    = sqrt(syst_total_dn*syst_total_dn+2.7*2.7);
    }

    float max_syst_totalEXP = 0;
    float max_syst_totalTH = 0;
    float max_syst_totalSYS = 0;
    float max_syst_total = 0;
    
    if (doAverageErr) max_syst_totalEXP = (syst_totalEXP_up+syst_totalEXP_dn)/2.0;
    if (doAverageErr) max_syst_totalTH  = (syst_totalTH_up+syst_totalTH_dn)/2.0;
    if (doAverageErr) max_syst_totalSYS = (syst_totalSYS_up+syst_totalSYS_dn)/2.0;
    if (doAverageErr) max_syst_total    = (syst_total_up+syst_total_dn)/2.0;
    
    if (!doAverageErr) max_syst_totalEXP = max(syst_totalEXP_up,syst_totalEXP_dn);
    if (!doAverageErr) max_syst_totalTH  = max(syst_totalTH_up,syst_totalTH_dn);
    if (!doAverageErr) max_syst_totalSYS = max(syst_totalSYS_up,syst_totalSYS_dn);
    if (!doAverageErr) max_syst_total    = max(syst_total_up,syst_total_dn);

    for (int is = 0; is < nSYST; is++){
      h_syst[is]->SetBinContent(i+1,max_systs[is]);
      h_syst[is]->SetBinError(i+1,0.001);
    }
    h_syst_stat->SetBinContent(i+1,syst_stat);  
    h_syst_stat->SetBinError(i+1,0.001);
    h_syst_tot->SetBinContent(i+1,max_syst_total);
    h_syst_tot->SetBinError(i+1,0.001);
    if (doNormalized == false) {
      h_lumi->SetBinContent(i+1,2.7);
      h_lumi->SetBinError(i+1,0.001);
    }

    // ----------------------------------------------------------------------------------------------------------------
    // print outs for cross-section table
    // ----------------------------------------------------------------------------------------------------------------

    int ibin = h_systEXP_up->GetBin(i+1);
    float lowedge = h_systEXP_up->GetXaxis()->GetBinLowEdge(ibin);
    float highedge = h_systEXP_up->GetXaxis()->GetBinUpEdge(ibin);

    if (toUnfold == "pt"){
      if (lowedge > 300 && highedge < 1300) {
	cout << (float)lowedge << "--" << (float)highedge << " & " << count[0] << " & " << syst_stat << " & " 
	     << max_syst_totalEXP << " & " << max_syst_totalTH << " & " << max_syst_total << " & " 
	     << h_true->GetBinContent(i+1) << endl;
      }
      
      if (lowedge > 300.) {
	xsec_meas += count[0]*h_true->GetBinWidth(i+1);
	xsec_true += h_true->GetBinContent(i+1)*h_true->GetBinWidth(i+1);
      }   
    }

    else if (toUnfold == "y"){
      cout << (float)lowedge << "--" << (float)highedge << " & " << count[0] << " & " << syst_stat << " & " 
	   << max_syst_totalEXP << " & " << max_syst_totalTH << " & " << max_syst_total << " & " 
	   << h_true->GetBinContent(i+1) << endl;
      
      xsec_meas += count[0]*h_true->GetBinWidth(i+1);
      xsec_true += h_true->GetBinContent(i+1)*h_true->GetBinWidth(i+1);
    }
  }
    
  cout << "TOTAL cross section (parton level), measured = " << xsec_meas << " true = " << xsec_true << endl;

  // Symmetrize uncertainties for rapidity
  if (toUnfold == "y"){
    for (int is = 0; is < nSYST; is++) symmetrize(h_syst[is]);
    symmetrize(h_syst_stat);
    symmetrize(h_syst_tot);
  }

  // ----------------------------------------------------------------------------------------------------------------
  // Covariance matrices
  // ----------------------------------------------------------------------------------------------------------------
  TH2F* h_cov_tot = (TH2F*) h_cov[0]->Clone("h_cov_tot");
  h_cov_tot->Reset();
  
  for (int is = 0; is < nSYST; is++){
    for (int ix = 0; ix < nBINS; ix++){
      for (int iy = 0; iy < nBINS; iy++){
	
	float cov = 0.0;
	
	if (h_unfolded[0]->GetXaxis()->GetBinLowEdge(ix+1) >= 400. && h_unfolded[0]->GetXaxis()->GetBinLowEdge(iy+1) >= 400. &&
	    h_unfolded[0]->GetXaxis()->GetBinLowEdge(ix+1) < 1199. && h_unfolded[0]->GetXaxis()->GetBinLowEdge(iy+1) < 1199.) {
	  float dxmax = max(fabs(h_unfolded[2*is+1]->GetBinContent(ix+1) - h_unfolded[0]->GetBinContent(ix+1)),
			    fabs(h_unfolded[2*is+2]->GetBinContent(ix+1) - h_unfolded[0]->GetBinContent(ix+1)));
	  float dymax =  max(fabs(h_unfolded[2*is+1]->GetBinContent(iy+1) - h_unfolded[0]->GetBinContent(iy+1)),
			     fabs(h_unfolded[2*is+2]->GetBinContent(iy+1) - h_unfolded[0]->GetBinContent(iy+1)));
	  float signx = (h_unfolded[2*is+1]->GetBinContent(ix+1)-h_unfolded[2*is+2]->GetBinContent(ix+1))
	    / fabs(h_unfolded[2*is+1]->GetBinContent(ix+1)-h_unfolded[2*is+2]->GetBinContent(ix+1));
	  float signy = (h_unfolded[2*is+1]->GetBinContent(iy+1)-h_unfolded[2*is+2]->GetBinContent(iy+1))
	    / fabs(h_unfolded[2*is+1]->GetBinContent(iy+1)-h_unfolded[2*is+2]->GetBinContent(iy+1));
	
	  cov = dxmax * dymax * signx * signy;
	}

	h_cov[is]->SetBinContent(ix+1,iy+1,cov);
      }
    }
    h_cov_tot->Add(h_cov[is]);
  }

  // ----------------------------------------------------------------------------------------------------------------
  // RATIO PLOTS !!!!!!!
  // ----------------------------------------------------------------------------------------------------------------
  
  Float_t f_textSize=0.057;
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetTextSize(f_textSize);
  gStyle->SetLabelSize(f_textSize,"x");
  gStyle->SetTitleSize(f_textSize,"x");
  gStyle->SetLabelSize(f_textSize,"y");
  gStyle->SetTitleSize(f_textSize,"y");
  gStyle->SetOptStat(0);

  TH1F* h_fullunc = (TH1F*) h_unfolded[0]->Clone("fullunc");

  for (int ib=0; ib<h_fullunc->GetNbinsX(); ib++) {
    if (h_fullunc->GetBinContent(ib+1) > 0) {
      float up = h_systTOT_up->GetBinContent(ib+1);
      float dn = h_systTOT_dn->GetBinContent(ib+1);
      h_fullunc->SetBinError(ib+1,(up-dn)/2);
    }
  }

  gStyle->SetEndErrorSize(5);

  TCanvas *c = new TCanvas("c", "", 700, 625);
  
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

  /*
  if (doLogscale && doNormalized) {
    h_dummy->SetAxisRange(0.000015,0.015,"Y");
    h_true->SetAxisRange(0.000015,0.015,"Y");
    h_trueMG->SetAxisRange(0.000015,0.015,"Y");
    h_trueMCNLO->SetAxisRange(0.000015,0.015,"Y");
    h_unfolded[0]->SetAxisRange(0.000015,0.015,"Y");
    h_fullunc->SetAxisRange(0.000015,0.015,"Y");
    gPad->SetLogy();
  }
  else if (doLogscale) {
    if (toUnfold == "y") {
      h_dummy->SetAxisRange(1,3000,"Y");
      h_true->SetAxisRange(1,3000,"Y");
      h_trueMG->SetAxisRange(1,3000,"Y");
      h_trueMCNLO->SetAxisRange(1,3000,"Y");
      h_unfolded[0]->SetAxisRange(1,3000,"Y");
      h_fullunc->SetAxisRange(1,3000,"Y");
      gPad->SetLogy();
    }
    else {
      h_dummy->SetAxisRange(0.015,30,"Y");
      h_true->SetAxisRange(0.015,30,"Y");
      h_trueMG->SetAxisRange(0.015,30,"Y");
      h_trueMCNLO->SetAxisRange(0.015,30,"Y");
      h_unfolded[0]->SetAxisRange(0.015,30,"Y");
      h_fullunc->SetAxisRange(0.015,30,"Y");
      gPad->SetLogy();
    }
  }
  */

  h_dummy->GetYaxis()->SetTitleSize(0.075);    
  if (doLogscale) h_dummy->GetYaxis()->SetTitleOffset(0.8);
  else if (toUnfold == "y") h_dummy->GetYaxis()->SetTitleOffset(0.9);
  else h_dummy->GetYaxis()->SetTitleOffset(1.1);
  h_dummy->GetYaxis()->SetLabelSize(0.065);
  h_dummy->GetXaxis()->SetTitleSize(0);
  h_dummy->GetXaxis()->SetLabelSize(0);

  h_dummy->Draw("hist");
  h_true->Draw("hist,same");
  TExec *setex1 = new TExec("setex1","gStyle->SetErrorX(0)");
  setex1->Draw();
  h_fullunc->Draw("hist,ep,same");
  h_unfolded[0]->Draw("hist,E1p,same");
  TExec *setex2 = new TExec("setex2","gStyle->SetErrorX(0.5)");
  setex2->Draw();
  h_dummy->Draw("hist,axis,same"); 
  
  // ----------------------------------------------------------------------------------------------------------------
  // making ratio part of plot
  // ----------------------------------------------------------------------------------------------------------------
  
  TString baselab = h_unfolded[0]->GetName();
  TH1F* h_ratio = (TH1F*) h_unfolded[0]->Clone(baselab+"_ratio");

  TH1F* h_ratioSTAT_up = (TH1F*) h_stat_up->Clone(baselab+"_ratioSTAT_up");
  TH1F* h_ratioSTAT_dn = (TH1F*) h_stat_dn->Clone(baselab+"_ratioSTAT_dn");

  TH1F* h_ratioEXP_up = (TH1F*) h_systEXP_up->Clone(baselab+"_ratioEXP_up");
  TH1F* h_ratioEXP_dn = (TH1F*) h_systEXP_dn->Clone(baselab+"_ratioEXP_dn");
  
  TH1F* h_ratioTOT_up = (TH1F*) h_systTOT_up->Clone(baselab+"_ratioTOT_up");
  TH1F* h_ratioTOT_dn = (TH1F*) h_systTOT_dn->Clone(baselab+"_ratioTOT_dn");

  TH1F* h_ratioTH_up = (TH1F*) h_systTH_up->Clone(baselab+"_ratioTH_up");
  TH1F* h_ratioTH_dn = (TH1F*) h_systTH_dn->Clone(baselab+"_ratioTH_dn");

  TH1F* h_ratioGEN = (TH1F*) h_true->Clone(baselab+"_ratioGEN");

  h_ratio->Divide(h_unfolded[0]);
  h_ratioSTAT_up->Divide(h_unfolded[0]);
  h_ratioSTAT_dn->Divide(h_unfolded[0]);
  h_ratioEXP_up->Divide(h_unfolded[0]);
  h_ratioEXP_dn->Divide(h_unfolded[0]);
  h_ratioTH_up->Divide(h_unfolded[0]);
  h_ratioTH_dn->Divide(h_unfolded[0]);
  h_ratioTOT_up->Divide(h_unfolded[0]);
  h_ratioTOT_dn->Divide(h_unfolded[0]);
  h_ratioGEN->Divide(h_unfolded[0]);
  
  // ----------------------------------------------------------------------------------------------------------------
  TH1F* blaSTAT = (TH1F*) h_ratioSTAT_up->Clone("blaSTAT");
  blaSTAT->Reset();
  for (int i=0; i<blaSTAT->GetNbinsX(); i++) {
    float up = h_ratioSTAT_up->GetBinContent(i+1);
    float dn = h_ratioSTAT_dn->GetBinContent(i+1);
    
    blaSTAT->SetBinContent(i+1,(up-dn)/2+dn);
    blaSTAT->SetBinError(i+1,(up-dn)/2);
  }
  
  TH1F* blaEXP = (TH1F*) h_ratioEXP_up->Clone("blaEXP");
  blaEXP->Reset();
  for (int i=0; i<blaEXP->GetNbinsX(); i++) {
    float up = h_ratioEXP_up->GetBinContent(i+1);
    float dn = h_ratioEXP_dn->GetBinContent(i+1);
    
    blaEXP->SetBinContent(i+1,(up-dn)/2+dn);
    blaEXP->SetBinError(i+1,(up-dn)/2);
  }
  
  TH1F* blaTH = (TH1F*) h_ratioTH_up->Clone("blaTH");
  blaTH->Reset();
  for (int i=0; i<blaTH->GetNbinsX(); i++) {
    float up = h_ratioTH_up->GetBinContent(i+1);
    float dn = h_ratioTH_dn->GetBinContent(i+1);
    
    blaTH->SetBinContent(i+1,(up-dn)/2+dn);
    blaTH->SetBinError(i+1,(up-dn)/2);
  }

  TH1F* blaTOT = (TH1F*) h_ratioTOT_up->Clone("blaTOT");
  blaTOT->Reset();
  for (int i=0; i<blaTOT->GetNbinsX(); i++) {
    float up = h_ratioTOT_up->GetBinContent(i+1);
    float dn = h_ratioTOT_dn->GetBinContent(i+1);
    
    blaTOT->SetBinContent(i+1,(up-dn)/2+dn);
    blaTOT->SetBinError(i+1,(up-dn)/2);
  }
  
  blaSTAT->SetMarkerSize(0);
  blaSTAT->SetLineColor(0);
  blaSTAT->SetLineWidth(0);
  blaSTAT->SetFillColor(18);
  blaSTAT->SetFillStyle(1001);

  blaEXP->SetMarkerSize(0);
  blaEXP->SetLineColor(0);
  blaEXP->SetFillColor(kAzure);
  blaEXP->SetFillStyle(1001);

  blaTH->SetMarkerSize(0);
  blaTH->SetLineColor(0);
  blaTH->SetFillColor(kBlue);
  blaTH->SetFillStyle(3453);

  blaTOT->SetMarkerSize(0);
  blaTOT->SetLineColor(0);
  blaTOT->SetFillColor(kAzure-9);
  blaTOT->SetFillStyle(1001);

  if (toUnfold == "y"){
    symmetrize(blaTOT);
    symmetrize(blaSTAT);
  }

  // ----------------------------------------------------------------------------------------------------------------
  
  float legxlow = 0.48;
  float legxhigh = 0.81;
  float legylow = 0.5;
  float legyhigh = 0.85;
  if (toUnfold == "y"){
    legxlow = 0.37;
    legylow = 0.05;
    legxhigh = 0.70;
    legyhigh = 0.38;
  }
  
  TLegend* leg = new TLegend(legxlow,legylow,legxhigh,legyhigh);  
  leg->AddEntry(h_unfolded[0],"Data","pe1");
  leg->AddEntry(h_true,"Powheg+Pythia8","l");
  leg->AddEntry(blaSTAT,"Stat. uncertainty","f");
  leg->AddEntry(blaTOT,"Stat. #oplus syst. uncertainties","f");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if (toUnfold == "pt") leg->SetTextSize(0.055);
  else if (toUnfold == "y") leg->SetTextSize(0.05);
  leg->SetTextFont(42);
  leg->Draw();

  if (toUnfold=="y") drawCMS(0.21,0.79,0.21,0.72,true,true,false);
  else drawCMS(0.36,0.79,0.47,0.79,true,true,false);

  // ----------------------------------------------------------------------------------------------------------------

  p2->cd();

  h_ratio->GetXaxis()->SetTitleSize(0.15);
  h_ratio->GetYaxis()->SetTitleSize(0.13);
  h_ratio->GetXaxis()->SetLabelSize(0.12);
  h_ratio->GetYaxis()->SetLabelSize(0.11);

  h_ratio->GetYaxis()->SetTitle("Theory / Data");
  if (toUnfold == "pt") h_ratio->GetXaxis()->SetTitle("Top quark p_{T} (GeV)");
  else if (toUnfold == "y") h_ratio->GetXaxis()->SetTitle("Top quark y");
  h_ratio->GetYaxis()->SetTitleOffset(0.8);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetYaxis()->SetTitleOffset(0.38);
  h_ratio->GetXaxis()->SetTitleOffset(1.0);

  h_ratio->SetAxisRange(0.0,2.0,"Y");
  h_ratioGEN->SetAxisRange(0.0,2.0,"Y");

  h_ratio->Draw("hist");
  blaTOT->Draw("same,e2");
  //blaEXP->Draw("same,e2");
  //blaTH->Draw("same,e2");
  //blaSTAT->Draw("same,ep");
  
  blaSTAT->Draw("same,e2");
  h_ratioGEN->Draw("same,hist");
  h_ratio->Draw("same,hist");
  h_ratio->Draw("same,axis");
  
  c->SaveAs("UnfoldingPlots/unfoldWithError_"+toUnfold+"_"+channel+normflag+".pdf");
  
  // ----------------------------------------------------------------------------------------------
  // Uncertainties plots

  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->SetTopMargin(0.08);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.14);
  c1->SetLeftMargin(0.16);

  if (toUnfold == "pt") h_dummy_r->GetXaxis()->SetTitle("Top quark p_{T} (GeV)");
  else if (toUnfold == "y") h_dummy_r->GetXaxis()->SetTitle("Top quark y");
  h_dummy_r->GetYaxis()->SetTitle("Uncertainty [%]");
  if (toUnfold == "pt") h_dummy_r->SetAxisRange(400,1150,"X");
  h_dummy_r->SetAxisRange(0,100,"Y");
  
  h_dummy_r->GetYaxis()->SetTitleSize(0.055);    
  h_dummy_r->GetYaxis()->SetTitleOffset(1.1);
  h_dummy_r->GetYaxis()->SetLabelSize(0.045);
  
  h_dummy_r->GetXaxis()->SetTitleSize(0.05);
  h_dummy_r->GetXaxis()->SetTitleOffset(1.2);
  h_dummy_r->GetXaxis()->SetLabelSize(0.0455);
  
  c1->cd();

  h_syst_tot->SetFillColor(17);
  h_syst_tot->SetFillStyle(3344);
  h_syst_tot->SetLineColor(16);
  h_syst_tot->SetLineWidth(2);

  h_syst_stat->SetLineColor(1);
  h_syst_stat->SetLineWidth(2);
  h_syst_stat->SetMarkerColor(1);
  h_syst_stat->SetMarkerStyle(20);

  int colors[nSYST] = {632,600,617,417,432,801,864,906,419,426};
  int markers[nSYST] = {20,21,22,23,33,24,25,26,32,27};
  for (int is = 0; is < nSYST; is++){
    h_syst[is]->SetLineColor(colors[is]);
    h_syst[is]->SetLineWidth(2);
    h_syst[is]->SetMarkerColor(colors[is]);
    h_syst[is]->SetMarkerStyle(markers[is]);
  }

  if (doNormalized == false){
    h_lumi->SetLineColor(40);
    h_lumi->SetLineWidth(2);
    h_lumi->SetMarkerColor(40);
    h_lumi->SetMarkerStyle(34);
  }

  TString longnames[nSYST] = {"Lepton ID","Pileup","Jet energy scale","Jet energy resolution",
			      "b tagging efficiency","t tagging efficiency","Background normalization",
			      "PDF Uncertainty","#mu_{R}, #mu_{F} scales","#alpha_{s} scale"};
  
  TLegend* leg2;
  if (toUnfold == "pt" && doNormalized) leg2 = new TLegend(0.2,0.65,0.45,0.8);
  else leg2 = new TLegend(0.2,0.65,0.45,0.88);
  leg2->AddEntry(h_syst_tot,"Total syst. uncertainty","f");
  leg2->AddEntry(h_syst_stat,"Statistical uncertainty","lp");
  if (doNormalized == false) leg2->AddEntry(h_lumi,"Int. luminosity","lp");
  for (int is = 0; is < 4; is++){
    leg2->AddEntry(h_syst[is],longnames[is],"lp");
  }
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);
  leg2->SetTextFont(42);

  TLegend* leg22;
  if (toUnfold == "pt" && doNormalized) leg22 = new TLegend(0.55,0.65,0.8,0.8);
  else leg22 = new TLegend(0.55,0.65,0.8,0.88);
  for (int is = 4; is < nSYST; is++){
    leg22->AddEntry(h_syst[is],longnames[is],"lp");
  }
  leg22->SetFillStyle(0);
  leg22->SetBorderSize(0);
  leg22->SetTextSize(0.04);
  leg22->SetTextFont(42);

  h_dummy_r->Draw("hist");
  h_syst_tot->Draw("hist,same");
  h_syst_stat->Draw("ep,same");
  if (doNormalized == false) h_lumi->Draw("ep,same");
  for (int is = 0; is < nSYST; is++){
    h_syst[is]->Draw("ep,same");
  }
  h_dummy_r->Draw("hist,axis,same");
  leg2->Draw(); 
  leg22->Draw(); 

  drawCMS(0.18,0.94,0.28,0.94,false,true,false);

  c1->SaveAs("UnfoldingPlots/unfold_relative_uncertainties_"+toUnfold+"_"+channel+normflag+".pdf");

  TFile fout = TFile("dataForRivet_"+toUnfold+"_"+channel+normflag+".root","recreate");
  if (toUnfold == "y") symmetrize_err(h_fullunc);
  h_fullunc->Write();

  fout.Close();

  // Normalize diagonal
  float diags[10] = {0.};
  for (int id = 0; id < nBINS; id++){
    int idbin = h_cov_tot->GetBin(id+1,id+1);
    diags[id] = sqrt(h_cov_tot->GetBinContent(idbin));
  }
  for (int ix = 1; ix < nBINS-1; ix++){
    for (int iy = 1; iy < nBINS-1; iy++){
      int ibin = h_cov_tot->GetBin(ix+1,iy+1);
      float rawval = h_cov_tot->GetBinContent(ibin);
      rawval *= 100.0 / diags[ix] / diags[iy];
      int roundval = rawval * 1000;
      float finalval = (float) roundval / 1000.;
      h_cov_tot->SetBinContent(ibin,finalval);
    }
  }

  TString binnames[5] = {"Bin 1","Bin 2","Bin 3","Bin 4","Bin 5"};
  for (int ii = 0; ii < nBINS - 2; ii++){
    h_cov_tot->GetXaxis()->SetBinLabel(ii+2,binnames[ii]);
    h_cov_tot->GetYaxis()->SetBinLabel(ii+2,binnames[ii]);
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",900,700);
  c2->SetLeftMargin(0.18);
  c2->SetRightMargin(0.12);
  c2->SetBottomMargin(0.1);
  h_cov_tot->SetAxisRange(1.0,nBINS-1.001,"X");
  h_cov_tot->SetAxisRange(1.0,nBINS-1.001,"Y");
  h_cov_tot->GetZaxis()->SetLabelSize(0.04);
  h_cov_tot->Draw("COLZ");
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)h_cov_tot->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.88);
  palette->SetX2NDC(0.93);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  gPad->Modified();
  gPad->Update();
  c2->SaveAs("UnfoldingPlots/total_covariance_"+toUnfold+"_"+channel+".pdf");

  return;

}

void symmetrize(TH1F* h_input){
  if (h_input->GetNbinsX() != 6) {
    cout << "Error! Wrong histogram being symmetrized" << endl;
    continue;
  }
  for (int ii = 0; ii < 3; ii++){
    float temp1 = h_input->GetBinContent(ii+1);
    float temp2 = h_input->GetBinContent(6-ii);
    float average = (temp1 + temp2) / 2;
    h_input->SetBinContent(ii+1,average);
    h_input->SetBinContent(6-ii,average);
  }

  for (int ii = 0; ii < 3; ii++){
    float temp1 = h_input->GetBinError(ii+1);
    float temp2 = h_input->GetBinError(6-ii);
    float average = (temp1 + temp2) / 2;
    h_input->SetBinError(ii+1,average);
    h_input->SetBinError(6-ii,average);
  }
  
  return;
}

void symmetrize_err(TH1F* h_input){
  if (h_input->GetNbinsX() != 6) {
    cout << "Error! Wrong histogram being symmetrized" << endl;
    continue;
  }
  for (int ii = 0; ii < 3; ii++){
    float temp1 = h_input->GetBinError(ii+1);
    float temp2 = h_input->GetBinError(6-ii);
    float average = (temp1 + temp2) / 2;
    h_input->SetBinError(ii+1,average);
    h_input->SetBinError(6-ii,average);
  }
  return;
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
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.18);
  
  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.2);
  
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
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
}


void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize, char *text) {
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void drawCMS(Double_t x1, Double_t y1, Double_t x2, Double_t y2, bool pad, bool prel, bool forPublic) {

  float cmsTextSize = 0.08;
  if (!pad) cmsTextSize = 0.058;
  float extraOverCmsTextSize = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;


  TLatex l;
  l.SetTextSize(cmsTextSize); 
  l.SetTextFont(61); 
  l.SetTextAngle(0);
  l.SetNDC();
  l.SetTextColor(1);
  if (forPublic) l.DrawLatex(x1,y1,"CMS");

  if (prel && forPublic) {
    TLatex lp;
    lp.SetTextSize(extraTextSize); 
    lp.SetTextFont(52); 
    lp.SetNDC();
    lp.SetTextColor(1);
    lp.DrawLatex(x2,y2,"Preliminary");
  }

  TLatex ll;
  ll.SetTextSize(extraTextSize); 
  ll.SetTextFont(42); 
  ll.SetTextAngle(0);
  ll.SetNDC();
  ll.SetTextColor(1);
  if (pad) ll.DrawLatex(0.7,0.93,"12.4 fb^{-1} (13 TeV)");
  else ll.DrawLatex(0.72,0.94,"12.4 fb^{-1} (13 TeV)");

}
