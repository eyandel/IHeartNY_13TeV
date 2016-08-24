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

#include "RooUnfold/src/RooUnfold.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldDagostini.h"
#include "RooUnfold/src/RooUnfoldErrors.h"
#include "RooUnfold/src/RooUnfoldInvert.h"
#include "RooUnfold/src/RooUnfoldParms.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
#include "RooUnfold/src/RooUnfoldTUnfold.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,Double_t tsize,char *text); 


void unfoldPlots(TString which="", bool y=false, bool doElectron=false) {
  
  gSystem->Load("RooUnfold/libRooUnfold");

  SetPlotStyle();

  TString muOrEl = "mu";
  if (doElectron) muOrEl = "el";

  TString what = "_pt";
  if (y) what = "_y";

  TFile* f = new TFile("UnfoldingPlots/closureTest"+which+"_full.root");

  TH2F* h_response = (TH2F*) f->Get("responseMatrix"+which);

  TH1F* h_stability = (TH1F*) h_response->ProjectionX()->Clone();
  h_stability->Reset();
  h_stability->SetName("stability");

  TH1F* h_purity = (TH1F*) h_response->ProjectionX()->Clone();
  h_purity->Reset();
  h_purity->SetName("purity");

  TH1F* h_responseX = (TH1F*) h_response->ProjectionX()->Clone();
  h_responseX->SetName("responseX");

  TH1F* h_responseY = (TH1F*) h_response->ProjectionY()->Clone();
  h_responseY->SetName("responseY");


  TH1F* h_reco = (TH1F*) f->Get("ptRecoTop"+which+"_measured");
  TH1F* h_gen  = (TH1F*) f->Get("ptGenTop"+which+"_true");
  TH1F* h_eff = (TH1F*) h_reco->Clone();
  h_eff->SetName("efficiency");
  h_eff->Divide(h_gen);
    

  // stability
  for (int ix=0; ix<(int)h_responseX->GetNbinsX(); ix++) {

    float nrec = h_responseX->GetBinContent(ix+1);
    float erec = h_responseX->GetBinError(ix+1);

    if (nrec==0) continue;
    
    for (int iy=0; iy<(int)h_responseY->GetNbinsX(); iy++) {

      if (ix==iy) {
	float nrecgen = h_response->GetBinContent(ix+1,iy+1);
	float erecgen = h_response->GetBinError(ix+1,iy+1);
	
	float eA = erecgen/nrec;
	float eB = (nrecgen/(nrec*nrec))*erec;
	float err = sqrt(eA*eA + eB*eB);      
	
	h_stability->SetBinContent(ix+1, nrecgen/nrec);
	h_stability->SetBinError(ix+1, err);
      }
      
    }//end iy loop
  }//end ix loop

  // purity
  for (int iy=0; iy<(int)h_responseY->GetNbinsX(); iy++) {
    
    float ngen = h_responseY->GetBinContent(iy+1);
    float egen = h_responseY->GetBinError(iy+1);
    
    if (ngen==0) continue;
    
    for (int ix=0; ix<(int)h_responseX->GetNbinsX(); ix++) {

      if (ix==iy) {
	float nrecgen = h_response->GetBinContent(ix+1,iy+1);
	float erecgen = h_response->GetBinError(ix+1,iy+1);

	float eA = erecgen/ngen;
	float eB = (nrecgen/(ngen*ngen))*egen;
	float err = sqrt(eA*eA + eB*eB);      
	
	h_purity->SetBinContent(ix+1, nrecgen/ngen);
	h_purity->SetBinError(ix+1, err);
      }
      
    }//end iy loop
  }//end ix loop


  // ----------------------------------------------------------------------------------------------------------------
  // plot etc...
  // ----------------------------------------------------------------------------------------------------------------

  h_purity->SetAxisRange(0,1.2,"Y");
  if (!y) h_purity->SetAxisRange(400.,1150.0,"X");
  h_purity->GetYaxis()->SetTitleOffset(1.1);
  h_purity->GetYaxis()->SetTitle("");  
  h_purity->GetYaxis()->SetTitle("Fractional");
  h_purity->GetXaxis()->SetTitle("Top quark p_{T} [GeV]");


  TCanvas c;
  h_purity->SetLineColor(4); 
  h_purity->SetLineWidth(2);
  h_purity->Draw("hist");
  h_stability->SetLineColor(2); 
  h_stability->SetLineWidth(2);
  h_stability->Draw("same,hist");
  h_eff->SetLineColor(8); 
  h_eff->SetLineWidth(2);
  h_eff->Draw("same,hist");

  TLegend* leg = new TLegend(0.22,0.76,0.47,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_purity, " Purity", "l");
  leg->AddEntry(h_stability, " Stability", "l");
  leg->AddEntry(h_eff, " Efficiency", "l");
  leg->Draw();

  c.SaveAs("purity-stability"+which+"_"+muOrEl+what+".png");
  c.SaveAs("purity-stability"+which+"_"+muOrEl+what+".eps");


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
  /*
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.16);
    gStyle->SetPadLeftMargin(0.16);
  */
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
