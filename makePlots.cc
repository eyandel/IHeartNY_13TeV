#include "makePlots.h"

#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "THStack.h"
#include "TExec.h"
#include "TPaletteAxis.h"

#include <iostream>
#include <iomanip>


void setStyle() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(000000);
  
  gStyle->SetTitleFont(43);
  gStyle->SetTitleFont(43, "XYZ");
  gStyle->SetTitleSize(28, "XYZ");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.0, "Y");
  gStyle->SetLabelFont(43, "XYZ");
  gStyle->SetLabelSize(24, "XYZ");

  //gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(0.);

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
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  
}

void myLargeText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.06); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void myText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.043); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void myItalicText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.05); 
  l.SetTextFont(52); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
void mySmallText(Double_t x,Double_t y,Color_t color,char const *text) {
  TLatex l;
  l.SetTextSize(0.042); 
  l.SetTextFont(42); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void drawCMS(Double_t x,Double_t y, bool prel) {

  float cmsTextSize = 0.07;
  float extraOverCmsTextSize = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  TLatex l;
  l.SetTextSize(cmsTextSize); 
  l.SetTextFont(61); 
  l.SetTextAngle(0);
  l.SetNDC();
  l.SetTextColor(1);
  l.DrawLatex(x,y,"CMS");

  if (prel) {
    TLatex lp;
    lp.SetTextSize(extraTextSize); 
    lp.SetTextFont(52); 
    lp.SetNDC();
    lp.SetTextColor(1);
    float offset = 0.09;
    lp.DrawLatex(x+offset,y,"Preliminary");
  }

  TLatex ll;
  ll.SetTextSize(extraTextSize); 
  ll.SetTextFont(42); 
  ll.SetTextAngle(0);
  ll.SetNDC();
  ll.SetTextColor(1);
  ll.DrawLatex(0.69,0.94,"19.7 fb^{-1} (8 TeV)");

}


// -------------------------------------------------------------------------------------
// make pretty plots
// -------------------------------------------------------------------------------------

void makePlots(TString DIR, TString DIRqcd, TString channel, TString var, TString region, bool useQCDMC = false, bool unBlind = false, bool usePost = false) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+region;
  
  // get histograms
  SummedHist* wjets = getWJets( DIR, var, region, channel, false, "nom", usePost );
  SummedHist* singletop = getSingleTop( DIR, var, region, channel, false, "nom", usePost );
  SummedHist* ttbar = getTTbar( DIR, var, region, channel, false, "nom", usePost );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, region, channel, false, "nom", usePost );
  SummedHist* data = getData( DIR, var, region, channel, false);

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  
  TH1F* h_qcd;
  if (useQCDMC) {
    SummedHist* qcd = getQCDMC(DIR, var, region, channel, false, "nom", usePost);
    h_qcd = (TH1F*) qcd->hist();
  }
  else {
    h_qcd = getQCDData( DIR, DIRqcd, var, region, channel, "nom", usePost); // Currently taking QCD from data sideband with normalization from MC in signal region
  }
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  TH1F* h_singletop = (TH1F*) singletop->hist();
  TH1F* h_data = (TH1F*) data->hist();

  // -----------------------------------------------------------------------------------------------------
  // Apply post-fit normalizations if desired
  // The scale factors are determined from the post-fit nuisance parameters, not the post-fit event yields

  if (usePost){
    h_qcd->Scale(0.89);
    h_singletop->Scale(0.95);
    h_wjets->Scale(1.06);
    h_ttbar->Scale(0.79);
    h_ttbar_nonSemiLep->Scale(0.79);
  }

  // -------------------------------------------------------------------------------------
  // various hist plotting edits
  if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))){
    if (h_qcd) h_qcd->Rebin(5);
    if (h_wjets) h_wjets->Rebin(5);
    if (h_singletop) h_singletop->Rebin(5);
    if (h_ttbar_nonSemiLep) h_ttbar_nonSemiLep->Rebin(5);
    if (h_ttbar) h_ttbar->Rebin(5);
    h_data->Rebin(5);
  }

  h_data->UseCurrentStyle();
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(1);

  // create stack & summed histogram for ratio plot
  THStack* h_stack = new THStack();    
  if (h_qcd) h_stack->Add(h_qcd);
  if (h_wjets) h_stack->Add(h_wjets);
  if (h_singletop) h_stack->Add(h_singletop);
  if (h_ttbar_nonSemiLep) h_stack->Add(h_ttbar_nonSemiLep);
  if (h_ttbar) h_stack->Add(h_ttbar);

  TH1F* h_totalbkg = (TH1F*) h_ttbar->Clone("totalbkg_"+hist);
  if (h_ttbar_nonSemiLep) h_totalbkg->Add(h_ttbar_nonSemiLep);
  if (h_wjets) h_totalbkg->Add(h_wjets);
  if (h_singletop) h_totalbkg->Add(h_singletop);
  if (h_qcd) h_totalbkg->Add(h_qcd);

  // -------------------------------------------------------------------------------------
  // poisson errors
  h_data->SetBinErrorOption(TH1::kPoisson);
  for (int ibindata=0; ibindata<(int)h_data->GetNbinsX(); ibindata++) {
    double err_low = h_data->GetBinErrorLow(ibindata);
    double err_up = h_data->GetBinErrorUp(ibindata);
  }
  // -------------------------------------------------------------------------------------
  
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {

    float stat_tt, stat_tt_non, stat_st, stat_wj, stat_qcd = 0.0;
    stat_tt = h_ttbar->GetBinError(ib+1);
    if (h_ttbar_nonSemiLep) stat_tt_non = h_ttbar_nonSemiLep->GetBinError(ib+1);
    if (h_singletop) stat_st = h_singletop->GetBinError(ib+1);
    if (h_wjets) stat_wj = h_wjets->GetBinError(ib+1);
    if (h_qcd) stat_qcd = h_qcd->GetBinError(ib+1);

    float binbkg = sqrt( stat_tt*stat_tt + // Stat error only for now
			 stat_tt_non*stat_tt_non + 
			 stat_st*stat_st + 
			 stat_wj*stat_wj + 
			 stat_qcd*stat_qcd 
			 );
    
    h_totalbkg->SetBinError(ib+1, binbkg);
  }

  if (hist.Contains("ak8jetPt")) h_data->GetXaxis()->SetRangeUser(400.,1200.);
  
  TH1F* h_ratio;
  TH1F* h_ratio2;
  h_ratio = (TH1F*) h_data->Clone("ratio_"+hist);  // Data / MC
  h_ratio->Sumw2();
  h_ratio->Divide(h_totalbkg);
  
  h_ratio2 = (TH1F*) h_totalbkg->Clone("ratio2_"+hist); // Uncertainty on Data / MC
  h_ratio2->Sumw2();
  for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {
    h_ratio2->SetBinContent(ib+1, 1.0);
    float tmperr = h_totalbkg->GetBinError(ib+1);
    float tmpcount = h_totalbkg->GetBinContent(ib+1);
    
    if (tmpcount==0) continue;
    if (ib==0) cout << "tmperr = " << tmperr << " tmpcount = " << tmpcount << endl;
    h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
  }
  
  float mymax = max(h_data->GetMaximum(),h_totalbkg->GetMaximum());
  h_data->SetAxisRange(0,mymax*1.05,"Y");

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1;
  TPad* p2;
  if (region == "Pre" || unBlind){
    p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
    p1->SetTopMargin(0.08);
    p1->SetBottomMargin(0.05);
    p1->SetNumber(1);
    p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
    p2->SetNumber(2);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.35);
    
    p1->Draw();
    p2->Draw();
    p1->cd();
    
    h_ratio2->SetMarkerSize(0);
    h_ratio2->SetLineColor(0);
    h_ratio2->SetFillColor(15);
    h_ratio2->SetFillStyle(1001);
    
    h_ratio->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    
    //h_data->GetXaxis()->SetLabelSize(26);
    h_data->GetYaxis()->SetLabelSize(26);
    h_data->GetYaxis()->SetTitleSize(32);
    h_data->GetYaxis()->SetTitleOffset(1.4);
    h_data->GetXaxis()->SetTitle("");
    
    h_data->Draw("LE0P");
  }

  h_totalbkg->UseCurrentStyle();
  h_totalbkg->SetFillColor(0);
  h_totalbkg->SetLineWidth(1);
  h_totalbkg->SetLineColor(1);
  
  if (region == "Pre" || unBlind) h_totalbkg->Draw("hist,same");
  else h_totalbkg->Draw("hist");
  
  h_stack->Draw("hist,same");
  if (region == "Pre" || unBlind) h_data->Draw("LE0P,same");

  float xmin = 0.71;
  float ymin = 0.52;

  float xwidth = 0.20;
  float ywidth = 0.38;

  if (region != "Pre" && !unBlind){
    ymin = 0.57;
    xwidth = 0.18;
    ywidth = 0.32;
  }

  //Legend top left
  if ((var == "ak8jetTau32" || var == "ak8jetTau21") || (var == "ak4jetCSV" && (region == "0t1b" || region == "1t0b" || region == "1t1b")) || ((var == "ak8jetMass" || var == "ak8jetSDm01" || var == "ak8jetSDmass" || var == "ak8jetSDsubjetMaxCSV") && (region == "1t0b" || region == "1t1b"))) xmin = 0.16;

  //Legend bottom center
  if (var.Contains("Phi")) {
    xmin = 0.40;
    ymin = 0.12;
    if (region != "Pre" && !unBlind) ymin = 0.25;
    if ((region == "0t1b" || region == "1t0b" || region == "1t1b" ) && var == "ak4jetPhi") ymin = 0.45;
  }

  //Legend center center
  if (var.Contains("Phi") && channel == "el" && (region == "1t0b" || region == "1t1b")){
    ymin = 0.35;
    if (var == "ak4jetPhi") ymin = 0.31;
  }
  
  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  if (region == "Pre" || unBlind) leg->SetTextSize(0.045);
  else leg->SetTextSize(0.037);
  if (region == "Pre" || unBlind) leg->AddEntry(h_data, "Data", "pe");
  leg->AddEntry(h_ttbar, "t#bar{t} signal", "f");
  leg->AddEntry(h_ttbar_nonSemiLep, "t#bar{t} other", "f");
  leg->AddEntry(h_singletop, "Single t", "f");
  leg->AddEntry(h_wjets, "W+jets", "f");
  leg->AddEntry(h_qcd, "Multijet" , "f");
  leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
  leg->Draw();

  myText(0.10,0.94,1,"#intLdt = 2.7 fb^{-1}");
  myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");

  // plot ratio part
  if (region == "Pre" || unBlind){
    p2->cd();
    p2->SetGridy();
    h_ratio->UseCurrentStyle();
    h_ratio->SetMarkerStyle(8);
    h_ratio->SetMarkerSize(1);
    h_ratio->Draw("le0p");
    h_ratio2->Draw("same,e2");
    h_ratio->Draw("le0p,same");
    h_ratio->SetMaximum(1.8);
    h_ratio->SetMinimum(0.2);
    h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetXaxis()->SetLabelSize(26);
    h_ratio->GetYaxis()->SetLabelSize(26);
    h_ratio->GetXaxis()->SetTitleOffset(2.8);
    h_ratio->GetYaxis()->SetTitleOffset(1.4);
    h_ratio->GetXaxis()->SetTitleSize(32);
    h_ratio->GetYaxis()->SetTitleSize(32);
  }

  // save output
  TString outname = "Plots/"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_stack->Delete();
  leg->Delete();
}

void compareShapes(TString DIR, TString DIRqcd, TString channel, TString var, TString region) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+region;
  
  // get histograms
  SummedHist* wjets = getWJets( DIR, var, region, channel, false, "nom" );
  SummedHist* ttbar = getTTbar( DIR, var, region, channel, false, "nom" );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, var, region, channel, false, "nom" );
  //SummedHist* singletop = getSingleTop( DIR, var, region, channel, false, "nom" );

  // -------------------------------------------------------------------------------------
  // get the TH1F versions
  
  TH1F* h_qcd = getQCDData( DIR, DIRqcd, var, region, channel, "nom"); // Currently taking QCD from data sideband with no triangular cut,
                                                                // with normalization from MC in signal region
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  h_ttbar->Add(h_ttbar_nonSemiLep);
  //TH1F* h_singletop = (TH1F*) singletop->hist();

  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  
  if (h_qcd){
    h_qcd->SetLineColor(h_qcd->GetFillColor());
    h_qcd->SetFillColor(0);
    h_qcd->Scale(1./h_qcd->Integral());
  }
  h_wjets->SetLineColor(h_wjets->GetFillColor());
  h_wjets->SetFillColor(0);
  h_wjets->Scale(1./h_wjets->Integral());
  h_ttbar->SetLineColor(h_ttbar->GetFillColor());
  h_ttbar->SetFillColor(0);
  h_ttbar->Scale(1./h_ttbar->Integral());
  //h_singletop->SetLineColor(h_singletop->GetFillColor());
  //h_singletop->SetFillColor(0);
  //h_singletop->Scale(1./h_singletop->Integral());

  if (!(hist.Contains("nAK4jet") || hist.Contains("nAK8jet") || hist.Contains("nBjet") || hist.Contains("nTjet"))){
    if (h_qcd) h_qcd->Rebin(5);
    h_wjets->Rebin(5);
    h_ttbar->Rebin(5);
    //h_singletop->Rebin(5);
  }

  h_ttbar->SetMaximum(h_ttbar->GetMaximum()*1.5);
  
  h_ttbar->Draw("hist");
  h_wjets->Draw("hist,same");
  if (h_qcd) h_qcd->Draw("hist,same");
  //h_singletop->Draw("hist,same");

  float xmin = 0.70;
  float ymin = 0.70;

  if (hist.Contains("Tau32") || hist.Contains("SDmass")) xmin = 0.2;

  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmin+0.2,ymin+0.2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.042);
  leg->AddEntry(h_ttbar, "t#bar{t}", "l");
  leg->AddEntry(h_wjets, "W+jets", "l");
  if (h_qcd) leg->AddEntry(h_qcd, "Multijet" , "l");
  //leg->AddEntry(h_singletop, "Single t", "l");
  leg->Draw();

  // save output
  TString outname = "Plots/compShapes_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  leg->Delete();
  if (h_qcd) h_qcd->Delete();
  h_ttbar->Delete();
  h_ttbar_nonSemiLep->Delete();
  h_wjets->Delete();
  //h_singletop->Delete();
  delete wjets;
  delete ttbar;
  delete ttbar_nonSemiLep;
  //delete singletop;
}


void makeCombineInputs(TString DIR, TString DIRqcd) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  const int nchannels = 2;
  TString channels[nchannels] = {"mu","el"};
  const int nhist = 5;
  TString histnames[nhist] = {"lepAbsEta","ak8jetTau21","lepAbsEta","ak8jetTau32","ak8jetSDmass"};
  TString regions[nhist] = {"0t","0t","1t0b","1t0b","1t1b"};
  int rebinby[nhist] = {5,5,5,5,5};
  const int nsys = 13;
  TString sysnames[nsys] = {"nom","puUp","puDown","JECUp","JECDown","JERUp","JERDown","lepUp","lepDown","BTagUp","BTagDown","TopTagUp","TopTagDown"};
  
  TH1F* h_qcd[nchannels][nhist][nsys];
  TH1F* h_wjets[nchannels][nhist][nsys];
  TH1F* h_ttbar[nchannels][nhist][nsys];
  TH1F* h_singletop[nchannels][nhist][nsys];
  TH1F* h_data[nchannels][nhist];
  
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){
      for (int kk = 0; kk < nsys; kk++){
	TString append = "";
	if (sysnames[kk] != "nom") append = "_" + sysnames[kk];
      
	// get histograms
	SummedHist* wjets = getWJets( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk] );
	SummedHist* singletop = getSingleTop( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk] );
	SummedHist* ttbar = getTTbar( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk] );
	SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( DIR, histnames[jj], regions[jj], channels[ii], false, sysnames[kk] );
	
	// -------------------------------------------------------------------------------------
	// get the TH1F versions
	
	h_qcd[ii][jj][kk] = (TH1F*) getQCDData( DIR, DIRqcd, histnames[jj], regions[jj], channels[ii], sysnames[kk])->Clone("QCD"+append);
	h_wjets[ii][jj][kk] = (TH1F*) wjets->hist()->Clone("WJets"+append);
	h_ttbar[ii][jj][kk] = (TH1F*) ttbar->hist()->Clone("TTbar"+append);
	h_ttbar[ii][jj][kk]->Add((TH1F*) ttbar_nonSemiLep->hist());
	h_singletop[ii][jj][kk] = (TH1F*) singletop->hist()->Clone("SingleTop"+append);

	if (sysnames[kk] == "TopTagUp"){
	  h_qcd[ii][jj][kk]->SetName("QCD_TopMisTagUp");
	  h_wjets[ii][jj][kk]->SetName("WJets_TopMisTagUp");
	  h_singletop[ii][jj][kk]->SetName("SingleTop_TopMisTagUp");
	}
	
	if (sysnames[kk] == "TopTagDown"){
	  h_qcd[ii][jj][kk]->SetName("QCD_TopMisTagDown");
	  h_wjets[ii][jj][kk]->SetName("WJets_TopMisTagDown");
	  h_singletop[ii][jj][kk]->SetName("SingleTop_TopMisTagDown");
	}
	
	h_qcd[ii][jj][kk]->Rebin(rebinby[jj]);
	h_wjets[ii][jj][kk]->Rebin(rebinby[jj]);
	h_ttbar[ii][jj][kk]->Rebin(rebinby[jj]);
	h_singletop[ii][jj][kk]->Rebin(rebinby[jj]);
      }

      SummedHist* data = getData( DIR, histnames[jj], regions[jj], channels[ii], false);
      h_data[ii][jj] = (TH1F*) data->hist()->Clone("data_obs");
      h_data[ii][jj]->Rebin(rebinby[jj]);
    }
  }

  // Write input file
  TFile* fout = new TFile("combineInputs.root","recreate");
  TDirectory* topdir = fout->GetDirectory("");
  TDirectory* dirs[nchannels][nhist];
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){
      dirs[ii][jj] = fout->mkdir(histnames[jj]+regions[jj]+"_"+channels[ii]);
      dirs[ii][jj]->cd();
      for (int kk = 0; kk < nsys; kk++){
	h_qcd[ii][jj][kk]->Write();
	h_wjets[ii][jj][kk]->Write();
	h_singletop[ii][jj][kk]->Write();
	h_ttbar[ii][jj][kk]->Write();
      }
      h_data[ii][jj]->Write();
      topdir->cd();
    }
  }
  fout->Close();

  // Make comparison plots
  int colors[nsys-1] = {2,2,3,3,4,4,5,5,6,6,7,7};
  int styles[nsys-1] = {1,2,1,2,1,2,1,2,1,2,1,2};
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj++){

      TCanvas* c = new TCanvas("c_"+histnames[jj]+regions[jj]+"_"+channels[ii],"c_"+histnames[jj]+regions[jj]+"_"+channels[ii],900,800);

      h_ttbar[ii][jj][0]->SetFillColor(0);
      h_ttbar[ii][jj][0]->SetMaximum(1.5*h_ttbar[ii][jj][0]->GetMaximum());
      h_ttbar[ii][jj][0]->Draw("hist");
    
      TLegend* leg;
      float xmin = 0.66;
      if (regions[jj] == "1t1b" || regions[jj] == "1t0b") xmin = 0.2;
      leg = new TLegend(xmin,0.65,xmin+0.2,0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      leg->AddEntry(h_ttbar[ii][jj][0], "TTbar nominal", "l");

      for (int kk = 1; kk < nsys; kk++){
	h_ttbar[ii][jj][kk]->SetFillColor(0);
	h_ttbar[ii][jj][kk]->SetLineColor(colors[kk-1]);
	h_ttbar[ii][jj][kk]->SetLineStyle(styles[kk-1]);
	h_ttbar[ii][jj][kk]->Draw("hist,same");
	if (kk%2 == 1) leg->AddEntry(h_ttbar[ii][jj][kk],sysnames[kk],"l");
      }

      leg->Draw();

      TString outname = "Plots/sysVar_"+channels[ii]+"_"+histnames[jj]+regions[jj]+".pdf";
      c->SaveAs(outname);
    }
  }
  
  // Cleanup
  for (int ii = 0; ii < nchannels; ii++){
    for (int jj = 0; jj < nhist; jj ++){
      for (int kk = 0; kk < nsys; kk++){
	h_qcd[ii][jj][kk]->Delete();
	h_wjets[ii][jj][kk]->Delete();
	h_singletop[ii][jj][kk]->Delete();
	h_ttbar[ii][jj][kk]->Delete();
      }
      h_data[ii][jj]->Delete();
    }
  }
}

void makeQCDComp(TString DIR, TString DIRqcd, TString channel, TString var) {
  
  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hist = var+"Pre";
  
  // get histograms
  SummedHist* sigQCDMC = getQCDMC(DIR, var, "Pre", channel, false, "nom");
  SummedHist* sideQCDMC = getQCDMC(DIRqcd, var, "Pre", channel, true, "nom");
  TH1F* h_sigQCDMC = (TH1F*) sigQCDMC->hist();
  TH1F* h_sideQCDMC = (TH1F*) sideQCDMC->hist();

  TH1F* h_sideQCDData = getQCDData(DIR,DIRqcd,var,"Pre",channel,"nom");

  h_sigQCDMC->SetLineColor(2);
  h_sigQCDMC->SetMarkerColor(2);
  h_sideQCDMC->SetLineColor(4);
  h_sideQCDMC->SetMarkerColor(4);
  h_sideQCDData->SetLineColor(1);
  h_sideQCDData->SetMarkerColor(1);
  h_sigQCDMC->Sumw2();
  h_sideQCDMC->Sumw2();
  h_sideQCDData->Sumw2();

  if (var != "nAK4jet" && var != "nAK8jet" && var != "nBjet" && var != "nTjet"){
    h_sigQCDMC->Rebin(5);
    h_sideQCDMC->Rebin(5);
    h_sideQCDData->Rebin(5);
  }
  
  h_sigQCDMC->Scale(1./h_sigQCDMC->Integral());
  h_sideQCDMC->Scale(1./h_sideQCDMC->Integral());
  h_sideQCDData->Scale(1./h_sideQCDData->Integral());

  TH1F* h_ratio = (TH1F*) h_sideQCDData->Clone();
  TH1F* h_ratio2 = (TH1F*) h_sigQCDMC->Clone();
  TH1F* h_ratio3 = (TH1F*) h_sideQCDMC->Clone();
  h_ratio->Divide(h_sideQCDData);
  h_ratio2->Divide(h_sideQCDData);
  h_ratio3->Divide(h_sideQCDData);
  
  // -------------------------------------------------------------------------------------
  // plotting!

  TCanvas* c = new TCanvas("c_"+hist,"c_"+hist,900,800);
  TPad* p1 = new TPad("datamcp1_"+hist,"datamcp1_"+hist,0,0.3,1,1);
  p1->SetTopMargin(0.08);
  p1->SetBottomMargin(0.05);
  p1->SetNumber(1);
  TPad* p2 = new TPad("datamcp2_"+hist,"datamcp2_"+hist,0,0,1,0.3);
  p2->SetNumber(2);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.35);

  p1->Draw();
  p2->Draw();
  p1->cd();

  h_ratio->GetXaxis()->SetTitle(h_sigQCDMC->GetXaxis()->GetTitle());

  h_sideQCDData->GetXaxis()->SetTitle("");
    
  h_sideQCDData->Draw();
  h_sideQCDMC->Draw("same");
  h_sigQCDMC->Draw("same");

  float xmin = 0.66;
  float xmax = 0.85;
  float ymin = 0.65;
  float ymax = 0.88;
  
  // legend
  TLegend* leg;
  leg = new TLegend(xmin,ymin,xmax,ymax);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->AddEntry(h_sideQCDData, "Data-driven QCD", "l");
  leg->AddEntry(h_sideQCDMC, "QCD MC, sideband", "l");
  leg->AddEntry(h_sigQCDMC, "QCD MC, sig. region", "l");
  leg->Draw();

  // plot ratio part
  p2->cd();
  //h_ratio->UseCurrentStyle();
  h_ratio->Draw();
  h_ratio2->Draw("same");
  h_ratio3->Draw("same");
  h_ratio->SetMaximum(2.0);
  h_ratio->SetMinimum(0.0);
  h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
  h_ratio->GetYaxis()->SetTitle("MC / Data-driven");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleOffset(1.0);
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleSize(0.1);

  // save output
  TString outname = "Plots/compQCD_"+channel+"_"+hist+".pdf";
  c->SaveAs(outname);

  // Cleanup, since we are using gDirectory(kFALSE)
  h_sideQCDData->Delete();
  h_sideQCDMC->Delete();
  h_sigQCDMC->Delete();
  h_ratio->Delete();
  h_ratio2->Delete();
  h_ratio3->Delete();
}

// -------------------------------------------------------------------------------------
// print cutflow latex table
// -------------------------------------------------------------------------------------
void makeTable(TString DIR, TString DIRqcd, TString channel, bool inSideband, bool useQCDMC) {

  const int nhist = 4;
  TString what = "lepPhi";
  TString regions[nhist] = {"Pre","0t","1t0b","1t1b"};

  // counts for cutflow
  float count_tt_semiLep[nhist] = {0};
  float count_tt_nonsemilep[nhist] = {0};
  float count_singletop[nhist] = {0};
  float count_wjets[nhist] = {0};
  float count_qcd[nhist] = {0};
  float count_tot[nhist] = {0};
  float count_data[nhist] = {0};
  
  // errors for cutflow
  float err_tt_semiLep[nhist] = {0};
  float err_tt_nonsemilep[nhist] = {0};
  float err_singletop[nhist] = {0};
  float err_wjets[nhist] = {0};
  float err_qcd[nhist] = {0};
  float err_tot[nhist] = {0};

  // Get count and error for each sample
  for (int ii = 0; ii < nhist; ii ++){
    // get histograms
    SummedHist* wjets = getWJets(DIR, what, regions[ii], channel, inSideband, "nom" );
    SummedHist* singletop = getSingleTop(DIR, what, regions[ii], channel, inSideband, "nom" );
    SummedHist* ttbar = getTTbar(DIR, what, regions[ii], channel, inSideband, "nom" );
    SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep(DIR, what, regions[ii], channel, inSideband, "nom" );
    SummedHist* data = getData(DIR, what, regions[ii], channel, inSideband);

    TH1F* h_wjets = (TH1F*) wjets->hist();
    TH1F* h_ttbar_semiLep = (TH1F*) ttbar->hist();
    TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
    TH1F* h_singletop = (TH1F*) singletop->hist();
    TH1F* h_data = (TH1F*) data->hist();

    TH1F* h_qcd;
    if(inSideband || useQCDMC) {
      SummedHist* qcd = getQCDMC(DIRqcd, what, regions[ii], channel, inSideband, "nom");
      h_qcd = (TH1F*) qcd->hist();
    }
    else h_qcd = getQCDData(DIR,DIRqcd,what, regions[ii],channel,"nom");
    
    // error on pre-fit counts
    for (int ib=0; ib<h_ttbar_semiLep->GetNbinsX(); ib++) {
      if (h_ttbar_semiLep) err_tt_semiLep[ii]    += h_ttbar_semiLep->GetBinError(ib+1)*h_ttbar_semiLep->GetBinError(ib+1);
      if (h_ttbar_nonSemiLep) err_tt_nonsemilep[ii] += h_ttbar_nonSemiLep->GetBinError(ib+1)*h_ttbar_nonSemiLep->GetBinError(ib+1);
      if (h_singletop) err_singletop[ii]     += h_singletop->GetBinError(ib+1)*h_singletop->GetBinError(ib+1);
      if (h_wjets) err_wjets[ii]         += h_wjets->GetBinError(ib+1)*h_wjets->GetBinError(ib+1);
      if (h_qcd) err_qcd[ii]           += h_qcd->GetBinError(ib+1)*h_qcd->GetBinError(ib+1);
    }

    err_tt_semiLep[ii]    = sqrt(err_tt_semiLep[ii]);
    err_tt_nonsemilep[ii] = sqrt(err_tt_nonsemilep[ii]);
    err_singletop[ii]     = sqrt(err_singletop[ii]);
    err_wjets[ii]         = sqrt(err_wjets[ii]);
    err_qcd[ii]          = sqrt(err_qcd[ii]);
    
    err_tot[ii] = err_tt_semiLep[ii]*err_tt_semiLep[ii] + err_tt_nonsemilep[ii]*err_tt_nonsemilep[ii] + err_singletop[ii]*err_singletop[ii] + err_wjets[ii]*err_wjets[ii] + err_qcd[ii]*err_qcd[ii];
    err_tot[ii] = sqrt(err_tot[ii]);

    if (h_ttbar_semiLep) count_tt_semiLep[ii] = h_ttbar_semiLep->GetSum();
    if (h_ttbar_nonSemiLep) count_tt_nonsemilep[ii] = h_ttbar_nonSemiLep->GetSum();
    if (h_singletop) count_singletop[ii] = h_singletop->GetSum();
    if (h_wjets) count_wjets[ii] = h_wjets->GetSum();
    if (h_qcd) count_qcd[ii] = h_qcd->GetSum();
    count_tot[ii] = count_tt_semiLep[ii] + count_tt_nonsemilep[ii] + count_singletop[ii] + count_wjets[ii] + count_qcd[ii];
    if (h_data) count_data[ii] = h_data->GetSum();
    
  }

  // Data counts (for preselection only)

  cout << endl << "--------------------------------------------------" << endl;
  if (channel == "el") cout << "*** " << DIR << ", electron+jets channel ***" ;
  else cout << "*** " << DIR << ", muon+jets channel ***" ;
  if (inSideband) cout << " (sideband counts) ***";
  cout << endl;
  cout         << "--------------------------------------------------" << endl;
  cout << setprecision(5);
  cout << "Sample               &     Preselection      &          0t          &         1t0b         &         1t1b         \\\\ " << endl;
  cout << "\\hline \\hline" << endl;
  cout << "\\ttbar (signal)      & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_tt_semiLep[ii] << " $\\pm$ " << setw(7) << err_tt_semiLep[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\ttbar (non-semilep) & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_tt_nonsemilep[ii] << " $\\pm$ " << setw(7) << err_tt_nonsemilep[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "Single top           & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_singletop[ii] << " $\\pm$ " << setw(7) << err_singletop[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "W+jets               & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_wjets[ii] << " $\\pm$ " << setw(7) << err_wjets[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "QCD                  & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_qcd[ii] << " $\\pm$ " << setw(7) << err_qcd[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline" << endl;
  cout << "Total                & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_tot[ii] << " $\\pm$ " << setw(7) << err_tot[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline \\hline" << endl;
  cout << "Data                 &            ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_data[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << "      &      ";
  }
}

void combineResults(TString channel, TString fit, TString preorpost) {

  TH1::AddDirectory(kFALSE); 
  setStyle();

  int nhist = 3;
  if (fit.Contains("C")) nhist = 4;
  TString whatA[3] = {"lepAbsEta0t","ak8jetTau321t0b","ak8jetSDmass1t1b"};
  TString whatB[3] = {"lepAbsEta0t","lepAbsEta1t0b","ak8jetSDmass1t1b"};
  TString whatC[4] = {"lepAbsEta0t","lepAbsEta1t0b","ak8jetTau321t0b","ak8jetSDmass1t1b"};
  const int nsample = 4;
  TString samples[nsample] = {"TTbar","SingleTop","WJets","QCD"};
  
  // counts for cutflow
  float count_tt[4] = {0};
  float count_singletop[4] = {0};
  float count_wjets[4] = {0};
  float count_qcd[4] = {0};
  float count_tot[4] = {0};
  float count_data[4] = {0};
  
  // errors for cutflow
  float err_tt[4] = {0};
  float err_singletop[4] = {0};
  float err_wjets[4] = {0};
  float err_qcd[4] = {0};
  float err_tot[4] = {0};

  TH1F* h_data[4];

  // Get data hists
  TFile* datafile = TFile::Open("combineInputs.root");
  for (int ii = 0; ii < nhist; ii++){
    if (fit.Contains("A")) h_data[ii] = (TH1F*) datafile->Get(whatA[ii]+"_"+channel+"/data_obs");
    if (fit.Contains("B")) h_data[ii] = (TH1F*) datafile->Get(whatB[ii]+"_"+channel+"/data_obs");
    if (fit.Contains("C")) h_data[ii] = (TH1F*) datafile->Get(whatC[ii]+"_"+channel+"/data_obs");
    count_data[ii] = h_data[ii]->GetSum();
  }
  datafile->Close();
  delete datafile;

  // Get sample hists
  TFile* MCfile = TFile::Open( "mlfit_"+fit+".root" );

  TString DIR = "shapes_fit_s/";
  if (preorpost == "pre") DIR = "shapes_prefit/";

  // Get count and error for each sample
  for (int ii = 0; ii < nhist; ii ++){
    TString thishist = "";
    if (fit.Contains("A")) thishist = whatA[ii];
    if (fit.Contains("B")) thishist = whatB[ii];
    if (fit.Contains("C")) thishist = whatC[ii];
    
    TH1F* h_wjets = (TH1F*) h_data[ii]->Clone();
    h_wjets->Reset();
    TH1F* h_ttbar = (TH1F*) h_data[ii]->Clone();
    h_ttbar->Reset();
    TH1F* h_singletop = (TH1F*) h_data[ii]->Clone();
    h_singletop->Reset();
    TH1F* h_qcd = (TH1F*) h_data[ii]->Clone();
    h_qcd->Reset();
    TH1F* h_totalbkg = (TH1F*) h_data[ii]->Clone();
    h_totalbkg->Reset();

    TH1F* h_wjets_tmp = (TH1F*) MCfile->Get(DIR+thishist+"_"+channel+"/WJets");
    TH1F* h_ttbar_tmp = (TH1F*) MCfile->Get(DIR+thishist+"_"+channel+"/TTbar");
    TH1F* h_singletop_tmp = (TH1F*) MCfile->Get(DIR+thishist+"_"+channel+"/SingleTop");
    TH1F* h_qcd_tmp = (TH1F*) MCfile->Get(DIR+thishist+"_"+channel+"/QCD");
    TH1F* h_totalbkg_tmp = (TH1F*) MCfile->Get(DIR+thishist+"_"+channel+"/total");

    // error on pre-fit counts
    for (int ib=0; ib<h_data[ii]->GetNbinsX(); ib++) {
      h_wjets->SetBinContent(ib+1,h_wjets_tmp->GetBinContent(ib+1));
      h_ttbar->SetBinContent(ib+1,h_ttbar_tmp->GetBinContent(ib+1));
      h_singletop->SetBinContent(ib+1,h_singletop_tmp->GetBinContent(ib+1));
      h_qcd->SetBinContent(ib+1,h_qcd_tmp->GetBinContent(ib+1));
      h_totalbkg->SetBinContent(ib+1,h_totalbkg_tmp->GetBinContent(ib+1));
      h_totalbkg->SetBinError(ib+1,h_totalbkg_tmp->GetBinError(ib+1));
      
      err_tt[ii]        += h_ttbar_tmp->GetBinError(ib+1)*h_ttbar_tmp->GetBinError(ib+1);
      err_singletop[ii] += h_singletop_tmp->GetBinError(ib+1)*h_singletop_tmp->GetBinError(ib+1);
      err_wjets[ii]     += h_wjets_tmp->GetBinError(ib+1)*h_wjets_tmp->GetBinError(ib+1);
      err_qcd[ii]       += h_qcd_tmp->GetBinError(ib+1)*h_qcd_tmp->GetBinError(ib+1);
    }

    err_tt[ii]        = sqrt(err_tt[ii]);
    err_singletop[ii] = sqrt(err_singletop[ii]);
    err_wjets[ii]     = sqrt(err_wjets[ii]);
    err_qcd[ii]       = sqrt(err_qcd[ii]);
    
    err_tot[ii] = err_tt[ii]*err_tt[ii] + err_singletop[ii]*err_singletop[ii] + err_wjets[ii]*err_wjets[ii] + err_qcd[ii]*err_qcd[ii];
    err_tot[ii] = sqrt(err_tot[ii]);

    count_tt[ii]        = h_ttbar->GetSum();
    count_singletop[ii] = h_singletop->GetSum();
    count_wjets[ii]     = h_wjets->GetSum();
    count_qcd[ii]       = h_qcd->GetSum();
    count_tot[ii]       = count_tt[ii] + count_singletop[ii] + count_wjets[ii] + count_qcd[ii];

    // create stack & summed histogram for ratio plot
    h_qcd->SetFillColor(kYellow);
    h_wjets->SetFillColor(kGreen-3);
    h_singletop->SetFillColor(6);
    h_ttbar->SetFillColor(kRed+1);
    THStack* h_stack = new THStack();    
    h_stack->Add(h_qcd);
    h_stack->Add(h_wjets);
    h_stack->Add(h_singletop);
    h_stack->Add(h_ttbar);

    h_data[ii]->SetBinErrorOption(TH1::kPoisson);

    TH1F* h_ratio;
    TH1F* h_ratio2;
    h_ratio = (TH1F*) h_data[ii]->Clone("ratio_"+thishist+"_"+channel);  // Data / MC
    h_ratio->Sumw2();
    h_ratio->Divide(h_totalbkg);

    h_ratio2 = (TH1F*) h_totalbkg->Clone("ratio2_"+thishist+"_"+channel); // Uncertainty on Data / MC
    h_ratio2->Sumw2();
    for (int ib=0; ib<h_totalbkg->GetNbinsX(); ib++) {
      h_ratio2->SetBinContent(ib+1, 1.0);
      float tmperr = h_totalbkg->GetBinError(ib+1);
      float tmpcount = h_totalbkg->GetBinContent(ib+1);
      if (tmpcount <= 0.00001) h_ratio2->SetBinError(ib+1,0.0);
      else h_ratio2->SetBinError(ib+1,tmperr/tmpcount);
    }

    float mymax = max(h_data[ii]->GetMaximum(),h_totalbkg->GetMaximum());
    h_data[ii]->SetAxisRange(0,mymax*1.05,"Y");

    // -------------------------------------------------------------------------------------
    // plotting!

    TCanvas* c = new TCanvas("c_"+thishist+"_"+channel,"c_"+thishist+"_"+channel,900,800);
    TPad* p1 = new TPad("datamcp1_"+thishist+"_"+channel,"datamcp1_"+thishist+"_"+channel,0,0.3,1,1);
    p1->SetTopMargin(0.08);
    p1->SetBottomMargin(0.05);
    p1->SetNumber(1);
    TPad* p2 = new TPad("datamcp2_"+thishist+"_"+channel,"datamcp2_"+thishist+"_"+channel,0,0,1,0.3);
    p2->SetNumber(2);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.35);
    
    p1->Draw();
    p2->Draw();
    p1->cd();
    
    h_totalbkg->UseCurrentStyle();
    h_totalbkg->SetFillColor(0);
    h_totalbkg->SetLineWidth(1);
    h_totalbkg->SetLineColor(1);
    
    h_ratio2->SetMarkerSize(0);
    h_ratio2->SetLineColor(0);
    h_ratio2->SetFillColor(15);
    h_ratio2->SetFillStyle(1001);
    
    h_ratio->GetXaxis()->SetTitle(h_data[ii]->GetXaxis()->GetTitle());
    
    h_data[ii]->UseCurrentStyle();
    h_data[ii]->SetMarkerColor(1);
    h_data[ii]->SetMarkerStyle(8);
    h_data[ii]->SetMarkerSize(1);
    h_data[ii]->GetYaxis()->SetLabelSize(26);
    h_data[ii]->GetYaxis()->SetTitleSize(32);
    h_data[ii]->GetYaxis()->SetTitleOffset(1.4);
    h_data[ii]->GetXaxis()->SetTitle("");
    
    h_data[ii]->Draw("LE0P");
    h_totalbkg->Draw("hist,same");
    h_stack->Draw("hist,same");
    h_data[ii]->Draw("LE0P,same");

    float xmin = 0.71;
    float ymin = 0.57;
    float xwidth = 0.18;
    float ywidth = 0.32;
    
    if (thishist == "ak8jetTau321t0b" || thishist == "ak8jetSDmass1t1b") xmin = 0.16;

    // legend
    TLegend* leg;
    leg = new TLegend(xmin,ymin,xmin+xwidth,ymin+ywidth);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045);
    leg->AddEntry(h_data[ii], "Data", "pe");
    leg->AddEntry(h_ttbar, "t#bar{t}", "f");
    leg->AddEntry(h_singletop, "Single t", "f");
    leg->AddEntry(h_wjets, "W+jets", "f");
    leg->AddEntry(h_qcd, "Multijet" , "f");
    leg->AddEntry(h_ratio2, "MC Stat. Unc.","f");
    leg->Draw();

    myText(0.10,0.94,1,"#intLdt = 2.7 fb^{-1}");
    myText(0.80,0.94,1,"#sqrt{s} = 13 TeV");

    // plot ratio part
    p2->cd();
    p2->SetGridy();
    h_ratio->UseCurrentStyle();
    h_ratio->SetMarkerStyle(8);
    h_ratio->SetMarkerSize(1);
    h_ratio->Draw("le0p");
    h_ratio2->Draw("same,e2");
    h_ratio->Draw("le0p,same");
    h_ratio->SetMaximum(1.8);
    h_ratio->SetMinimum(0.2);
    h_ratio->GetYaxis()->SetNdivisions(2,4,0,false);
    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetXaxis()->SetLabelSize(26);
    h_ratio->GetYaxis()->SetLabelSize(26);
    h_ratio->GetXaxis()->SetTitleOffset(2.8);
    h_ratio->GetYaxis()->SetTitleOffset(1.4);
    h_ratio->GetXaxis()->SetTitleSize(32);
    h_ratio->GetYaxis()->SetTitleSize(32);

    // save output
    TString outname = "Plots/"+channel+"_"+thishist+"_"+fit+"_"+preorpost+".pdf";
    c->SaveAs(outname);

  }

  // Now plot correlation matrix
  if (preorpost == "post"){
    TH2F* h_mcorr = (TH2F*) MCfile->Get("covariance_fit_s");
    TCanvas* c2 = new TCanvas("c2","c2",900,700);
    c2->SetLeftMargin(0.18);
    c2->SetRightMargin(0.12);
    c2->SetBottomMargin(0.1);
    h_mcorr->GetZaxis()->SetLabelSize(0.04);
    h_mcorr->Draw("COLZ");
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)h_mcorr->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.88);
    palette->SetX2NDC(0.93);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
    gPad->Modified();
    gPad->Update();
    c2->SaveAs("Plots/correlation_"+fit+".pdf");
  }

  // Finally make table
  cout << endl << "--------------------------------------------------" << endl;
  if (preorpost == "post"){
    if (channel == "el") cout << "*** Postfit event counts, electron+jets channel ***" ;
    else cout << "*** Postfit event counts, muon+jets channel ***" ;
  }
  else {
    if (channel == "el") cout << "*** Prefit event counts, electron+jets channel ***" ;
    else cout << "*** Prefit event counts, muon+jets channel ***" ;
  }
  cout << endl;
  cout         << "--------------------------------------------------" << endl;
  cout << setprecision(5);
  cout << "Sample               &          0t          &         1t0b         &         1t1b         \\\\ " << endl;
  cout << "\\hline \\hline" << endl;
  cout << "\\ttbar               & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_tt[ii] << " $\\pm$ " << setw(7) << err_tt[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "Single top           & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_singletop[ii] << " $\\pm$ " << setw(7) << err_singletop[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "W+jets               & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_wjets[ii] << " $\\pm$ " << setw(7) << err_wjets[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "QCD                  & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_qcd[ii] << " $\\pm$ " << setw(7) << err_qcd[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline" << endl;
  cout << "Total                & ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_tot[ii] << " $\\pm$ " << setw(7) << err_tot[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << " & ";
  }
  cout << "\\hline \\hline" << endl;
  cout << "Data                 &            ";
  for (int ii = 0; ii < nhist; ii++){
    cout << setw(7) << count_data[ii];
    if (ii == nhist - 1) cout << " \\\\ " << endl;
    else cout << "      &      ";
  }
}

void make2DScans(TString DIR, TString channel) {

  TH1::AddDirectory(kFALSE); 
  setStyle();

  TString hists[7] = {"2DisoPt15","2DisoPt20","2DisoPt25","2DisoPt30","2DisoPt35","2DisoPt40","2DisoPt45"};
  
  // get histograms
  for (int ii = 0; ii < 7; ii ++){
    TH2F* h_signal = (TH2F*) getSignal( DIR, hists[ii], "", channel, false );
    TH2F* h_bkg = (TH2F*) getBackground( DIR, hists[ii], "", channel, false);

    //Construct test statistic per bin
    TH2F* h_sigma = (TH2F*) h_signal->Clone();
    h_sigma->Reset();
    int nbinsx = h_sigma->GetNbinsX();
    int nbinsy = h_sigma->GetNbinsY();
    for (int jj = 1; jj < nbinsx; jj++){
      for (int kk = 1; kk < nbinsy; kk++){
	float n_sig = h_signal->GetSum() - h_signal->Integral(1,jj,1,kk); //Number of signal events
	float n_bkg = h_bkg->GetSum() - h_bkg->Integral(1,jj,1,kk); //Number of background events
	float sigma = n_sig / sqrt(n_sig + n_bkg);
	h_sigma->SetBinContent(jj,kk,sigma);
      }
    }
    
    // Plot
    TCanvas* c = new TCanvas("c_"+hists[ii],"c_"+hists[ii],900,700);
    c->SetRightMargin(1.2);
    h_signal->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_signal.pdf");
    h_bkg->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_bkg.pdf");
    h_sigma->Draw("COLZ");
    c->SaveAs("Plots/"+channel+hists[ii]+"_sigma.pdf");

  }
}

void make1DScans(TString DIR, TString channel) {

  TH1::AddDirectory(kFALSE); 
  setStyle();
  
  TString hists[5] = {"metPt","ht","htLep","lepBJetdR","lepTJetdR"};
  
  // get histograms
  for (int ii = 0; ii < 5; ii ++){
    TH1F* h_signal = (TH1F*) getSignal( DIR, hists[ii], "Pre", channel, false );
    cout << "Got signal" << endl;
    TH1F* h_bkg = (TH1F*) getBackground( DIR, hists[ii], "Pre", channel, false);
    cout << "Got background" << endl;

    //Construct test statistic per bin
    TH1F* h_sigma = (TH1F*) h_signal->Clone();
    h_sigma->Reset();
    int Nbins = h_sigma->GetNbinsX();
    for (int jj = 1; jj < Nbins+1; jj++){
      float n_sig = h_signal->Integral(jj,Nbins+1); //Number of signal events
      float n_bkg = h_bkg->Integral(jj,Nbins+1); //Number of background events
      if (hists[ii] == "lepBJetdRPre"){
	n_sig = h_signal->Integral(0,jj-1); //Number of signal events
	n_bkg = h_bkg->Integral(0,jj-1); //Number of background events
      }

      float sigma = 0.0;
      if (n_sig > 0.0 && n_bkg > 0.0) sigma = n_sig / sqrt(n_sig + n_bkg);
      h_sigma->SetBinContent(jj,sigma);
    }
    
    // Plot
    TCanvas* c = new TCanvas("c_"+hists[ii],"c_"+hists[ii],900,700);
    h_sigma->GetYaxis()->SetTitle("S / sqrt(S+B)");
    h_sigma->GetYaxis()->SetTitleOffset(1.0);
    h_sigma->Draw("hist");
    c->SaveAs("Plots/"+channel+hists[ii]+"_sigma.pdf");

  }
}

void makeEffPlots(TString channel, TString which) {
  setStyle();
  TH1::AddDirectory(kFALSE); 

  const int niso = 6;
  TString isos[niso] = {"2DisoPt25","2DisoPt45","2DisoB2G","2DisoIHNY","MiniIso10","MiniIso20"};
  int colors[niso] = {1,2,3,4,6,7};
  TH1F* h_signal_eff[niso];
  TH1F* h_bkg_eff[niso];
  
  for (int ii = 0; ii < niso; ii ++){
    // Get signal
    TH1F* h_signal_pre = (TH1F*) getSignal("histfiles_mMu_mEl_Loose/", which, "Pre", channel, false ); //These should be for no iso
    TH1F* h_signal_post = (TH1F*) getSignal("histfiles_mMu_mEl_"+isos[ii]+"/", which, "Pre", channel, false ); //These should be with iso

    // Get QCD
    SummedHist* bkg_pre = getQCDMC("histfiles_mMu_mEl_Loose/",which, "Pre",channel,false,"nom");
    SummedHist* bkg_post = getQCDMC("histfiles_mMu_mEl_"+isos[ii]+"/",which, "Pre",channel,false,"nom");
    TH1F* h_bkg_pre = (TH1F*) bkg_pre->hist();
    TH1F* h_bkg_post = (TH1F*) bkg_post->hist();

    h_bkg_pre->SetFillColor(0);
    h_bkg_post->SetFillColor(0);

    h_signal_pre->Rebin(5);
    h_bkg_pre->Rebin(5);
    h_signal_post->Rebin(5);
    h_bkg_post->Rebin(5);
     
    //Construct signal efficiency
    h_signal_eff[ii] = (TH1F*) h_signal_post->Clone();
    h_signal_eff[ii]->Divide(h_signal_pre);
    if (which == "ak8jetPt") h_signal_eff[ii]->GetXaxis()->SetRangeUser(400.,1200.);
    if (which == "leadLepPt") h_signal_eff[ii]->GetXaxis()->SetRangeUser(50.,500.);
    h_bkg_eff[ii] = (TH1F*) h_bkg_post->Clone();
    h_bkg_eff[ii]->Divide(h_bkg_pre);
    if (which == "ak8jetPt") h_bkg_eff[ii]->GetXaxis()->SetRangeUser(400.,1200.);
    if (which == "leadLepPt") h_bkg_eff[ii]->GetXaxis()->SetRangeUser(50.,500.);
  }

  // Plot
  TCanvas* c = new TCanvas("c","c",900,700);

  float xcor = 0.2;
  if (which == "leadLepPt") xcor = 0.7;  
  TLegend* leg1 = new TLegend(xcor,0.15,xcor+0.2,0.45);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.045);
  leg1->AddEntry(h_signal_eff[0], isos[0], "l");

  h_signal_eff[0]->SetLineColor(colors[0]);
  h_signal_eff[0]->GetYaxis()->SetRangeUser(0.5,1.0);
  h_signal_eff[0]->GetYaxis()->SetTitle("Signal (semilep ttbar) efficiency");
  h_signal_eff[0]->GetXaxis()->SetTitleOffset(1.1);
  h_signal_eff[0]->Draw("hist");
  
  for (int jj = 1; jj < niso; jj ++){
    h_signal_eff[jj]->SetLineColor(colors[jj]);
    h_signal_eff[jj]->Draw("hist,same");
    leg1->AddEntry(h_signal_eff[jj], isos[jj], "l");
  }
  leg1->Draw();
  c->SaveAs("Plots/eff_sig_"+which+"_"+channel+".pdf");

  TLegend* leg2 = new TLegend(xcor,0.58,xcor+0.2,0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.045);
  leg2->AddEntry(h_bkg_eff[0], isos[0], "l");

  h_bkg_eff[0]->SetLineColor(colors[0]);
  h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,1.0);
  h_bkg_eff[0]->GetYaxis()->SetTitle("Background (QCD) efficiency");
  h_bkg_eff[0]->GetXaxis()->SetTitleOffset(1.1);
  if (channel == "mu" && which == "ak8jetPt")  h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,0.1);
  h_bkg_eff[0]->Draw("hist");
  
  for (int jj = 1; jj < niso; jj ++){
    h_bkg_eff[jj]->SetLineColor(colors[jj]);
    h_bkg_eff[jj]->Draw("hist,same");
    leg2->AddEntry(h_bkg_eff[jj], isos[jj], "l");
  }
  leg2->Draw();
  c->SaveAs("Plots/eff_bkg_"+which+"_"+channel+".pdf");
}

void makeEffPlotsFinal() {
  setStyle();
  TH1::AddDirectory(kFALSE); 

  TString channels[2] = {"mu","el"};
  TString hists[2] = {"ak8jetPt","leadLepPt"};
  
  for (int ii = 0; ii < 2; ii++){
    TH1F* h_signal_eff[2];
    TH1F* h_bkg_eff[2];
    for (int jj = 0; jj < 2; jj++){
      // get histograms
      TH1F* h_signal_pre = (TH1F*) getSignal("histfiles_mMu_mEl_Loose/", hists[ii], "Pre", channels[jj], false );
      h_signal_pre->Sumw2();
      TH1F* h_signal_post = (TH1F*) getSignal("histfiles_mMu_tEl_MiniIso10/", hists[ii], "Pre", channels[jj], false );
      h_signal_post->Sumw2();

      //This version uses QCD only
      SummedHist* bkg_pre = getQCDMC("histfiles_mMu_mEl_Loose/",hists[ii],"Pre",channels[jj],false,"nom");
      SummedHist* bkg_post = getQCDMC("histfiles_mMu_tEl_MiniIso10/",hists[ii],"Pre",channels[jj],false,"nom");
      TH1F* h_bkg_pre = (TH1F*) bkg_pre->hist();
      h_bkg_pre->Sumw2();
      TH1F* h_bkg_post = (TH1F*) bkg_post->hist();
      h_bkg_post->Sumw2();

      h_bkg_pre->SetFillColor(0);
      h_bkg_post->SetFillColor(0);

      h_signal_pre->Rebin(5);
      h_bkg_pre->Rebin(5);
      h_signal_post->Rebin(5);
      h_bkg_post->Rebin(5);
     
      //Construct signal efficiency
      h_signal_eff[jj] = (TH1F*) h_signal_post->Clone();
      h_signal_eff[jj]->SetName("signal_eff_"+channels[jj]);
      h_signal_eff[jj]->GetYaxis()->SetTitle("Signal Efficiency");
      h_signal_eff[jj]->Divide(h_signal_post,h_signal_pre,1.0,1.0,"B");
      if (hists[ii] == "ak8jetPt") h_signal_eff[jj]->GetXaxis()->SetRangeUser(400.,1200.);
      if (hists[ii] == "leadLepPt") h_signal_eff[jj]->GetXaxis()->SetRangeUser(50.,500.);
      h_bkg_eff[jj] = (TH1F*) h_bkg_post->Clone();
      h_bkg_eff[jj]->SetName("bkg_eff_"+channels[jj]);
      h_bkg_eff[jj]->GetYaxis()->SetTitle("Background Efficiency");
      h_bkg_eff[jj]->Divide(h_bkg_post,h_bkg_pre,1.0,1.0,"B");
      if (hists[ii] == "ak8jetPt") h_bkg_eff[jj]->GetXaxis()->SetRangeUser(400.,1200.);
      if (hists[ii] == "leadLepPt") h_bkg_eff[jj]->GetXaxis()->SetRangeUser(50.,500.);
    }//End channel loop
    
    // Plot
    TCanvas* c = new TCanvas("c","c",900,700);

    float xcor = 0.2;
    if (hists[ii] == "leadLepPt") xcor = 0.7;  
    TLegend* leg1 = new TLegend(xcor,0.15,xcor+0.2,0.45);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.045);
    leg1->AddEntry(h_signal_eff[0], "Eff_sig, mu", "l");
    leg1->AddEntry(h_signal_eff[1], "Eff_sig, el", "l");
    
    h_signal_eff[0]->SetLineColor(1);
    h_signal_eff[0]->GetYaxis()->SetRangeUser(0.5,1.0);
    h_signal_eff[0]->GetXaxis()->SetTitleOffset(1.1);
    h_signal_eff[0]->Draw("e");
  
    h_signal_eff[1]->SetLineColor(2);
    h_signal_eff[1]->Draw("e,same");
    leg1->Draw();
    c->SaveAs("Plots/eff_sig_final_"+hists[ii]+".pdf");

    TLegend* leg2 = new TLegend(xcor,0.58,xcor+0.2,0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.045);
    leg2->AddEntry(h_bkg_eff[0], "Eff_bkg, mu", "l");
    leg2->AddEntry(h_bkg_eff[1], "Eff_bkg, el", "l");

    h_bkg_eff[0]->SetLineColor(1);
    h_bkg_eff[0]->GetYaxis()->SetRangeUser(0.0,0.5);
    h_bkg_eff[0]->GetXaxis()->SetTitleOffset(1.1);
    h_bkg_eff[0]->Draw("e");
  
    h_bkg_eff[1]->SetLineColor(2);
    h_bkg_eff[1]->Draw("e,same");
    leg2->Draw();
    c->SaveAs("Plots/eff_bkg_final_"+hists[ii]+".pdf");
  }
}

void drawROCCurve(TString channel){

  setStyle();
  
  SummedHist* Iso2DScanPoints_sig = getTTbar("histfiles_mMu_mEl_Loose/","2DisoScanPoints","",channel,false,"nom");
  SummedHist* Iso2DScanPoints_bkg = getQCDMC("histfiles_mMu_mEl_Loose/","2DisoScanPoints","",channel,false,"nom");
  SummedHist* MiniIsoScanPoints_sig = getTTbar("histfiles_mMu_mEl_Loose/","MiniIsoScanPoints","",channel,false,"nom");
  SummedHist* MiniIsoScanPoints_bkg = getQCDMC("histfiles_mMu_mEl_Loose/","MiniIsoScanPoints","",channel,false,"nom");
  SummedHist* Iso2DScanPoints_sig_tight = getTTbar("histfiles_tMu_tEl/","2DisoScanPoints","",channel,false,"nom");
  SummedHist* Iso2DScanPoints_bkg_tight = getQCDMC("histfiles_tMu_tEl/","2DisoScanPoints","",channel,false,"nom");
  SummedHist* MiniIsoScanPoints_sig_tight = getTTbar("histfiles_tMu_tEl/","MiniIsoScanPoints","",channel,false,"nom");
  SummedHist* MiniIsoScanPoints_bkg_tight = getQCDMC("histfiles_tMu_tEl/","MiniIsoScanPoints","",channel,false,"nom");
    
  TH1F* h_Iso2DScanPoints_sig = (TH1F*) Iso2DScanPoints_sig->hist();
  TH1F* h_Iso2DScanPoints_bkg = (TH1F*) Iso2DScanPoints_bkg->hist();
  TH1F* h_MiniIsoScanPoints_sig = (TH1F*) MiniIsoScanPoints_sig->hist();
  TH1F* h_MiniIsoScanPoints_bkg = (TH1F*) MiniIsoScanPoints_bkg->hist();
  TH1F* h_Iso2DScanPoints_sig_tight = (TH1F*) Iso2DScanPoints_sig_tight->hist();
  TH1F* h_Iso2DScanPoints_bkg_tight = (TH1F*) Iso2DScanPoints_bkg_tight->hist();
  TH1F* h_MiniIsoScanPoints_sig_tight = (TH1F*) MiniIsoScanPoints_sig_tight->hist();
  TH1F* h_MiniIsoScanPoints_bkg_tight = (TH1F*) MiniIsoScanPoints_bkg_tight->hist();
  
  const int n2D = 168;
  const int nMini = 30;
  
  float Sig_eff_2D[n2D];
  float Bkg_eff_2D[n2D];
  float Sig_eff_2D_tight[n2D];
  float Bkg_eff_2D_tight[n2D];
  float Mini_sig_eff[nMini];
  float Mini_bkg_eff[nMini];
  float Mini_sig_eff_tight[nMini];
  float Mini_bkg_eff_tight[nMini];
  
  for (int i=0; i<n2D; i++) {
    Sig_eff_2D[i] = h_Iso2DScanPoints_sig->GetBinContent(i+2) / h_Iso2DScanPoints_sig->GetBinContent(1);
    Bkg_eff_2D[i] = h_Iso2DScanPoints_bkg->GetBinContent(i+2) / h_Iso2DScanPoints_bkg->GetBinContent(1);
    Sig_eff_2D_tight[i] = h_Iso2DScanPoints_sig_tight->GetBinContent(i+2) / h_Iso2DScanPoints_sig_tight->GetBinContent(1);
    Bkg_eff_2D_tight[i] = h_Iso2DScanPoints_bkg_tight->GetBinContent(i+2) / h_Iso2DScanPoints_bkg_tight->GetBinContent(1);
  }
  
  for (int i=0; i<nMini; i++) {
    Mini_sig_eff[i] = h_MiniIsoScanPoints_sig->GetBinContent(i+2) / h_MiniIsoScanPoints_sig->GetBinContent(1);
    Mini_bkg_eff[i] = h_MiniIsoScanPoints_bkg->GetBinContent(i+2) / h_MiniIsoScanPoints_bkg->GetBinContent(1);
    Mini_sig_eff_tight[i] = h_MiniIsoScanPoints_sig_tight->GetBinContent(i+2) / h_MiniIsoScanPoints_sig_tight->GetBinContent(1);
    Mini_bkg_eff_tight[i] = h_MiniIsoScanPoints_bkg_tight->GetBinContent(i+2) / h_MiniIsoScanPoints_bkg_tight->GetBinContent(1);
  }
    
  TGraph* g_eff2D = new TGraph(n2D,Sig_eff_2D,Bkg_eff_2D);
  TGraph* g_effMini = new TGraph(nMini,Mini_sig_eff,Mini_bkg_eff);
  TGraph* g_eff2D_tight = new TGraph(n2D,Sig_eff_2D_tight,Bkg_eff_2D_tight);
  TGraph* g_effMini_tight = new TGraph(nMini,Mini_sig_eff_tight,Mini_bkg_eff_tight);

  TH2F* h_dummy_mu = new TH2F("dummy_mu", "; efficiency (ttbar); efficiency (QCD)",50,0.5,1.0,100,0.0,0.1);
  TH2F* h_dummy_el = new TH2F("dummy_el", "; efficiency (ttbar); efficiency (QCD)",50,0.5,1.0,100,0.0,1.0);
  
  g_eff2D->SetMarkerStyle(8);
  g_effMini->SetMarkerColor(2);
  g_effMini->SetLineColor(2);
  g_effMini->SetMarkerStyle(8);
  g_eff2D_tight->SetMarkerColor(3);
  g_eff2D_tight->SetLineColor(3);
  g_eff2D_tight->SetMarkerStyle(8);
  g_effMini_tight->SetMarkerColor(4);
  g_effMini_tight->SetLineColor(4);
  g_effMini_tight->SetMarkerStyle(8);
  
  TCanvas* c = new TCanvas("c","c",900,700);
  if (channel == "mu") h_dummy_mu->Draw();
  if (channel == "el") h_dummy_el->Draw();
  g_eff2D->Draw("same,ep");
  g_effMini->Draw("same,ep");
  g_eff2D_tight->Draw("same,ep");
  g_effMini_tight->Draw("same,ep");
  
  TLegend* lh = new TLegend(0.20,0.60,0.45,0.8);
  lh->SetFillStyle(0);
  lh->SetBorderSize(0);
  lh->SetTextSize(0.045);
  lh->AddEntry(g_eff2D,"2D iso, medium ID","lp");
  lh->AddEntry(g_effMini,"MiniIso, medium ID","lp");
  lh->AddEntry(g_eff2D_tight,"2D iso, tight ID","lp");
  lh->AddEntry(g_effMini_tight,"MiniIso, tight ID","lp");
  lh->SetTextFont(42);
  lh->Draw();	
  
  c->SaveAs("ROCCurve_"+channel+".pdf");
  
}



