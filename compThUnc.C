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

void compThUnc() {
  
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  const int nHIST = 2;
  TString histnames[nHIST] = {"ptGenTop","ptRecoTop"};

  const int nCH = 2;
  TString channels[nCH] = {"mu","el"};

  const int nSYS = 7;
  TString systnames[nSYS] = {"nom","PDFUp","PDFDown","Q2Up","Q2Down","ASUp","ASDown"};
  int colors[nSYS] = {1,4,4,6,6,8,8};
  int styles[nSYS] = {1,1,2,1,2,1,2};

  TH1D* hists[nHIST][nCH][nSYS];

  for (int ii = 0; ii < nCH; ii++){
    for (int jj = 0; jj < nSYS; jj++){
      TFile* f_sys = TFile::Open("histfiles_80X/hists_PowhegPythia8_fullTruth_" + channels[ii] + "_" + systnames[jj] + ".root");
      for (int kk = 0; kk < nHIST; kk++){
	hists[kk][ii][jj] = (TH1D*) f_sys->Get(histnames[kk]);
      }
      f_sys->Close();
      delete f_sys;
    }
  }

  double normfactors[nCH][nSYS];
  for (int ii = 0; ii < nCH; ii++){
    for (int jj = 0; jj < nSYS; jj++){
      normfactors[ii][jj] = hists[0][ii][0]->Integral() / hists[0][ii][jj]->Integral();
      cout << "Acceptance for " << systnames[jj] << " in " << channels[ii] << " channel is " << hists[1][ii][jj]->Integral() / hists[0][ii][jj]->Integral() << endl;
    }
  }

  for (int ii = 0; ii < nCH; ii++){
    for (int jj = 0; jj < nHIST; jj++){
      TH1D* h_ratio[nSYS];
      for (int kk = 0; kk < nSYS; kk++){
	hists[jj][ii][kk]->SetLineColor(colors[kk]);
	hists[jj][ii][kk]->SetLineStyle(styles[kk]);
	hists[jj][ii][kk]->SetFillColor(0);
	hists[jj][ii][kk]->Sumw2();
	hists[jj][ii][kk]->Scale(normfactors[ii][kk]);
	h_ratio[kk] = (TH1D*) hists[jj][ii][kk]->Clone();
	h_ratio[kk]->Divide(hists[jj][ii][0]);
	h_ratio[kk]->SetLineColor(colors[kk]);
	h_ratio[kk]->SetLineStyle(styles[kk]);
	h_ratio[kk]->SetFillColor(0);
      }
      
      TCanvas* c = new TCanvas("c_"+histnames[jj]+"_"+channels[ii],"c_"+histnames[jj]+"_"+channels[ii],900,800);
      TPad* p1 = new TPad("datamcp1_"+histnames[jj]+"_"+channels[ii],"datamcp1_"+histnames[jj]+"_"+channels[ii],0,0.3,1,1);
      p1->SetTopMargin(0.08);
      p1->SetBottomMargin(0.05);
      p1->SetNumber(1);
      TPad* p2 = new TPad("datamcp2_"+histnames[jj]+"_"+channels[ii],"datamcp2_"+histnames[jj]+"_"+channels[ii],0,0,1,0.3);
      p2->SetNumber(2);
      p2->SetTopMargin(0.05);
      p2->SetBottomMargin(0.35);
      
      p1->Draw();
      p2->Draw();
      p1->cd();
      p1->SetLogy();
      
      h_ratio[0]->GetXaxis()->SetTitle(hists[jj][ii][0]->GetXaxis()->GetTitle());
      hists[jj][ii][0]->GetXaxis()->SetTitle("");
      hists[jj][ii][0]->SetMaximum(1.5 * hists[jj][ii][0]->GetMaximum());
      hists[jj][ii][0]->GetYaxis()->SetTitleOffset(1.0);
      if (histnames[jj] == "ptRecoTop") hists[jj][ii][0]->GetXaxis()->SetRangeUser(400.,1200.);
      hists[jj][ii][0]->Draw("hist");

      float xmin = 0.66;
      float xmax = 0.85;
      float ymin = 0.60;
      float ymax = 0.88;
      
      // legend
      TLegend* leg;
      leg = new TLegend(xmin,ymin,xmax,ymax);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.045);
      leg->AddEntry(hists[jj][ii][0], "Nominal", "l");

      for (int kk = 1; kk < nSYS; kk++){
	hists[jj][ii][kk]->Draw("hist,same");
	leg->AddEntry(hists[jj][ii][kk], systnames[kk], "l");
      }
	
      leg->Draw();

      // plot ratio part
      p2->cd();

      h_ratio[0]->Draw("hist");
      for (int kk = 1; kk < nSYS; kk++){
	h_ratio[kk]->Draw("hist,same");
      }

      h_ratio[0]->SetMaximum(1.2);
      h_ratio[0]->SetMinimum(0.8);
      h_ratio[0]->GetYaxis()->SetNdivisions(2,4,0,false);
      h_ratio[0]->GetYaxis()->SetTitle("Sys Var / Nom");
      h_ratio[0]->GetXaxis()->SetLabelSize(0.1);
      h_ratio[0]->GetYaxis()->SetLabelSize(0.1);
      h_ratio[0]->GetXaxis()->SetTitleOffset(1.0);
      h_ratio[0]->GetYaxis()->SetTitleOffset(0.45);
      h_ratio[0]->GetXaxis()->SetTitleSize(0.12);
      h_ratio[0]->GetYaxis()->SetTitleSize(0.1);
      if (histnames[jj] == "ptRecoTop") h_ratio[0]->GetXaxis()->SetRangeUser(400.,1200.);

      // save output
      TString outname = "Plots/compThUnc_"+channels[ii]+"_"+histnames[jj]+".pdf";
      c->SaveAs(outname);

      for (int kk = 0; kk < nSYS; kk++){
	hists[jj][ii][kk]->Delete();
	h_ratio[kk]->Delete();
      }
    }
  }
}
