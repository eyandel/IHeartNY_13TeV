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
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void calcPUweight(){
  TFile * F_data_PUnom       = new TFile( "MyDataPileupHistogram.root"                      );
  TFile * F_MC               = new TFile( "test_iheartNY.root" );

  TH1D * NPU_data_true       = (TH1D*) F_data_PUnom->Get("pileup");
  TH1D * NPU_MC_true         = (TH1D*) F_MC->Get("h_NtrueIntPU");

  NPU_data_true       ->Sumw2();
  NPU_MC_true         ->Sumw2();

  // scale to have area 1
  NPU_data_true       ->Scale( 1.0 / NPU_data_true       ->Integral() );
  NPU_MC_true         ->Scale( 1.0 / NPU_MC_true         ->Integral() );

  // Draw NPU data and MC
  TCanvas *c1 = new TCanvas("c1","",200,10,950,700);
  gStyle->SetOptStat(0);
  gStyle->SetHistLineWidth(2);

  NPU_data_true ->SetLineColor(1);
  NPU_MC_true   ->SetLineColor(2);

  NPU_data_true ->SetMarkerColor(1);
  NPU_MC_true   ->SetMarkerColor(2);

  NPU_data_true ->SetMarkerStyle(20);
  NPU_MC_true   ->SetMarkerStyle(21);

  NPU_data_true->SetTitle(";N_{PU} true;Number of events");
  NPU_data_true->GetXaxis()->SetRangeUser(0,30);

  NPU_data_true ->Draw();
  NPU_MC_true   ->Draw("same");

  TLegend * leg = new TLegend(0.7,0.7,0.85,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle( 4050 );
  leg->AddEntry(NPU_data_true ,"data","LP");
  leg->AddEntry(NPU_MC_true   ,"MC","LP");
  leg->Draw("same");

  c1->RedrawAxis();
  c1->SaveAs("NPU.pdf");

  // calculate weight = data / MC
  NPU_data_true       ->Divide( NPU_MC_true );
  NPU_data_true      ->SetName( "PUweight_true"      ) ;
  NPU_data_true      ->SetTitle( ";N_{PU} true; Weight"      ) ;

  // Save weight in output file
  TFile *Out;
  Out = new TFile("pileup_reweight.root","RECREATE");
  Out->cd();

  NPU_data_true       ->Write();

  Out->ls();
  Out->Write();
}
