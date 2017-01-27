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
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void calcSF() {
  
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  const int nelPtbins = 8;
  float elPtbins[nelPtbins+1] = {50.,75.,100.,125.,150.,200.,250.,300.,500.};
  const int nleadJetPtbins = 7;
  float leadJetPtbins[nleadJetPtbins+1] = {250.,300.,350.,400.,500.,600.,800.,1200.};
  const int nsubJetPtbins = 6;
  float subJetPtbins[nsubJetPtbins+1] = {50.,100.,150.,200.,250.,350.,550.};
  
  TH1F* pass_elPt      = new TH1F("pass_elPt"     ,"",nelPtbins,elPtbins);
  TH1F* pass_elAbsEta  = new TH1F("pass_elAbsEta" ,"",10,0.,2.5);
  TH1F* pass_leadJetPt = new TH1F("pass_leadJetPt","",nleadJetPtbins,leadJetPtbins);
  TH1F* pass_subJetPt  = new TH1F("pass_subJetPt" ,"",nsubJetPtbins,subJetPtbins);
  TH1F* pass_dRelLead  = new TH1F("pass_dRelLead" ,"",10,0.,3.5);
  TH1F* pass_dRelSub   = new TH1F("pass_dRelSub"  ,"",10,0.,3.5);

  TH1F* tot_elPt       = new TH1F("tot_elPt"      ,"",nelPtbins,elPtbins);
  TH1F* tot_elAbsEta   = new TH1F("tot_elAbsEta"  ,"",10,0.,2.5);
  TH1F* tot_leadJetPt  = new TH1F("tot_leadJetPt" ,"",nleadJetPtbins,leadJetPtbins);
  TH1F* tot_subJetPt   = new TH1F("tot_subJetPt"  ,"",nsubJetPtbins,subJetPtbins);
  TH1F* tot_dRelLead   = new TH1F("tot_dRelLead"  ,"",10,0.,3.5);
  TH1F* tot_dRelSub    = new TH1F("tot_dRelSub"   ,"",10,0.,3.5);

  TChain* tree = new TChain("recoTree");
  tree->Add("skimTrees_80X/TnP_mu_new.root");
    
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
    
  vector<int>*   passElTrig  = 0;
  vector<float>* muPt        = 0;
  vector<float>* muEta       = 0;
  vector<float>* muPhi       = 0;
  vector<float>* elPt        = 0;
  vector<float>* elEta       = 0;
  vector<float>* elPhi       = 0;
  vector<float>* leadJetPt   = 0;
  vector<float>* leadJetEta  = 0;
  vector<float>* leadJetPhi  = 0;
  vector<float>* leadJetMass = 0;
  vector<float>* subJetPt    = 0;
  vector<float>* subJetEta   = 0;
  vector<float>* subJetPhi   = 0;
  vector<float>* subJetMass  = 0;
  
  TBranch* b_passElTrig;
  TBranch* b_muPt;
  TBranch* b_muEta;
  TBranch* b_muPhi;
  TBranch* b_elPt;
  TBranch* b_elEta;
  TBranch* b_elPhi;
  TBranch* b_leadJetPt;
  TBranch* b_leadJetEta;
  TBranch* b_leadJetPhi;
  TBranch* b_leadJetMass;
  TBranch* b_subJetPt;
  TBranch* b_subJetEta;
  TBranch* b_subJetPhi;
  TBranch* b_subJetMass;
  
  tree->SetBranchAddress("passElTrig" , &passElTrig , &b_passElTrig  );     
  tree->SetBranchAddress("muPt"       , &muPt       , &b_muPt        );     
  tree->SetBranchAddress("muEta"      , &muEta      , &b_muEta       );     
  tree->SetBranchAddress("muPhi"      , &muPhi      , &b_muPhi       );     
  tree->SetBranchAddress("elPt"       , &elPt       , &b_elPt        );     
  tree->SetBranchAddress("elEta"      , &elEta      , &b_elEta       );     
  tree->SetBranchAddress("elPhi"      , &elPhi      , &b_elPhi       );     
  tree->SetBranchAddress("leadJetPt"  , &leadJetPt  , &b_leadJetPt   );     
  tree->SetBranchAddress("leadJetEta" , &leadJetEta , &b_leadJetEta  );     
  tree->SetBranchAddress("leadJetPhi" , &leadJetPhi , &b_leadJetPhi  );     
  tree->SetBranchAddress("leadJetMass", &leadJetMass, &b_leadJetMass );     
  tree->SetBranchAddress("subJetPt"   , &subJetPt   , &b_subJetPt    );     
  tree->SetBranchAddress("subJetEta"  , &subJetEta  , &b_subJetEta   );     
  tree->SetBranchAddress("subJetPhi"  , &subJetPhi  , &b_subJetPhi   );     
  tree->SetBranchAddress("subJetMass" , &subJetMass , &b_subJetMass  );     
  
  for (int i=0; i<tree->GetEntries(); i++) {
    
    tree->GetEntry(i,0);
    
    TLorentzVector muP4;
    muP4.SetPtEtaPhiM(muPt->at(0),muEta->at(0),muPhi->at(0),0.105);
    TLorentzVector elP4;
    elP4.SetPtEtaPhiM(elPt->at(0),elEta->at(0),elPhi->at(0),0.0);
    TLorentzVector leadJetP4;
    leadJetP4.SetPtEtaPhiM(leadJetPt->at(0),leadJetEta->at(0),leadJetPhi->at(0),leadJetMass->at(0));
    TLorentzVector subJetP4;
    subJetP4.SetPtEtaPhiM(subJetPt->at(0),subJetEta->at(0),subJetPhi->at(0),subJetMass->at(0));

    tot_elPt->Fill(elP4.Perp());
    tot_elAbsEta->Fill(fabs(elP4.Eta()));
    tot_leadJetPt->Fill(leadJetP4.Perp());
    tot_subJetPt->Fill(subJetP4.Perp());
    tot_dRelLead->Fill(elP4.DeltaR(leadJetP4));
    tot_dRelSub->Fill(elP4.DeltaR(subJetP4));

    if (passElTrig->at(0) == 1){
      pass_elPt->Fill(elP4.Perp());
      pass_elAbsEta->Fill(fabs(elP4.Eta()));
      pass_leadJetPt->Fill(leadJetP4.Perp());
      pass_subJetPt->Fill(subJetP4.Perp());
      pass_dRelLead->Fill(elP4.DeltaR(leadJetP4));
      pass_dRelSub->Fill(elP4.DeltaR(subJetP4));
    }
  }

  tree->Delete();

  pass_elPt->Sumw2();
  pass_elAbsEta->Sumw2();
  pass_leadJetPt->Sumw2();
  pass_subJetPt->Sumw2();
  pass_dRelLead->Sumw2();
  pass_dRelSub->Sumw2();
  tot_elPt->Sumw2();
  tot_elAbsEta->Sumw2();
  tot_leadJetPt->Sumw2();
  tot_subJetPt->Sumw2();
  tot_dRelLead->Sumw2();
  tot_dRelSub->Sumw2();

  TH1F* eff_elPt = (TH1F*) pass_elPt->Clone();
  eff_elPt->SetName("eff_elPt");
  eff_elPt->GetYaxis()->SetTitle("Efficiency");
  eff_elPt->GetXaxis()->SetTitle("Electron p_{T}");
  eff_elPt->Divide(pass_elPt, tot_elPt, 1.0, 1.0, "B");
  TH1F* eff_elAbsEta = (TH1F*) pass_elAbsEta->Clone();
  eff_elAbsEta->SetName("eff_elAbsEta");
  eff_elAbsEta->GetYaxis()->SetTitle("Efficiency");
  eff_elAbsEta->GetXaxis()->SetTitle("Electron |eta|");
  eff_elAbsEta->Divide(pass_elAbsEta, tot_elAbsEta, 1.0, 1.0, "B");
  TH1F* eff_leadJetPt = (TH1F*) pass_leadJetPt->Clone();
  eff_leadJetPt->SetName("eff_leadJetPt");
  eff_leadJetPt->GetYaxis()->SetTitle("Efficiency");
  eff_leadJetPt->GetXaxis()->SetTitle("Leading jet p_{T}");
  eff_leadJetPt->Divide(pass_leadJetPt, tot_leadJetPt, 1.0, 1.0, "B");
  TH1F* eff_subJetPt = (TH1F*) pass_subJetPt->Clone();
  eff_subJetPt->SetName("eff_subJetPt");
  eff_subJetPt->GetYaxis()->SetTitle("Efficiency");
  eff_subJetPt->GetXaxis()->SetTitle("Subleading jet p_{T}");
  eff_subJetPt->Divide(pass_subJetPt, tot_subJetPt, 1.0, 1.0, "B");
  TH1F* eff_dRelLead = (TH1F*) pass_dRelLead->Clone();
  eff_dRelLead->SetName("eff_dRelLead");
  eff_dRelLead->GetYaxis()->SetTitle("Efficiency");
  eff_dRelLead->Divide(pass_dRelLead, tot_dRelLead, 1.0, 1.0, "B");
  TH1F* eff_dRelSub = (TH1F*) pass_dRelSub->Clone();
  eff_dRelSub->SetName("eff_dRelSub");
  eff_dRelSub->GetYaxis()->SetTitle("Efficiency");
  eff_dRelSub->Divide(pass_dRelSub, tot_dRelSub, 1.0, 1.0, "B");

  TCanvas* c = new TCanvas("c","c",900,800);
  eff_elPt->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_elPt->Draw();
  c->SaveAs("Plots/eff_elPt.pdf");
  eff_elAbsEta->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_elAbsEta->Draw();
  c->SaveAs("Plots/eff_elAbsEta.pdf");
  eff_leadJetPt->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_leadJetPt->Draw();
  c->SaveAs("Plots/eff_leadJetPt.pdf");
  eff_subJetPt->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_subJetPt->Draw();
  c->SaveAs("Plots/eff_subJetPt.pdf");
  eff_dRelLead->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_dRelLead->Draw();
  c->SaveAs("Plots/eff_dRelLead.pdf");
  eff_dRelSub->GetYaxis()->SetRangeUser(0.5,1.5);
  eff_dRelSub->Draw();
  c->SaveAs("Plots/eff_dRelSub.pdf");

  TFile* f_SF = new TFile("elTrigEff.root","recreate");
  eff_elPt->Write();
  f_SF->Close();
  delete f_SF;
  
  pass_elPt->Delete();
  pass_elAbsEta->Delete();
  pass_leadJetPt->Delete();
  pass_subJetPt->Delete();
  pass_dRelLead->Delete();
  pass_dRelSub->Delete();
  tot_elPt->Delete();
  tot_elAbsEta->Delete();
  tot_leadJetPt->Delete();
  tot_subJetPt->Delete();
  tot_dRelLead->Delete();
  tot_dRelSub->Delete();
  eff_elPt->Delete();
  eff_elAbsEta->Delete();
  eff_leadJetPt->Delete();
  eff_subJetPt->Delete();
  eff_dRelLead->Delete();
  eff_dRelSub->Delete();


}
