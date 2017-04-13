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

void calcSF(TString channel) {
  
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  const int nLepPtbins = 8;
  float lepPtbins[nLepPtbins+1] = {50.,75.,100.,125.,150.,200.,250.,300.,500.};
  const int nleadJetPtbins = 7;
  float leadJetPtbins[nleadJetPtbins+1] = {250.,300.,350.,400.,500.,600.,800.,1200.};
  const int nsubJetPtbins = 6;
  float subJetPtbins[nsubJetPtbins+1] = {50.,100.,150.,200.,250.,350.,550.};

  TH1F* h_eff[2][6];

  TString vars[6] = {"lepPt","lepAbsEta","leadJetPt","subJetPt","dRlepLead","dRlepSub"};
  TString titles[6] = {"Lepton p_{T}","Lepton |#eta|","Leading jet p_{T}","Subleading jet p_{T}","dR(Lepton, leading jet)","dR(Lepton, subleading jet)"};
  
  TString samples[2] = {"MC","Data"};
  
  for (int ii = 0; ii < 2; ii++){
    
    TH1F* h_pass[6];
    TH1F* h_tot[6];
    
    h_pass[0] = new TH1F("pass_lepPt"    ,"",nLepPtbins,lepPtbins);
    h_pass[1] = new TH1F("pass_lepAbsEta","",10,0.,2.5);
    h_pass[2] = new TH1F("pass_leadJetPt","",nleadJetPtbins,leadJetPtbins);
    h_pass[3] = new TH1F("pass_subJetPt" ,"",nsubJetPtbins,subJetPtbins);
    h_pass[4] = new TH1F("pass_dRlepLead","",10,0.,3.5);
    h_pass[5] = new TH1F("pass_dRlepSub" ,"",10,0.,3.5);
    
    h_tot[0] = new TH1F("tot_lepPt"     ,"",nLepPtbins,lepPtbins);
    h_tot[1] = new TH1F("tot_lepAbsEta" ,"",10,0.,2.5);
    h_tot[2] = new TH1F("tot_leadJetPt" ,"",nleadJetPtbins,leadJetPtbins);
    h_tot[3] = new TH1F("tot_subJetPt"  ,"",nsubJetPtbins,subJetPtbins);
    h_tot[4] = new TH1F("tot_dRlepLead" ,"",10,0.,3.5);
    h_tot[5] = new TH1F("tot_dRlepSub"  ,"",10,0.,3.5);
  
    TChain* tree = new TChain("recoTree");
    tree->Add("skimTrees_full2016/TnP_"+samples[ii]+"_"+channel+".root");
    
    if (tree->GetEntries() == 0) {
      cout << "File doesn't exist or is empty, returning..." << endl;
      return;
    }

    int nentries = tree->GetEntries();
    
    vector<int>*   passProbeTrig = 0;
    vector<float>* prescaleTag   = 0;
    vector<float>* eventWeight   = 0;
    vector<float>* muPt          = 0;
    vector<float>* muEta         = 0;
    vector<float>* muPhi         = 0;
    vector<float>* elPt          = 0;
    vector<float>* elEta         = 0;
    vector<float>* elPhi         = 0;
    vector<float>* leadJetPt     = 0;
    vector<float>* leadJetEta    = 0;
    vector<float>* leadJetPhi    = 0;
    vector<float>* leadJetMass   = 0;
    vector<float>* subJetPt      = 0;
    vector<float>* subJetEta     = 0;
    vector<float>* subJetPhi     = 0;
    vector<float>* subJetMass    = 0;
  
    TBranch* b_passProbeTrig;
    TBranch* b_prescaleTag;
    TBranch* b_eventWeight;
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
    
    tree->SetBranchAddress("passProbeTrig", &passProbeTrig, &b_passProbeTrig);
    tree->SetBranchAddress("prescaleTag"  , &prescaleTag  , &b_prescaleTag  );
    tree->SetBranchAddress("eventWeight"  , &eventWeight  , &b_eventWeight  );
    tree->SetBranchAddress("muPt"         , &muPt         , &b_muPt         );     
    tree->SetBranchAddress("muEta"        , &muEta        , &b_muEta        );     
    tree->SetBranchAddress("muPhi"        , &muPhi        , &b_muPhi        );     
    tree->SetBranchAddress("elPt"         , &elPt         , &b_elPt         );     
    tree->SetBranchAddress("elEta"        , &elEta        , &b_elEta        );     
    tree->SetBranchAddress("elPhi"        , &elPhi        , &b_elPhi        );     
    tree->SetBranchAddress("leadJetPt"    , &leadJetPt    , &b_leadJetPt    );     
    tree->SetBranchAddress("leadJetEta"   , &leadJetEta   , &b_leadJetEta   );     
    tree->SetBranchAddress("leadJetPhi"   , &leadJetPhi   , &b_leadJetPhi   );     
    tree->SetBranchAddress("leadJetMass"  , &leadJetMass  , &b_leadJetMass  );     
    tree->SetBranchAddress("subJetPt"     , &subJetPt     , &b_subJetPt     );     
    tree->SetBranchAddress("subJetEta"    , &subJetEta    , &b_subJetEta    );     
    tree->SetBranchAddress("subJetPhi"    , &subJetPhi    , &b_subJetPhi    );     
    tree->SetBranchAddress("subJetMass"   , &subJetMass   , &b_subJetMass   );     
  
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

      TLorentzVector lepP4;
      if (channel == "mu") lepP4 = muP4;
      else lepP4 = elP4;

      h_tot[0]->Fill(lepP4.Perp());
      h_tot[1]->Fill(fabs(lepP4.Eta()));
      h_tot[2]->Fill(leadJetP4.Perp());
      h_tot[3]->Fill(subJetP4.Perp());
      h_tot[4]->Fill(lepP4.DeltaR(leadJetP4));
      h_tot[5]->Fill(lepP4.DeltaR(subJetP4));
      
      if (passProbeTrig->at(0) == 1){
	h_pass[0]->Fill(lepP4.Perp());
	h_pass[1]->Fill(fabs(lepP4.Eta()));
	h_pass[2]->Fill(leadJetP4.Perp());
	h_pass[3]->Fill(subJetP4.Perp());
	h_pass[4]->Fill(lepP4.DeltaR(leadJetP4));
	h_pass[5]->Fill(lepP4.DeltaR(subJetP4));
      }
    }

    tree->Delete();

    for (int jj = 0; jj < 6; jj++){
      h_pass[jj]->Sumw2();
      h_tot[jj]->Sumw2();
      h_eff[ii][jj] = (TH1F*) h_pass[jj]->Clone();
      h_eff[ii][jj]->SetName(vars[jj]);
      h_eff[ii][jj]->GetYaxis()->SetTitle("Efficiency");
      h_eff[ii][jj]->GetXaxis()->SetTitle(titles[jj]);
      h_eff[ii][jj]->Divide(h_pass[jj], h_tot[jj], 1.0, 1.0, "B");

      h_pass[jj]->Delete();
      h_tot[jj]->Delete();
    }
  }

  for (int jj = 0; jj < 6; jj++){
    
    TH1F* h_SF = (TH1F*) h_eff[1][jj]->Clone();
    h_SF->Sumw2();
    h_SF->Divide(h_eff[0][jj]);

    h_eff[0][jj]->GetYaxis()->SetRangeUser(0.5,1.2);
    h_eff[0][jj]->SetLineColor(kBlue);
    h_SF->SetLineColor(kRed);
    
    TCanvas* c = new TCanvas("c","c",900,800);

    h_eff[0][jj]->Draw("");
    h_eff[1][jj]->Draw("same");
    h_SF->Draw("same");
    
    TLegend* leg;
    leg = new TLegend(0.3,0.2,0.6,0.4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_eff[0][jj], "MC efficiency", "l");
    leg->AddEntry(h_eff[1][jj], "Data efficiency", "l");
    leg->AddEntry(h_SF, "Data / MC scale factor", "l");
    leg->Draw();
    
    c->SaveAs("Plots/eff_"+vars[jj]+"_"+channel+".pdf");

    if (jj == 0){
      TFile* f_SF = new TFile(channel+"TrigEff.root","recreate");
      h_SF->Write();
      f_SF->Close();
      delete f_SF;
    }

    h_SF->Delete();
  }
        
}
