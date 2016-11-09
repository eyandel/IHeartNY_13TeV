#ifndef makePlots_h
#define makePlots_h

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <THStack.h>
#include <TColor.h>
#include <TFile.h>
#include <TROOT.h>
#include <Rtypes.h>
#include <vector>
#include <iostream>


// -------------------------------------------------------------------------------------
// various declarations
// -------------------------------------------------------------------------------------

//const double LUM[2] = {12358.75,12295.65}; //pb-1
// Missing small part of mu, el datasets -- rescale lumi by reduced fraction of events
const double LUM[2] = {12337.98,12267.67}; //pb-1

// -------------------------------------------------------------------------------------
// helper class for summed, weighted histograms (e.g. single top)
// -------------------------------------------------------------------------------------

class SummedHist {
 public : 
 SummedHist( TString const & name, int color ) : name_(name), color_(color) {
    summedHist_ = 0;
  };
  
  // return the summed histogram if created, else create it (summing up histograms in vector hists_) 
  TH1F* hist() { 
    if (summedHist_ != 0) {
      return summedHist_; 
    }
    else if (hists_.size() == 0) {
      return 0; 
    } 
    else {
      summedHist_ = (TH1F*)hists_[0]->Clone();
      summedHist_->SetName( name_ );
      summedHist_->SetFillColor( color_ );
      for (unsigned int j = 1; j<hists_.size(); ++j) {
	summedHist_->Add( hists_[j], 1.0 );
      }
      return summedHist_; 
    };
  }
  
  // return the vector of input histograms  
  std::vector<TH1F*> const & hists() const {
    return hists_;
  }
  
  // add histogram to the vector of inputs
  void push_back( TH1F const * ihist, double norm ) {
    TH1F* clone = (TH1F*) ihist->Clone();
    TString iname( name_ );
    iname += hists_.size();
    clone->SetName( iname );
    clone->Scale( norm );
    hists_.push_back( clone );
    norms_.push_back( norm );
  };
  

 protected : 
  
  std::vector<TH1F*> hists_;
  std::vector<double> norms_;
  
  TString name_; 
  int color_;
  TH1F* summedHist_; 
  
};

TH1* getHist(TString filename, TString histname, TString region){
  TH1::AddDirectory(kFALSE); 

  TFile* infile = TFile::Open( filename );
  TH1* hist;
  if (region == "1t1b" || region == "1t" || region == "1b" || region == "Pre" || region == "") {
    hist = (TH1*) infile->Get(histname+region);
    hist->Sumw2();
  }
  if (region == "1t0b"){
    hist = (TH1*) infile->Get(histname+"1t");
    hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b");
    h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t1b"){
    hist = (TH1*) infile->Get(histname+"1b");
    hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b");
    h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t"){
    hist = (TH1*) infile->Get(histname+"Pre");
    hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t");
    h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0b"){
    hist = (TH1*) infile->Get(histname+"Pre");
    hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1b");
    h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t0b"){
    hist = (TH1*) infile->Get(histname+"Pre");
    hist->Sumw2();
    TH1* h_tmp1 = (TH1*) infile->Get(histname+"1t");
    TH1* h_tmp2 = (TH1*) infile->Get(histname+"1b");
    TH1* h_tmp3 = (TH1*) infile->Get(histname+"1t1b");
    h_tmp1->Sumw2();
    h_tmp2->Sumw2();
    h_tmp3->Sumw2();
    hist->Add(h_tmp3);
    hist->Add(h_tmp2,-1.0);
    hist->Add(h_tmp1,-1.0);
    delete h_tmp1;
    delete h_tmp2;
    delete h_tmp3;
  }
  delete infile;
  return hist;
  
}

// -------------------------------------------------------------------------------------
// W+jets
// -------------------------------------------------------------------------------------

SummedHist * getWJets( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  const int nwjets = 7;
  
  TString wjets_names[nwjets] = {
    "WJets_HT100to200",
    "WJets_HT200to400",
    "WJets_HT400to600",
    "WJets_HT600to800",
    "WJets_HT800to1200",
    "WJets_HT1200to2500",
    "WJets_HT2500toInf",
  };
  
  double wjets_norms[nwjets] = {
    1345.0 * 1.21 * LUM[ich] / 27529599.,  // Note: these are from AN-15-107, may need updating
    359.7 * 1.21 * LUM[ich] / 4963240.,  
    48.91 * 1.21 * LUM[ich] / 1963464.,  
    12.05 * 1.21 * LUM[ich] / 3722395.,
    5.501 * 1.21 * LUM[ich] / 6314257.,
    1.329 * 1.21 * LUM[ich] / 6768156.,
    0.03216 * 1.21 * LUM[ich] / 253561.,
  };
  
  SummedHist* wjets = new SummedHist( histname, kGreen-3 );
  
  for (int i=0 ; i<nwjets; i++) {
    TString iname = DIR + "hists_" + wjets_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1F* hist = (TH1F*) getHist(iname,histname,region);
    if (hist->Integral() > 0.0) wjets->push_back( hist, wjets_norms[i] );
  }
  
  return wjets;
  
}

// -------------------------------------------------------------------------------------
// single top
// -------------------------------------------------------------------------------------

SummedHist * getSingleTop( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nsingletop = 4;
  
  TString singletop_names[nsingletop] = {
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
  };
  
  double singletop_norms[nsingletop] = {
    136.02*0.322 * LUM[ich] / 3279200., // xsec from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    80.95*0.322 * LUM[ich] / 1682400.,  // BR from https://twiki.cern.ch/twiki/bin/view/Main/EdbrBackup (second-hand, but can't find original source)
    35.9 * LUM[ich] / 998400.,          // BR are needed because t-channel samples are leptonic final state only
    35.9 * LUM[ich] / 985000.,
  };
  
  SummedHist* singletop = new SummedHist( histname, 6 );
  
  for (int i=0; i<nsingletop; i++) {
    TString iname = DIR + "hists_" + singletop_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1F* hist = (TH1F*) getHist(iname,histname,region);
    if (hist->Integral() > 0.0) singletop->push_back( hist, singletop_norms[i] );
  }
  
  return singletop;
  
}


// -------------------------------------------------------------------------------------
// non-semileptonic ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbarNonSemiLep( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost ) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  TString ttbar_name = "PowhegPythia8_nonsemilep";
  double ttbar_norm = 831.76 * LUM[ich] / 182075200.; //182123200 - 48000
  
  SummedHist* ttbar = new SummedHist( histname, kRed-7);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1F* hist = (TH1F*) getHist(iname,histname,region);
  if (hist->Integral() > 0.0) ttbar->push_back( hist, ttbar_norm );
  
  return ttbar;
  
}


// -------------------------------------------------------------------------------------
// signal ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbar( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost) {
  
  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  TString ttbar_name = "PowhegPythia8_fullTruth";
  double ttbar_norm = 831.76 * LUM[ich] / 182123200.;
  
  SummedHist* ttbar = new SummedHist( histname, kRed+1);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1F* hist = (TH1F*) getHist(iname,histname,region);
  if (hist->Integral() > 0.0) ttbar->push_back( hist, ttbar_norm );
  
  return ttbar;
  
}

// -------------------------------------------------------------------------------------
// QCD
// -------------------------------------------------------------------------------------

SummedHist * getQCDMC( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nqcd = 6;
  
  TString qcd_names[nqcd] = {
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
  };
  
  double qcd_norms[nqcd] = {
    347700. * LUM[ich] / 37828442., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM[ich] / 44058594.,  // Hopefully they are close enough...
    6831. * LUM[ich] / 29832311.,  
    1207. * LUM[ich] / 4963881., //4980387 - 16506 (missing file)
    119.9 * LUM[ich] / 7803965.,
    25.24 * LUM[ich] / 4047532.,
  };
  
  SummedHist* qcd = new SummedHist( histname, kYellow );
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1F* hist = (TH1F*) getHist(iname,histname,region);
    if (hist->Integral() > 0.0) qcd->push_back( hist, qcd_norms[i] );
  }
  
  return qcd;
  
}

SummedHist * getData( TString DIR, TString histname, TString region, TString channel, bool isQCD) {
  
  TString append = "";
  if (isQCD) append = "_qcd";
    
  SummedHist* data = new SummedHist( histname, 0 );
  
  TString iname = DIR + "hists_Data_" + channel + append + ".root";
  TH1F* hist = (TH1F*) getHist(iname,histname,region);
  if (hist->Integral() > 0.0) data->push_back( hist, 1.0 );
  
  return data;
  
}

float getQCDnorm( TString DIR, TString channel, TString region, TString syst, bool usePost) {

  TString append = "";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nqcd = 6;
  
  TString qcd_names[nqcd] = {
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf"
  };

  double qcd_norms[nqcd] = {
    347700. * LUM[ich] / 37828442., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM[ich] / 44058594.,  // Hopefully they are close enough...
    6831. * LUM[ich] / 29832311.,  
    1207. * LUM[ich] / 4980387.,
    119.9 * LUM[ich] / 7803965.,
    25.24 * LUM[ich] / 4047532.,
  };

  float n_qcd = 0.0;
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1F* hist = (TH1F*) getHist(iname,"lepPhi",region);
    if (hist->Integral() > 0.0) hist->Scale(qcd_norms[i]);
    n_qcd += hist->Integral();
  }
  
  return n_qcd;
  
}

TH1F * getQCDData(TString sigDIR, TString sideDIR, TString histname, TString region, TString channel, TString syst, bool usePost) {

  SummedHist* wjets = getWJets( sideDIR, histname, region, channel, true, syst, usePost);
  SummedHist* singletop = getSingleTop( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* ttbar = getTTbar( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* data = getData( sideDIR, histname, region, channel, true);
  
  TH1F* h_wjets = (TH1F*) wjets->hist();
  TH1F* h_singletop = (TH1F*) singletop->hist();
  TH1F* h_ttbar = (TH1F*) ttbar->hist();
  TH1F* h_ttbar_nonSemiLep = (TH1F*) ttbar_nonSemiLep->hist();
  TH1F* h_data = (TH1F*) data->hist();

  if (h_data) {
    TH1F* h_qcd = (TH1F*) h_data->Clone("QCD");
    if (h_wjets) h_qcd->Add(h_wjets,-1.0);
    if (h_singletop) h_qcd->Add(h_singletop,-1.0);
    if (h_ttbar) h_qcd->Add(h_ttbar,-1.0);
    if (h_ttbar_nonSemiLep) h_qcd->Add(h_ttbar_nonSemiLep,-1.0);

    for (int ii = 1; ii < h_qcd->GetNbinsX()+1; ii++){
      if (h_qcd->GetBinContent(ii) < 0.0) h_qcd->SetBinContent(ii,0.0);
    }
    
    float n_qcd = getQCDnorm(sigDIR,channel,region,syst,usePost);
    h_qcd->Scale(n_qcd / h_qcd->Integral());
    h_qcd->SetFillColor(kYellow);
    
    //wjets->Delete();
    //singletop->Delete();
    //ttbar->Delete();
    //ttbar_nonSemiLep->Delete();
    //data->Delete();
    
    return h_qcd;
  }
  else return 0;
  
}

TObject * getSignal( TString DIR, TString histname, TString region, TString channel, bool isQCD ) {

  TString append = "";
  if (isQCD) append = "_qcd";

  int ich = 0;
  if (channel == "el") ich = 1;

  TString ttbar_name = "PowhegPythia8_semilep";
  double ttbar_norm = 831.76 * LUM[ich] / 182123200.;
    
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + append + "_nom.root";
  TH2F* ttbar = (TH2F*) getHist(iname,histname,region);
  if (ttbar->Integral() > 0.0) ttbar->Scale(ttbar_norm);
  
  return ttbar;
  
}

TObject * getBackground( TString DIR, TString histname, TString region, TString channel, bool isQCD ) {

  TString append = "";
  if (isQCD) append = "_qcd";

  int ich = 0;
  if (channel == "el") ich = 1;
  
  const int nbkg = 18;
  
  TString bkg_names[nbkg] = {
    "PowhegPythia8_nonsemilep",
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
    "WJets_HT100to200",
    "WJets_HT200to400",
    "WJets_HT400to600",
    "WJets_HT600to800",
    "WJets_HT800to1200",
    "WJets_HT1200to2500",
    "WJets_HT2500toInf",
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf"
  };
  
  double bkg_norms[nbkg] = {
    831.76 * LUM[ich] / 182123200.,       //TTbar nonsignal
    136.02*0.322 * LUM[ich] / 3279200.,   //SingleTop
    80.95*0.322 * LUM[ich] / 1682400.,  
    35.9 * LUM[ich] / 998400.,          
    35.9 * LUM[ich] / 985000.,
    1345.0 * 1.21 * LUM[ich] / 27529599., //WJets
    359.7 * 1.21 * LUM[ich] / 4963240.,  
    48.91 * 1.21 * LUM[ich] / 1963464.,  
    12.05 * 1.21 * LUM[ich] / 3722395.,
    5.501 * 1.21 * LUM[ich] / 6314257.,
    1.329 * 1.21 * LUM[ich] / 6768156.,
    0.03216 * 1.21 * LUM[ich] / 253561.,
    347700. * LUM[ich] / 37828442.,     //QCD
    32100. * LUM[ich] / 44058594.,  
    6831. * LUM[ich] / 29832311.,  
    1207. * LUM[ich] / 4980387.,
    119.9 * LUM[ich] / 7803965.,
    25.24 * LUM[ich] / 4047532.
  };
  
  TH2F* bkg_hists[nbkg];
  
  for (int i=0; i<nbkg; i++) {
    TString iname = DIR + "hists_" + bkg_names[i] + "_" + channel + "_nom" + append + ".root";
    bkg_hists[i] = (TH2F*) getHist(iname,histname,region);
    if (bkg_hists[i]->Integral() > 0.0) bkg_hists[i]->Scale(bkg_norms[i]);
  }
  
  TH2F* bkg = (TH2F*) bkg_hists[0]->Clone();
  for (int j=1; j<nbkg; j++){
    bkg->Add(bkg_hists[j]);
  }

  return bkg;
  
}


#endif
