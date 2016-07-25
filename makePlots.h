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

const double LUM = 2689; //pb-1

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

// -------------------------------------------------------------------------------------
// W+jets
// -------------------------------------------------------------------------------------

SummedHist * getWJets( TString DIR, TString histname, TString channel, bool isQCD, TString syst) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
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
    1345.0 * 1.21 * LUM / 10205377.,  // Note: these are from AN-15-107, may need updating
    359.7 * 1.21 * LUM / 4949568.,  
    48.91 * 1.21 * LUM / 1943664.,  
    12.05 * 1.21 * LUM / 3767766.,
    5.501 * 1.21 * LUM / 1568277.,
    1.329 * 1.21 * LUM / 246239.,
    0.03216 * 1.21 * LUM / 251982.,
  };
  
  SummedHist* wjets = new SummedHist( histname, kGreen-3 );
  
  for (int i=0 ; i<nwjets; i++) {
    TString iname = DIR + "hists_" + wjets_names[i] + "_" + channel + "_" + syst + append + ".root";
    TFile* infile = TFile::Open( iname );
    TH1F* hist = (TH1F*) infile->Get(histname);
    if (hist->Integral() > 0.0){
      hist->Sumw2();
      wjets->push_back( hist, wjets_norms[i] );
    }
    delete infile;
  }
  
  return wjets;
  
}

// -------------------------------------------------------------------------------------
// single top
// -------------------------------------------------------------------------------------

SummedHist * getSingleTop( TString DIR, TString histname, TString channel, bool isQCD, TString syst) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
  const int nsingletop = 5;
  
  TString singletop_names[nsingletop] = {
    "SingleTop_t_s",
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
  };
  
  double singletop_norms[nsingletop] = {
    3.36 * LUM / 998400.,     // Use xsec from AN-15-107 since this sample was leptonDecays only (76X B2G samples)
    136.02 * LUM / 64925700., // For all following, use xsec from 
    80.95 * LUM / 38932192.,  // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    35.9 * LUM / 1000000.,
    35.9 * LUM / 999400.,
  };
  
  SummedHist* singletop = new SummedHist( histname, 6 );
  
  for (int i=0; i<nsingletop; i++) {
    TString iname = DIR + "hists_" + singletop_names[i] + "_" + channel + "_" + syst + append + ".root";
    TFile* infile = TFile::Open( iname );
    TH1F* hist = (TH1F*) infile->Get(histname);
    if (hist->Integral() > 0.0){
      hist->Sumw2();
      singletop->push_back( hist, singletop_norms[i] );
    }
    delete infile;
  }
  
  return singletop;
  
}


// -------------------------------------------------------------------------------------
// non-semileptonic ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbarNonSemiLep( TString DIR, TString histname, TString channel, bool isQCD, TString syst ) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
  TString ttbar_name = "PowhegPythia8_nonsemilep";
  double ttbar_norm = 831.76 * LUM / 187626200.;
  
  SummedHist* ttbar = new SummedHist( histname, kRed-7);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TFile* infile = TFile::Open( iname );
  TH1F* hist = (TH1F*) infile->Get(histname);
  if (hist->Integral() > 0.0){
    hist->Sumw2();
    ttbar->push_back( hist, ttbar_norm );
  }
  delete infile;
  
  return ttbar;
  
}


// -------------------------------------------------------------------------------------
// signal ttbar
// -------------------------------------------------------------------------------------

SummedHist * getTTbar( TString DIR, TString histname, TString channel, bool isQCD, TString syst) {
  
  TString append = "";
  if (isQCD) append = "_qcd";
  
  TString ttbar_name = "PowhegPythia8_semilep";
  double ttbar_norm = 831.76 * LUM / 187626200.;
  
  SummedHist* ttbar = new SummedHist( histname, kRed+1);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TFile* infile = TFile::Open( iname );
  TH1F* hist = (TH1F*) infile->Get(histname);
  if (hist->Integral() > 0.0){
    hist->Sumw2();
    ttbar->push_back( hist, ttbar_norm );
  }
  delete infile;
  
  return ttbar;
  
}

// -------------------------------------------------------------------------------------
// QCD
// -------------------------------------------------------------------------------------

SummedHist * getQCDMC( TString DIR, TString histname, TString channel, bool isQCD, TString syst) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
  //const int nqcd = 6;
  const int nqcd = 5;
  
  TString qcd_names[nqcd] = {
    //"QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
  };
  
  double qcd_norms[nqcd] = {
    //347700. * LUM / 16909004., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM / 19665695.,  // Hopefully they are close enough...
    6831. * LUM / 15547962.,  
    1207. * LUM / 5049267.,
    119.9 * LUM / 3939077.,
    25.24 * LUM / 1981228.,
  };
  
  SummedHist* qcd = new SummedHist( histname, kYellow );
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    TFile* infile = TFile::Open( iname );
    TH1F* hist = (TH1F*) infile->Get(histname);
    if (hist->Integral() > 0.0){
      hist->Sumw2();
      qcd->push_back( hist, qcd_norms[i] );
    }
    delete infile;
  }
  
  return qcd;
  
}

SummedHist * getData( TString DIR, TString histname, TString channel, bool isQCD) { //TODO: decide whether to get data like this or more simply in makePlots.cc
  
  TString append = "";
  if (isQCD) append = "_qcd";
  
  const int ndata = 2;
  
  TString data_names[ndata] = {
    "Data_2015C",
    "Data_2015D",
  };
  
  SummedHist* data = new SummedHist( histname, 0 );
  
  for (int i=0; i<ndata; i++) {
    TString iname = DIR + "hists_" + data_names[i] + "_" + channel + append + ".root";
    TFile* infile = TFile::Open( iname );
    TH1F* hist = (TH1F*) infile->Get(histname);
    if (hist->Integral() > 0.0){
      hist->Sumw2();
      data->push_back( hist, 1.0 );
    }
    delete infile;
  }
  
  return data;
  
}

float getQCDnorm( TString DIR, TString channel, TString region, TString syst) {
  
  //const int nqcd = 6;
  const int nqcd = 5;
  
  TString qcd_names[nqcd] = {
    //"QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf"
  };
  
  double qcd_norms[nqcd] = {
    //347700. * LUM / 16909004., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM / 19665695., // Hopefully they are close enough...
    6831. * LUM / 15547962.,  
    1207. * LUM / 5049267.,
    119.9 * LUM / 3939077.,
    25.24 * LUM / 1981228.
  };

  float n_qcd = 0.0;
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + ".root";
    TFile* infile = TFile::Open( iname );
    TH1F* hist = (TH1F*) infile->Get("metPt"+region);
    if (hist->Integral() > 0.0){
      hist->Sumw2();
      hist->Scale(qcd_norms[i]);
    }
    n_qcd += hist->Integral();
    delete infile;
  }
  
  return n_qcd;
  
}

TH1F * getQCDData(TString sigDIR, TString sideDIR, TString histname, TString channel, TString syst) {

  SummedHist* wjets = getWJets( sideDIR, histname, channel, true, syst );
  SummedHist* singletop = getSingleTop( sideDIR, histname, channel, true, syst );
  SummedHist* ttbar = getTTbar( sideDIR, histname, channel, true, syst );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( sideDIR, histname, channel, true, syst );
  SummedHist* data = getData( sideDIR, histname, channel, true);
  
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
    
    float n_qcd = 0.0;
    if (histname.Contains("1t1b")) n_qcd = getQCDnorm(sigDIR,channel,"1t1b",syst);
    else if (histname.Contains("1t0b")) n_qcd = getQCDnorm(sigDIR,channel,"1t0b",syst);
    else if (histname.Contains("0t0b")) n_qcd = getQCDnorm(sigDIR,channel,"0t0b",syst);
    else if (histname.Contains("0t1b")) n_qcd = getQCDnorm(sigDIR,channel,"0t1b",syst);
    else if (histname.Contains("0t")) n_qcd = getQCDnorm(sigDIR,channel,"0t",syst);
    else n_qcd = getQCDnorm(sigDIR,channel,"Pre",syst);
    
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

TObject * getSignal( TString DIR, TString histname, TString channel, bool isQCD ) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
  TString ttbar_name = "PowhegPythia8_semilep";
  double ttbar_norm = 831.76 * LUM / 187626200.;
    
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + append + "_nom.root";
  TFile* infile = TFile::Open( iname );
  TH2F* ttbar = (TH2F*) infile->Get(histname);
  if (ttbar->Integral() > 0.0){
    ttbar->Sumw2();
    ttbar->Scale(ttbar_norm);
  }
  delete infile;
  
  return ttbar;
  
}

TObject * getBackground( TString DIR, TString histname, TString channel, bool isQCD ) {

  TString append = "";
  if (isQCD) append = "_qcd";
  
  //const int nbkg = 19;
  const int nbkg = 18;
  
  TString bkg_names[nbkg] = {
    "PowhegPythia8_nonsemilep",
    "SingleTop_t_s",
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
    //"QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf"
  };
  
  double bkg_norms[nbkg] = {
    831.76 * LUM / 187626200.,       //TTbar nonsignal
    3.36 * LUM / 998400.,            // Single Top
    136.02 * LUM / 64925700., 
    80.95 * LUM / 38932192.,  
    35.9 * LUM / 1000000.,
    35.9 * LUM / 999400.,
    1345.0 * 1.21 * LUM / 10205377., // WJets
    359.7 * 1.21 * LUM / 4949568.,  
    48.91 * 1.21 * LUM / 1943664.,  
    12.05 * 1.21 * LUM / 3767766.,
    5.501 * 1.21 * LUM / 1568277.,
    1.329 * 1.21 * LUM / 246239.,
    0.03216 * 1.21 * LUM / 251982.,
    //347700. * LUM / 16909004.,       // QCD
    32100. * LUM / 19665695., 
    6831. * LUM / 15547962.,  
    1207. * LUM / 5049267.,
    119.9 * LUM / 3939077.,
    25.24 * LUM / 1981228.
  };
  
  TH2F* bkg_hists[nbkg];
  
  for (int i=0; i<nbkg; i++) {
    TString iname = DIR + "hists_" + bkg_names[i] + "_" + channel + "_nom" + append + ".root";
    TFile* infile = TFile::Open( iname );
    bkg_hists[i] = (TH2F*) infile->Get(histname);
    if (bkg_hists[i]->Integral() > 0.0){
      bkg_hists[i]->Sumw2();
      bkg_hists[i]->Scale(bkg_norms[i]);
    }

    delete infile;
  }
  
  TH2F* bkg = (TH2F*) bkg_hists[0]->Clone();
  for (int j=1; j<nbkg; j++){
    bkg->Add(bkg_hists[j]);
  }

  return bkg;
  
}


#endif
