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

const double LUM[2] = {35867.0,35867.0}; //pb-1
// TODO: use correct lumis for B2G datasets (rather than lumi for total golden 2016 dataset)

// -------------------------------------------------------------------------------------
// helper class for summed, weighted histograms (e.g. single top)
// -------------------------------------------------------------------------------------

class SummedHist {
 public : 
 SummedHist( TString const & name, int color ) : name_(name), color_(color) {
    summedHist_ = 0;
  };
  
  // return the summed histogram if created, else create it (summing up histograms in vector hists_) 
  TH1* hist(bool noFill = false) { 
    if (summedHist_ != 0) {
      if (noFill) summedHist_->SetFillColor(0);
      return summedHist_; 
    }
    else if (hists_.size() == 0) {
      return 0; 
    } 
    else {
      summedHist_ = (TH1*)hists_[0]->Clone();
      summedHist_->SetName( name_ );
      if (!noFill) summedHist_->SetFillColor( color_ );
      for (unsigned int j = 1; j<hists_.size(); ++j) {
	summedHist_->Add( hists_[j], 1.0 );
      }
      return summedHist_; 
    };
  }
  
  // return the vector of input histograms  
  std::vector<TH1*> const & hists() const {
    return hists_;
  }
  
  // add histogram to the vector of inputs
  void push_back( TH1 const * ihist, double norm ) {
    TH1* clone = (TH1*) ihist->Clone();
    TString iname( name_ );
    iname += hists_.size();
    clone->SetName( iname );
    clone->Scale( norm );
    hists_.push_back( clone );
    norms_.push_back( norm );
  };
  

 protected : 
  
  std::vector<TH1*> hists_;
  std::vector<double> norms_;
  
  TString name_; 
  int color_;
  TH1* summedHist_; 
  
};

TH1* getHist(TString filename, TString histname, TString region, TString split = ""){
  TH1::AddDirectory(kFALSE); 

  TString append = "";
  if (split != "") append = "_"+split;
  
  TFile* infile = TFile::Open( filename );
  TH1* hist;
  if (region == "1t1b" || region == "1t" || region == "1b" || region == "Pre" || region == "") {
    hist = (TH1*) infile->Get(histname+region+append);
    //hist->Sumw2();
  }
  if (region == "1t0b"){
    hist = (TH1*) infile->Get(histname+"1t"+append);
    //hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b"+append);
    //h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t1b"){
    hist = (TH1*) infile->Get(histname+"1b"+append);
    //hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t1b"+append);
    //h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t"){
    hist = (TH1*) infile->Get(histname+"Pre"+append);
    //hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1t"+append);
    //h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0b"){
    hist = (TH1*) infile->Get(histname+"Pre"+append);
    //hist->Sumw2();
    TH1* h_tmp = (TH1*) infile->Get(histname+"1b"+append);
    //h_tmp->Sumw2();
    hist->Add(h_tmp,-1.0);
    delete h_tmp;
  }
  if (region == "0t0b"){
    hist = (TH1*) infile->Get(histname+"Pre"+append);
    //hist->Sumw2();
    TH1* h_tmp1 = (TH1*) infile->Get(histname+"1t"+append);
    TH1* h_tmp2 = (TH1*) infile->Get(histname+"1b"+append);
    TH1* h_tmp3 = (TH1*) infile->Get(histname+"1t1b"+append);
    //h_tmp1->Sumw2();
    //h_tmp2->Sumw2();
    //h_tmp3->Sumw2();
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

SummedHist * getWJets( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, TString split = "") {

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
    1345.0 * 1.21 * LUM[ich] / 39617787.,  // Note: cross sections are from AN-15-107, may need updating
    359.7 * 1.21 * LUM[ich] / 19914590.,  
    48.91 * 1.21 * LUM[ich] / 5796237.,  
    12.05 * 1.21 * LUM[ich] / 14822888.,
    5.501 * 1.21 * LUM[ich] / 6200954.,
    1.329 * 1.21 * LUM[ich] / 6324934.,
    0.03216 * 1.21 * LUM[ich] / 2384260.,
  };

  int plotcolor = kGreen-3;
  if (split == "bb") plotcolor = kRed+1;
  if (split == "b") plotcolor = kRed-7;
  if (split == "cc") plotcolor = 6;
  if (split == "c") plotcolor = kGreen-3;
  if (split == "l") plotcolor = kYellow;
  
  SummedHist* wjets = new SummedHist( histname, plotcolor );
  
  for (int i=0 ; i<nwjets; i++) {
    TString iname = DIR + "hists_" + wjets_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region,split);
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

  const int nsingletop = 5;
  
  TString singletop_names[nsingletop] = {
    "SingleTop_t_t",
    "SingleTop_tbar_t",
    "SingleTop_t_tW",
    "SingleTop_tbar_tW",
    "SingleTop_s",
  };
  
  double singletop_norms[nsingletop] = {
    136.02 * LUM[ich] / 5993676., // xsec from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    80.95 * LUM[ich] / 3928063.,  // BR from https://twiki.cern.ch/twiki/bin/view/Main/EdbrBackup (second-hand, but can't find original source)
    35.9 * LUM[ich] / 992024.,    // BR is needed because s-channel sample is leptonic final state only
    35.9 * LUM[ich] / 998276.,
    10.32 * 0.322 * LUM[ich] / 1000000.,
  };
  
  SummedHist* singletop = new SummedHist( histname, 6 );
  
  for (int i=0; i<nsingletop; i++) {
    TString iname = DIR + "hists_" + singletop_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,histname,region);
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
  double ttbar_norm = 831.76 * LUM[ich] / 77229341.; //TODO: technically this is incorrect, as there is a job missing -- but can't tell what files the job corresponds to
  
  SummedHist* ttbar = new SummedHist( histname, kRed-7);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region);
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
  double ttbar_norm = 831.76 * LUM[ich] / 77229341.;
  
  SummedHist* ttbar = new SummedHist( histname, kRed+1);
  TString iname = DIR + "hists_" + ttbar_name + "_" + channel + "_" + syst + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region);
  if (hist->Integral() > 0.0) ttbar->push_back( hist, ttbar_norm );
  
  return ttbar;
  
}

// -------------------------------------------------------------------------------------
// QCD
// -------------------------------------------------------------------------------------

SummedHist * getQCDMC( TString DIR, TString histname, TString region, TString channel, bool isQCD, TString syst, bool usePost, bool elID = false) {

  TString append = "";
  if (isQCD) append = "_qcd";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

  const int nqcd = 5;
  
  TString qcd_names[nqcd] = {
    //"QCD_HT300to500", //Empty currently
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf",
  };
  
  double qcd_norms[nqcd] = {
    //347700. * LUM[ich] / 17035891., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM[ich] / 18929951.,  // Technically part of this is missing, but can't tell which part -- ~1% though
    6831. * LUM[ich] / 15629253.,  
    1207. * LUM[ich] / 4767100., 
    119.9 * LUM[ich] / 3970819.,
    25.24 * LUM[ich] / 1991645.,
  };
  
  SummedHist* qcd = new SummedHist( histname, kYellow );
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    if (elID) iname = DIR + "hists_" + qcd_names[i] + "_elID.root";
    TH1* hist = (TH1*) getHist(iname,histname,region);
    if (hist->Integral() > 0.0) qcd->push_back( hist, qcd_norms[i] );
  }
  
  return qcd;
  
}

SummedHist * getData( TString DIR, TString histname, TString region, TString channel, bool isQCD) {
  
  TString append = "";
  if (isQCD) append = "_qcd";
    
  SummedHist* data = new SummedHist( histname, 0 );
  
  TString iname = DIR + "hists_Data_" + channel + append + ".root";
  TH1* hist = (TH1*) getHist(iname,histname,region);
  if (hist->Integral() > 0.0) data->push_back( hist, 1.0 );
  
  return data;
  
}

float getQCDnorm( TString DIR, TString channel, TString region, TString syst, bool usePost) {

  TString append = "";
  if (usePost) append += "_post";

  int ich = 0;
  if (channel == "el") ich = 1;

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
    //347700. * LUM[ich] / 17035891., // Cross sections are from AN-15-136, which is a bit random and not actually correct for the samples I use
    32100. * LUM[ich] / 18929951.,  // Hopefully they are close enough...
    6831. * LUM[ich] / 15629253.,  
    1207. * LUM[ich] / 4767100.,
    119.9 * LUM[ich] / 3970819.,
    25.24 * LUM[ich] / 1991645.,
  };

  float n_qcd = 0.0;
  
  for (int i=0; i<nqcd; i++) {
    TString iname = DIR + "hists_" + qcd_names[i] + "_" + channel + "_" + syst + append + ".root";
    TH1* hist = (TH1*) getHist(iname,"lepPhi",region);
    if (hist->Integral() > 0.0) hist->Scale(qcd_norms[i]);
    n_qcd += hist->Integral();
  }
  
  return n_qcd;
  
}

TH1 * getQCDData(TString sigDIR, TString sideDIR, TString histname, TString region, TString channel, TString syst, bool usePost) {

  SummedHist* wjets = getWJets( sideDIR, histname, region, channel, true, syst, usePost);
  SummedHist* singletop = getSingleTop( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* ttbar = getTTbar( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* ttbar_nonSemiLep = getTTbarNonSemiLep( sideDIR, histname, region, channel, true, syst, usePost );
  SummedHist* data = getData( sideDIR, histname, region, channel, true);
  
  TH1* h_wjets = (TH1*) wjets->hist();
  TH1* h_singletop = (TH1*) singletop->hist();
  TH1* h_ttbar = (TH1*) ttbar->hist();
  TH1* h_ttbar_nonSemiLep = (TH1*) ttbar_nonSemiLep->hist();
  TH1* h_data = (TH1*) data->hist();

  if (h_data) {
    TH1* h_qcd = (TH1*) h_data->Clone("QCD");
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
    "SingleTop_s",
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
    831.76 * LUM[ich] / 77229341.,       //TTbar nonsignal
    136.02 * LUM[ich] / 5993676.,        //Single top
    80.95 * LUM[ich] / 3928063.,  
    35.9 * LUM[ich] / 992024.,    
    35.9 * LUM[ich] / 998276.,
    10.32 * 0.322 * LUM[ich] / 1000000.,
    1345.0 * 1.21 * LUM[ich] / 39617787., // WJets
    359.7 * 1.21 * LUM[ich] / 19914590.,  
    48.91 * 1.21 * LUM[ich] / 5796237.,  
    12.05 * 1.21 * LUM[ich] / 14822888.,
    5.501 * 1.21 * LUM[ich] / 6200954.,
    1.329 * 1.21 * LUM[ich] / 6324934.,
    0.03216 * 1.21 * LUM[ich] / 2384260.,
    //347700. * LUM[ich] / 17035891.,       // QCD
    32100. * LUM[ich] / 18929951.,  
    6831. * LUM[ich] / 15629253.,  
    1207. * LUM[ich] / 4767100., 
    119.9 * LUM[ich] / 3970819.,
    25.24 * LUM[ich] / 1991645.,
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
