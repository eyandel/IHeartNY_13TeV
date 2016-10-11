// Script to run makeHists.C
// Syntax is makeHists(TString INDIR, TString OUTDIR, TString sample, TString channel, bool isData = false, bool isSignal = false, TString lepID = "Medium", TString iso = "None", bool doHemiCuts = false, float metCut = 0.0, bool doTriangular = false, bool isQCD = false, TString systematic = "nom", int oddOrEven = 0, bool usePost = false)

#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void runMakeHists(TString toMake = "prefit"){

  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L BTagCalibrationStandalone.cc++");
  gROOT->ProcessLine(".include RooUnfold/src");
  gSystem->Load("RooUnfold/libRooUnfold");
  gROOT->ProcessLine(".L makeHists.C++");

  const int nBKG = 18;
  TString bkgMCnames[nBKG] = {
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

  const int nSYS = 13;
  TString sysnames[nSYS] = {"nom","puUp","puDown","JECUp","JECDown","JERUp","JERDown","lepUp","lepDown","BTagUp","BTagDown","TopTagUp","TopTagDown"};
  const int nTHSYS = 6;
  TString thsysnames[nTHSYS] = {"PDFUp","PDFDown","Q2Up","Q2Down","ASUp","ASDown"};
  const int nISO = 7;
  TString isoWPs[nISO] = {"MiniIso10","MiniIso20","2DisoPt25","2DisoPt45","2DisoB2G","2DisoIHNY","Loose"};

  // ----------------------------------------------
  // Make histfiles for lepton optimization studies

  if (toMake == "all" || toMake == "lepOpt"){

    for (int ii = 0; ii < nISO; ii++){
      // Run signal
      makeHists("skimTrees_80X/","histfiles_80X_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_fullTruth","mu",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
      makeHists("skimTrees_80X/","histfiles_80X_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_fullTruth","el",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
      
      //Run QCD
      for (int jj = nBKG-6; jj < nBKG; jj ++){
	makeHists("skimTrees_80X/","histfiles_80X_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"mu",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);
	makeHists("skimTrees_80X/","histfiles_80X_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"el",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"nom",0,false);    
      }
    }

    // Run signal
    makeHists("skimTrees_80X/","histfiles_80X_tMu_tEl","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
    makeHists("skimTrees_80X/","histfiles_80X_tMu_tEl","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
    
    //Run QCD
    for (int jj = nBKG-6; jj < nBKG; jj ++){
      makeHists("skimTrees_80X/","histfiles_80X_tMu_tEl",bkgMCnames[jj],"mu",false,false,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);
      makeHists("skimTrees_80X/","histfiles_80X_tMu_tEl",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",false,0.0,false,false,"nom",0,false);    
    }
  }

  // ----------------------------------------------------------------------
  // Make histfiles with final lepton selection, pre-selection-optimization

  if (toMake == "all" || toMake == "selOpt"){
    //run data
    makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10","Data_mu","mu",true,false,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
    makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10","Data_el","el",true,false,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    
    // Run signal
    makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
    makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    
    //run other MCs
    for (int jj = 0; jj < nBKG; jj ++){
      makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,0.0,true,false,"nom",0,false);
      makeHists("skimTrees_80X/","histfiles_80X_mMu_tEl_MiniIso10",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,0.0,true,false,"nom",0,false);
    }
  }
    
  // ----------------------------------------------
  // Make final histfiles

  if (toMake == "all" || toMake == "prefit" || toMake == "postfit"){

    bool postTopTagSF = false;
    if (toMake == "postfit") postTopTagSF = true;
    //run data
    makeHists("skimTrees_80X/","histfiles_80X","Data_mu","mu",true,false,"Tight","MiniIso10",true,35.0,true,false,"nom",0,false); //Data
    makeHists("skimTrees_80X/","histfiles_80X","Data_el","el",true,false,"Tight","MiniIso10",true,35.0,true,false,"nom",0,false);
    
    makeHists("skimTrees_80X/","histfiles_80X","Data_mu","mu",true,false,"Medium","MiniIso10",true,35.0,true,true,"nom",0,false); //QCD
    makeHists("skimTrees_80X/","histfiles_80X","Data_el","el",true,false,"Medium","MiniIso10",true,35.0,true,true,"nom",0,false);
    
    //makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",1,postTopTagSF); //Odd (for unfolding closure)
    //makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",1,postTopTagSF);
    //makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",2,postTopTagSF); //Even (for unfolding closure)
    //makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",2,postTopTagSF);
   
    for (int ii = 0; ii < nSYS; ii++){
      // Run signal
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,sysnames[ii],0,postTopTagSF); //Signal
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,sysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,true,true,sysnames[ii],0,postTopTagSF); //QCD
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Medium","MiniIso10",true,35.0,true,true,sysnames[ii],0,postTopTagSF);
      
      //run other MCs
      for (int jj = 0; jj < nBKG; jj ++){
	makeHists("skimTrees_80X/","histfiles_80X",bkgMCnames[jj],"mu",false,false,"Tight","MiniIso10",true,35.0,true,false,sysnames[ii],0,postTopTagSF); //Signal
	makeHists("skimTrees_80X/","histfiles_80X",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,35.0,true,false,sysnames[ii],0,postTopTagSF);    
	makeHists("skimTrees_80X/","histfiles_80X",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,true,true,sysnames[ii],0,postTopTagSF); //QCD
	makeHists("skimTrees_80X/","histfiles_80X",bkgMCnames[jj],"el",false,false,"Medium","MiniIso10",true,35.0,true,true,sysnames[ii],0,postTopTagSF);    
      }
    }
   
    for (int ii = 0; ii < nTHSYS; ii++){
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,thsysnames[ii],0,postTopTagSF); //Signal
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,thsysnames[ii],0,postTopTagSF);
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Medium","MiniIso10",true,35.0,true,true,thsysnames[ii],0,postTopTagSF); //QCD
      makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Medium","MiniIso10",true,35.0,true,true,thsysnames[ii],0,postTopTagSF);
      }  
  }

  if (toMake == "unfold") {
    
    bool postTopTagSF = true;

    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",1,postTopTagSF); //Odd (for unfolding closure)
    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","mu",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",2,postTopTagSF); //Even (for unfolding closure)

    /*
    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",0,postTopTagSF);
    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",1,postTopTagSF); //Odd (for unfolding closure)
    makeHists("skimTrees_80X/","histfiles_80X","PowhegPythia8_fullTruth","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"nom",2,postTopTagSF); //Even (for unfolding closure)
    */
  }

}
