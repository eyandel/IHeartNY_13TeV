// Script to run makeHists.C
// Syntax is makeHists(TString INDIR, TString OUTDIR, TString sample, TString channel, bool isData = false, bool isSignal = false, TString lepID = "Medium", TString iso = "None", bool doHemiCuts = false, float metCut = 0.0, bool doTriangular = false, bool isQCD = false, TString regions = "fit", TString systematic = "nom")

#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void runMakeHists(TString toMake = "final"){

  //gSystem->CompileMacro("makeHists.C");
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine(".L BTagCalibrationStandalone.cc++");
  gROOT->ProcessLine(".L makeHists.C++");

  TString bkgMCnames[19] = {
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
    "QCD_HT300to500",
    "QCD_HT500to700",
    "QCD_HT700to1000",
    "QCD_HT1000to1500",
    "QCD_HT1500to2000",
    "QCD_HT2000toInf"
  };

  TString sysnames[13] = {"nom","puUp","puDown","JECUp","JECDown","JERUp","JERDown","lepUp","lepDown","BTagUp","BTagDown","TopTagUp","TopTagDown"};
  TString isoWPs[7] = {"MiniIso10","MiniIso20","2DisoPt25","2DisoPt45","2DisoB2G","2DisoIHNY","Loose"};

  // ----------------------------------------------
  // Make histfiles for lepton optimization studies

  if (toMake == "all" || toMake == "lepOpt"){

    for (int ii = 0; ii < 7; ii++){
      // Run signal
      makeHists("skimTrees/","histfiles_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_semilep","mu",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"fit","nom");
      makeHists("skimTrees/","histfiles_mMu_mEl_"+isoWPs[ii],"PowhegPythia8_semilep","el",false,true,"Medium",isoWPs[ii],false,0.0,false,false,"fit","nom");
      
      //Run QCD
      for (int jj = 13; jj < 19; jj ++){
	makeHists("skimTrees/","histfiles_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"mu",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"fit","nom");
	makeHists("skimTrees/","histfiles_mMu_mEl_"+isoWPs[ii],bkgMCnames[jj],"el",false,false,"Medium",isoWPs[ii],false,0.0,false,false,"fit","nom");    
      }
    }

    // Run signal
    makeHists("skimTrees/","histfiles_tMu_tEl","PowhegPythia8_semilep","mu",false,true,"Tight","MiniIso10",false,0.0,false,false,"fit","nom");
    makeHists("skimTrees/","histfiles_tMu_tEl","PowhegPythia8_semilep","el",false,true,"Tight","MiniIso10",false,0.0,false,false,"fit","nom");
    
    //Run QCD
    for (int jj = 13; jj < 19; jj ++){
      makeHists("skimTrees/","histfiles_tMu_tEl",bkgMCnames[jj],"mu",false,false,"Tight","MiniIso10",false,0.0,false,false,"fit","nom");
      makeHists("skimTrees/","histfiles_tMu_tEl",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",false,0.0,false,false,"fit","nom");    
    }
  }

  // ----------------------------------------------------------------------
  // Make histfiles with final lepton selection, pre-selection-optimization

  if (toMake == "all" || toMake == "selOpt"){
    //run data
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","Data_2015C_mu","mu",true,false,"Medium","MiniIso10",true,0.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","Data_2015C_el","el",true,false,"Tight","MiniIso10",true,0.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","Data_2015D_mu","mu",true,false,"Medium","MiniIso10",true,0.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","Data_2015D_el","el",true,false,"Tight","MiniIso10",true,0.0,true,false,"fit","nom");
    
    // Run signal
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","PowhegPythia8_semilep","mu",false,true,"Medium","MiniIso10",true,0.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10","PowhegPythia8_semilep","el",false,true,"Tight","MiniIso10",true,0.0,true,false,"fit","nom");
    
    //run other MCs
    for (int jj = 0; jj < 19; jj ++){
      makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,0.0,true,false,"fit","nom");
      makeHists("skimTrees/","histfiles_mMu_tEl_MiniIso10",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,0.0,true,false,"fit","nom");
    }
  }
    
  // ----------------------------------------------
  // Make final histfiles

  if (toMake == "all" || toMake == "final"){
    
    //run data
    makeHists("skimTrees/","histfiles","Data_2015C_mu","mu",true,false,"Medium","MiniIso10",true,35.0,true,false,"fit","nom"); //Signal
    makeHists("skimTrees/","histfiles","Data_2015C_el","el",true,false,"Tight","MiniIso10",true,35.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles","Data_2015D_mu","mu",true,false,"Medium","MiniIso10",true,35.0,true,false,"fit","nom");
    makeHists("skimTrees/","histfiles","Data_2015D_el","el",true,false,"Tight","MiniIso10",true,35.0,true,false,"fit","nom");
    
    makeHists("skimTrees/","histfiles","Data_2015C_mu","mu",true,false,"Medium","MiniIso10",true,35.0,true,true,"fit","nom"); //QCD
    makeHists("skimTrees/","histfiles","Data_2015C_el","el",true,false,"Medium","MiniIso10",true,35.0,true,true,"fit","nom");
    makeHists("skimTrees/","histfiles","Data_2015D_mu","mu",true,false,"Medium","MiniIso10",true,35.0,true,true,"fit","nom");
    makeHists("skimTrees/","histfiles","Data_2015D_el","el",true,false,"Medium","MiniIso10",true,35.0,true,true,"fit","nom");
    
    for (int ii = 0; ii < 13; ii++){
      // Run signal
      makeHists("skimTrees/","histfiles","PowhegPythia8_semilep","mu",false,true,"Medium","MiniIso10",true,35.0,true,false,"fit",sysnames[ii]); //Signal
      makeHists("skimTrees/","histfiles","PowhegPythia8_semilep","el",false,true,"Tight","MiniIso10",true,35.0,true,false,"fit",sysnames[ii]);
      makeHists("skimTrees/","histfiles","PowhegPythia8_semilep","mu",false,true,"Medium","MiniIso10",true,35.0,true,true,"fit",sysnames[ii]); //QCD
      makeHists("skimTrees/","histfiles","PowhegPythia8_semilep","el",false,true,"Medium","MiniIso10",true,35.0,true,true,"fit",sysnames[ii]);
      
      //run other MCs
      for (int jj = 0; jj < 19; jj ++){
	makeHists("skimTrees/","histfiles",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,true,false,"fit",sysnames[ii]);
	makeHists("skimTrees/","histfiles",bkgMCnames[jj],"el",false,false,"Tight","MiniIso10",true,35.0,true,false,"fit",sysnames[ii]);    
	makeHists("skimTrees/","histfiles",bkgMCnames[jj],"mu",false,false,"Medium","MiniIso10",true,35.0,true,true,"fit",sysnames[ii]);
	makeHists("skimTrees/","histfiles",bkgMCnames[jj],"el",false,false,"Medium","MiniIso10",true,35.0,true,true,"fit",sysnames[ii]);    
      }
    }
  }
}
