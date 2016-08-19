// Script to run makePlots.cc
// Syntax is makePlots(TString channel, TString var, TString region)

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void runMakePlots(TString toPlot = "final"){

  gSystem->CompileMacro("makePlots.cc");
  cout << "Compilation successful!" << endl;

  TString channels[2] = {"mu","el"};
  const int nregion = 4;
  TString regions[nregion] = {"Pre","0t","1t0b","1t1b"};
  const int nhist = 28;
  TString hists[nhist] = {"metPt",
			  "ht",
			  "htLep",
			  "ak4jetPt",
			  "ak4jetEta",
			  "ak4jetPhi",
			  "ak4jetMass",
			  "ak4jetCSV",
			  "ak4jetVtxMass",
			  "ak8jetPt",
			  "ak8jetEta",
			  "ak8jetPhi",
			  "ak8jetY",
			  "ak8jetMass",
			  "ak8jetTau1",
			  "ak8jetTau2",
			  "ak8jetTau3",
			  "ak8jetTau32",
			  "ak8jetTau21",
			  "ak8jetCSV",
			  "ak8jetSDmass",
			  "ak8jetSDsubjetMaxCSV",
			  "ak8jetSDm01",
			  "lepPt",
			  "lepEta",
			  "lepAbsEta",
			  "lepPhi",
			  //"lepBJetdR",
			  //"lepTJetdR",
			  //"lepBJetPtRel",
			  "lepMiniIso"
  };

  if (toPlot == "all" || toPlot == "lepSelOpt"){
    makeEffPlots("mu","ak8jetPt");
    makeEffPlots("el","ak8jetPt");
    makeEffPlots("mu","leadLepPt");
    makeEffPlots("el","leadLepPt");
    drawROCCurve("mu");
    drawROCCurve("el");
    makeEffPlotsFinal();
    
    cout << endl << "Finished with lepton selection optimization plots!" << endl << endl;
  }

  if (toPlot == "all" || toPlot == "selOpt"){

    make1DScans("histfiles_80X_mMu_tEl_MiniIso10/","mu");
    make1DScans("histfiles_80X_mMu_tEl_MiniIso10/","el");
    
    TString opthists[5] = {"metPt","ht","htLep","lepBJetdR","lepTJetdR"};
    
    for (int ii = 0; ii < 5; ii++){
      makePlots("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","mu",opthists[ii],"Pre",true,false,false);
      makePlots("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","el",opthists[ii],"Pre",true,false,false);
    }

    makeTable("histfiles_80X_mMu_tEl_MiniIso10/","histfiles_80X_mMu_tEl_MiniIso10/","el",false,true,false);
    
    cout << endl << "Finished with selection optimization plots!" << endl << endl;
  }

  if (toPlot == "all" || toPlot == "final"){
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  makePlots("histfiles_80X/","histfiles_80X/",channels[ii],hists[kk],regions[jj],false,false,false);
	}
      }
    }
    
    makePlots("histfiles_80X/","histfiles_80X/","mu","nAK4jet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","el","nAK4jet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nBjet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","el","nBjet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nAK8jet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","el","nAK8jet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nTjet","Pre",false,false,false);
    makePlots("histfiles_80X/","histfiles_80X/","el","nTjet","Pre",false,false,false);
    
    cout << endl << "Finished with regular plots!" << endl << endl;

    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nhist; jj++){
	makeQCDComp("histfiles_80X/","histfiles_80X/",channels[ii],hists[jj]);
      }
      makeQCDComp("histfiles_80X/","histfiles_80X/",channels[ii],"nAK4jet");
      makeQCDComp("histfiles_80X/","histfiles_80X/",channels[ii],"nBjet");
      makeQCDComp("histfiles_80X/","histfiles_80X/",channels[ii],"nAK8jet");
      makeQCDComp("histfiles_80X/","histfiles_80X/",channels[ii],"nTjet");
    }
    
    cout << endl << "Finished with QCD comparison plots!" << endl << endl;
    
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  compareShapes("histfiles_80X/","histfiles_80X/",channels[ii],hists[kk],regions[jj]);
	}
      }
    }

    cout << endl << "Finished making shape comparisons!" << endl << endl;
    
    //makeCombineInputs("histfiles_80X/","histfiles_80X/");
    
    //cout << endl << "Finished making combine inputs!" << endl << endl;

    makeTable("histfiles_80X/","histfiles_80X/","mu",true,true,false);
    cout << endl;
    makeTable("histfiles_80X/","histfiles_80X/","el",true,true,false);
    cout << endl;
    makeTable("histfiles_80X/","histfiles_80X/","mu",false,false,false);
    cout << endl;
    makeTable("histfiles_80X/","histfiles_80X/","el",false,false,false);
  }

  if (toPlot == "post"){
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  makePlots("histfiles_80X/","histfiles_80X/",channels[ii],hists[kk],regions[jj],false,true,true);
	}
      }
    }
    
    makePlots("histfiles_80X/","histfiles_80X/","mu","nAK4jet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","el","nAK4jet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nBjet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","el","nBjet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nAK8jet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","el","nAK8jet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","mu","nTjet","Pre",false,true,true);
    makePlots("histfiles_80X/","histfiles_80X/","el","nTjet","Pre",false,true,true);
    
    cout << endl << "Finished with postfit plots!" << endl << endl;
  }

  if (toPlot == "combine"){
    combineResults("mu","muA","pre");
    combineResults("mu","muA","post");
    combineResults("mu","muB","pre");
    combineResults("mu","muB","post");
    combineResults("mu","muC","pre");
    combineResults("mu","muC","post");
    combineResults("el","elA","pre");
    combineResults("el","elA","post");
    combineResults("el","elB","pre");
    combineResults("el","elB","post");
    combineResults("el","elC","pre");
    combineResults("el","elC","post");
    combineResults("mu","combA","pre");
    combineResults("mu","combA","post");
    combineResults("mu","combB","pre");
    combineResults("mu","combB","post");
    combineResults("mu","combC","pre");
    combineResults("mu","combC","post");
    combineResults("el","combA","pre");
    combineResults("el","combA","post");
    combineResults("el","combB","pre");
    combineResults("el","combB","post");
    combineResults("el","combC","pre");
    combineResults("el","combC","post");
  }  
}
