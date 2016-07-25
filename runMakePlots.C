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
  const int nhist = 31;
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
			  "lepBJetdR",
			  "lepTJetdR",
			  "lepBJetPtRel",
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

    make1DScans("histfiles_mMu_tEl_MiniIso10/","mu");
    make1DScans("histfiles_mMu_tEl_MiniIso10/","el");
    
    TString opthists[5] = {"metPt","ht","htLep","lepBJetdR","lepTJetdR"};
    
    for (int ii = 0; ii < 5; ii++){
      makePlots("histfiles_mMu_tEl_MiniIso10/","histfiles_mMu_tEl_MiniIso10/","mu",opthists[ii],"Pre",true,false);
      makePlots("histfiles_mMu_tEl_MiniIso10/","histfiles_mMu_tEl_MiniIso10/","el",opthists[ii],"Pre",true,false);
    }

    makeTable("histfiles_mMu_tEl_MiniIso10/","histfiles_mMu_tEl_MiniIso10/","el",false,true);
    
    cout << endl << "Finished with selection optimization plots!" << endl << endl;
  }

  if (toPlot == "all" || toPlot == "final"){
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  makePlots("histfiles/","histfiles/",channels[ii],hists[kk],regions[jj],false,true);
	}
      }
    }
    
    makePlots("histfiles/","histfiles/","mu","nAK4jet","Pre",false,true);
    makePlots("histfiles/","histfiles/","el","nAK4jet","Pre",false,true);
    makePlots("histfiles/","histfiles/","mu","nBjet","Pre",false,true);
    makePlots("histfiles/","histfiles/","el","nBjet","Pre",false,true);
    makePlots("histfiles/","histfiles/","mu","nAK8jet","Pre",false,true);
    makePlots("histfiles/","histfiles/","el","nAK8jet","Pre",false,true);
    makePlots("histfiles/","histfiles/","mu","nTjet","Pre",false,true);
    makePlots("histfiles/","histfiles/","el","nTjet","Pre",false,true);
    
    cout << endl << "Finished with regular plots!" << endl << endl;

    /*
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < 30; jj++){
	makeQCDComp("histfiles/","histfiles/",channels[ii],hists[jj]);
      }
      makeQCDComp("histfiles/","histfiles/",channels[ii],"nAK4jet");
      makeQCDComp("histfiles/","histfiles/",channels[ii],"nBjet");
      makeQCDComp("histfiles/","histfiles/",channels[ii],"nAK8jet");
      makeQCDComp("histfiles/","histfiles/",channels[ii],"nTjet");
    }
    
    cout << endl << "Finished with QCD comparison plots!" << endl << endl;
    */
    
    for (int ii = 0; ii < 2; ii++){
      for (int jj = 0; jj < nregion; jj++){
	for (int kk = 0; kk < nhist; kk++){
	  compareShapes("histfiles/","histfiles/",channels[ii],hists[kk],regions[jj]);
	}
      }
    }

    cout << endl << "Finished maing shape comparisons!" << endl << endl;
    
    makeCombineInputs("histfiles/","histfiles/","mu");
    //makeCombineInputs("histfiles/","histfiles/","el");
    
    cout << endl << "Finished making combine inputs!" << endl << endl;

    makeTable("histfiles/","histfiles/","mu",true,true);
    cout << endl;
    makeTable("histfiles/","histfiles/","el",true,true);
    cout << endl;
    makeTable("histfiles/","histfiles/","mu",false,false);
    cout << endl;
    makeTable("histfiles/","histfiles/","el",false,false);
  }

  if (toPlot == "combine"){
    combineResults("mu","mu","pre");
    combineResults("mu","mu","post");
  }
  
}
