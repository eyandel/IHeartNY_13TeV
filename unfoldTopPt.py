#!/usr/bin/env python

import math

# -------------------------------------------------------------------------------------
# Script for doing RooUnfold on the ttbar differential cross secion
# -------------------------------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()


# -------------------------------------------------------------------------------------
# input options
# -------------------------------------------------------------------------------------

parser.add_option('--closureTest', metavar='F', action='store_true',
                  default=False,
                  dest='closureTest',
                  help='Run closure test')

parser.add_option('--normalize', metavar='F', action='store_true',
                  default=False,
                  dest='normalize',
                  help='Do normalized differential cross section')

parser.add_option('--systVariation', metavar='F', type='string', action='store',
                  default='nom',
                  dest='syst',
                  help='Run nominal or systematic variation?')

parser.add_option('--bkgSyst', metavar='F', type='string', action='store',
                  default="nom",
                  dest='bkgSyst',
                  help='Use nom/up/down bkg normalization?')

parser.add_option('--lepType', metavar='F', type='string', action='store',
                  default='muon',
                  dest='lepType',
                  help='Lepton type (ele or muon)')

parser.add_option('--usePost', metavar='F', action='store_true',
                  default=False,
                  dest='usePost',
                  help='Use posterior top-tag SF and bkg norm')

parser.add_option('--nbr', metavar='F', type='string', action='store',
                  default='5',
                  dest='nbr',
                  help='')

parser.add_option('--toUnfold', metavar='F', type='string', action='store',
                  default='pt',
                  dest='toUnfold',
                  help='Distribution to unfold (pt or y)')


# -------------------------------------------------------------------------------------
# load options & set plot style
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

import sys

from ROOT import gRandom, TH1, TH1D, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TUnfold
from array import array

gROOT.Macro("rootlogon.C")
gROOT.SetBatch(True)

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0)

gStyle.SetTitleFont(43)
gStyle.SetTitleFont(43, "XYZ")
gStyle.SetTitleSize(30, "XYZ")
gStyle.SetLabelFont(43, "XYZ")
gStyle.SetLabelSize(24, "XYZ")

gStyle.SetPadTopMargin(0.07)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.18)
  
gSystem.Load("RooUnfold/libRooUnfold.so")

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold

# -------------------------------------------------------------------------------------
# cross sections, efficiencies, total number of events
# -------------------------------------------------------------------------------------

# luminosity
lum = 12337.98 # was 12358.75 (pb-1), but missing some data
if options.lepType == "ele":
    lum = 12267.67 # was 12295.65 (pb-1), but missing some data

PowhegPythia8_norm         = 831.76         * lum / 182123200.
PowhegPythia8_nonsemi_norm = 831.76         * lum / 182075200. #missing file
SingleTop_t_t_norm         = 136.02 * 0.322 * lum / 3279200.
SingleTop_tbar_t_norm      = 80.95 * 0.322  * lum / 1682400.
SingleTop_t_tW_norm        = 35.9           * lum / 998400.
SingleTop_tbar_tW_norm     = 35.9           * lum / 985000.
WJets_HT100to200_norm      = 1345.0 * 1.21  * lum / 27529599.
WJets_HT200to400_norm      = 359.7 * 1.21   * lum / 4963240.
WJets_HT400to600_norm      = 48.91 * 1.21   * lum / 1963464.
WJets_HT600to800_norm      = 12.05 * 1.21   * lum / 3722395.
WJets_HT800to1200_norm     = 5.501 * 1.21   * lum / 6314257.
WJets_HT1200to2500_norm    = 1.329 * 1.21   * lum / 6768156.
WJets_HT2500toInf_norm     = 0.03216 * 1.21 * lum / 253561.

## hack to account for that when doing closure test (unfold 1/2 sample with other 1/2), the ttbar sample is split in two 
eff_closure = 1.0
if options.closureTest == True :
    eff_closure = 2.0

    
# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

append = ""
if options.usePost:
    append = "_post"


muOrEl = "mu"
if options.lepType=="ele":
    print ""
    print "UNFOLDING FOR ELECTRON CHANNEL !!!" 
    print ""
    muOrEl = "el"
else:
    print ""
    print "UNFOLDING FOR MUON CHANNEL !!!" 
    print ""
    
if options.lepType=="ele" and not options.closureTest:
    f_data = TFile("histfiles_80X/hists_Data_el.root")
    f_QCD  = TFile("histfiles_80X/hists_Data_el_qcd.root")
elif not options.closureTest:
    f_data = TFile("histfiles_80X/hists_Data_mu.root")
    f_QCD  = TFile("histfiles_80X/hists_Data_mu_qcd.root")
    
# In the below, file named f_..._odd will be the one from which response matrix is extracted from (if closureTest == True) 
if options.closureTest == True : 
    f_ttbar     = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+options.syst+"_even"+append+".root")
    f_ttbar_odd = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+options.syst+"_odd"+append+".root")
else :
    f_ttbar     = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+options.syst+append+".root")

if not options.closureTest:
    f_ttbar_nonsemilep = TFile("histfiles_80X/hists_PowhegPythia8_nonsemilep_"+muOrEl+"_"+options.syst+append+".root")

    f_T_t     = TFile("histfiles_80X/hists_SingleTop_t_t_"+muOrEl+"_nom"+append+".root")
    f_Tbar_t  = TFile("histfiles_80X/hists_SingleTop_tbar_t_"+muOrEl+"_nom"+append+".root")
    f_T_tW    = TFile("histfiles_80X/hists_SingleTop_t_tW_"+muOrEl+"_nom"+append+".root")
    f_Tbar_tW = TFile("histfiles_80X/hists_SingleTop_tbar_tW_"+muOrEl+"_nom"+append+".root")
    
    f_WJets_HT100to200   = TFile("histfiles_80X/hists_WJets_HT100to200_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT200to400   = TFile("histfiles_80X/hists_WJets_HT200to400_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT400to600   = TFile("histfiles_80X/hists_WJets_HT400to600_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT600to800   = TFile("histfiles_80X/hists_WJets_HT600to800_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT800to1200  = TFile("histfiles_80X/hists_WJets_HT800to1200_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT1200to2500 = TFile("histfiles_80X/hists_WJets_HT1200to2500_"+muOrEl+"_nom"+append+".root")
    f_WJets_HT2500toInf  = TFile("histfiles_80X/hists_WJets_HT2500toInf_"+muOrEl+"_nom"+append+".root")
    
    f_qcd_ttbar              = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+options.syst+"_qcd"+append+".root")
    f_qcd_ttbar_nonsemilep   = TFile("histfiles_80X/hists_PowhegPythia8_nonsemilep_"+muOrEl+"_"+options.syst+"_qcd"+append+".root")
    f_qcd_T_t                = TFile("histfiles_80X/hists_SingleTop_t_t_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_Tbar_t             = TFile("histfiles_80X/hists_SingleTop_tbar_t_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_T_tW               = TFile("histfiles_80X/hists_SingleTop_t_tW_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_Tbar_tW            = TFile("histfiles_80X/hists_SingleTop_tbar_tW_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT100to200   = TFile("histfiles_80X/hists_WJets_HT100to200_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT200to400   = TFile("histfiles_80X/hists_WJets_HT200to400_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT400to600   = TFile("histfiles_80X/hists_WJets_HT400to600_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT600to800   = TFile("histfiles_80X/hists_WJets_HT600to800_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT800to1200  = TFile("histfiles_80X/hists_WJets_HT800to1200_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT1200to2500 = TFile("histfiles_80X/hists_WJets_HT1200to2500_"+muOrEl+"_nom_qcd"+append+".root")
    f_qcd_WJets_HT2500toInf  = TFile("histfiles_80X/hists_WJets_HT2500toInf_"+muOrEl+"_nom_qcd"+append+".root")

    
# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

syst_flag = options.syst
if options.bkgSyst != "nom" :
    syst_flag = "Bkg"+options.bkgSyst

if options.closureTest == True:
    response = f_ttbar_odd.Get("response_"+options.toUnfold+options.nbr+"fine")
    response.SetName("response_"+options.toUnfold+"_"+syst_flag)
else :
    response = f_ttbar.Get("response_"+options.toUnfold+options.nbr+"fine")
    response.SetName("response_"+options.toUnfold+"_"+syst_flag)

TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# read actual histograms
# -------------------------------------------------------------------------------------

hTrue = f_ttbar.Get(options.toUnfold+"GenTop"+options.nbr)
hTrue.Sumw2()

if not options.closureTest: 
    hMeas = f_data.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
else :
    hMeas = f_ttbar.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")

hMeas.Sumw2()

if not options.closureTest: 
    hMeas_tt_nonsemi         = f_ttbar_nonsemilep.Get(options.toUnfold+"RecoTop"+options.nbr+"fine").Clone()
    hMeas_T_t                = f_T_t.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_Tbar_t             = f_Tbar_t.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_T_tW               = f_T_tW.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_Tbar_tW            = f_Tbar_tW.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT100to200   = f_WJets_HT100to200.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT200to400   = f_WJets_HT200to400.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT400to600   = f_WJets_HT400to600.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT600to800   = f_WJets_HT600to800.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT800to1200  = f_WJets_HT800to1200.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT1200to2500 = f_WJets_HT1200to2500.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_WJets_HT2500toInf  = f_WJets_HT2500toInf.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    
    hMeas_tt_nonsemi.Sumw2()
    hMeas_T_t.Sumw2()
    hMeas_Tbar_t.Sumw2()
    hMeas_T_tW.Sumw2()
    hMeas_Tbar_tW.Sumw2()
    hMeas_WJets_HT100to200.Sumw2()
    hMeas_WJets_HT200to400.Sumw2()
    hMeas_WJets_HT400to600.Sumw2()
    hMeas_WJets_HT600to800.Sumw2()
    hMeas_WJets_HT800to1200.Sumw2()
    hMeas_WJets_HT1200to2500.Sumw2()
    hMeas_WJets_HT2500toInf.Sumw2()
    
    hMeas_qcd                    = f_QCD.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_tt_nonsemi         = f_qcd_ttbar_nonsemilep.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_T_t                = f_qcd_T_t.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_Tbar_t             = f_qcd_Tbar_t.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_T_tW               = f_qcd_T_tW.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_Tbar_tW            = f_qcd_Tbar_tW.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT100to200   = f_qcd_WJets_HT100to200.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT200to400   = f_qcd_WJets_HT200to400.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT400to600   = f_qcd_WJets_HT400to600.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT600to800   = f_qcd_WJets_HT600to800.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT800to1200  = f_qcd_WJets_HT800to1200.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT1200to2500 = f_qcd_WJets_HT1200to2500.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hMeas_qcd_WJets_HT2500toInf  = f_qcd_WJets_HT2500toInf.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    
    hMeas_qcd.Sumw2()
    hMeas_qcd_tt_nonsemi.Sumw2()
    hMeas_qcd_T_t.Sumw2()
    hMeas_qcd_Tbar_t.Sumw2()
    hMeas_qcd_T_tW.Sumw2()
    hMeas_qcd_Tbar_tW.Sumw2()
    hMeas_qcd_WJets_HT100to200.Sumw2()
    hMeas_qcd_WJets_HT200to400.Sumw2()
    hMeas_qcd_WJets_HT400to600.Sumw2()
    hMeas_qcd_WJets_HT600to800.Sumw2()
    hMeas_qcd_WJets_HT800to1200.Sumw2()
    hMeas_qcd_WJets_HT1200to2500.Sumw2()
    hMeas_qcd_WJets_HT2500toInf.Sumw2()

    
# -------------------------------------------------------------------------------------
# Normalize histograms
# -------------------------------------------------------------------------------------

### normalize measured "data" for closure tests
if options.closureTest: 
    hMeas.Scale(PowhegPythia8_norm * eff_closure)
    hMeas.SetName(options.toUnfold+"RecoTop_measured")
    
hTrue.Scale(PowhegPythia8_norm * eff_closure)
hTrue.SetName(options.toUnfold+"GenTop_true")


if not options.closureTest: 
    hMeas_tt_nonsemi.Scale(PowhegPythia8_nonsemi_norm)
    
    hMeas_T_t.Scale( SingleTop_t_t_norm)
    hMeas_Tbar_t.Scale( SingleTop_tbar_t_norm)
    hMeas_T_tW.Scale(SingleTop_t_tW_norm)
    hMeas_Tbar_tW.Scale(SingleTop_tbar_tW_norm)
    
    hMeas_WJets_HT100to200.Scale(WJets_HT100to200_norm)
    hMeas_WJets_HT200to400.Scale(WJets_HT200to400_norm)
    hMeas_WJets_HT400to600.Scale(WJets_HT400to600_norm)
    hMeas_WJets_HT600to800.Scale(WJets_HT600to800_norm)
    hMeas_WJets_HT800to1200.Scale(WJets_HT800to1200_norm)
    hMeas_WJets_HT1200to2500.Scale(WJets_HT1200to2500_norm)
    hMeas_WJets_HT2500toInf.Scale(WJets_HT2500toInf_norm)
    
    hMeas_qcd_tt_nonsemi.Scale(PowhegPythia8_nonsemi_norm)
    hMeas_qcd_T_t.Scale( SingleTop_t_t_norm)
    hMeas_qcd_Tbar_t.Scale( SingleTop_tbar_t_norm)
    hMeas_qcd_T_tW.Scale(SingleTop_t_tW_norm)
    hMeas_qcd_Tbar_tW.Scale(SingleTop_tbar_tW_norm)
    hMeas_qcd_WJets_HT100to200.Scale(WJets_HT100to200_norm)
    hMeas_qcd_WJets_HT200to400.Scale(WJets_HT200to400_norm)
    hMeas_qcd_WJets_HT400to600.Scale(WJets_HT400to600_norm)
    hMeas_qcd_WJets_HT600to800.Scale(WJets_HT600to800_norm)
    hMeas_qcd_WJets_HT800to1200.Scale(WJets_HT800to1200_norm)
    hMeas_qcd_WJets_HT1200to2500.Scale(WJets_HT1200to2500_norm)
    hMeas_qcd_WJets_HT2500toInf.Scale(WJets_HT2500toInf_norm)
    
    hMeas_tt_nonsemi.SetName(options.toUnfold+"RecoTop_TTnonSemi")
    
    hMeas_SingleTop = hMeas_T_t.Clone()
    hMeas_SingleTop.SetName(options.toUnfold+"RecoTop_SingleTop")
    
    hMeas_WJets = hMeas_WJets_HT100to200.Clone()
    hMeas_WJets.SetName(options.toUnfold+"RecoTop_WJets")
    
    hMeas_QCD = hMeas_qcd.Clone()
    hMeas_QCD.SetName(options.toUnfold+"RecoTop_QCD")
    
    for hist in [hMeas_Tbar_t, hMeas_T_tW, hMeas_Tbar_tW] :
        hMeas_SingleTop.Add( hist )
        
    for hist in [hMeas_WJets_HT200to400,hMeas_WJets_HT400to600,hMeas_WJets_HT600to800,hMeas_WJets_HT800to1200,hMeas_WJets_HT1200to2500,hMeas_WJets_HT2500toInf] :
        hMeas_WJets.Add( hist )

    for hist in [hMeas_qcd_tt_nonsemi,
                 hMeas_qcd_T_t,hMeas_qcd_Tbar_t,hMeas_qcd_T_tW,hMeas_qcd_Tbar_tW
                 ,hMeas_qcd_WJets_HT100to200,hMeas_qcd_WJets_HT200to400,hMeas_qcd_WJets_HT400to600,
                 hMeas_qcd_WJets_HT600to800,hMeas_qcd_WJets_HT800to1200,hMeas_qcd_WJets_HT1200to2500,hMeas_qcd_WJets_HT2500toInf] :
        hMeas_QCD.Add(hist,-1.0)

    # -------------------------------
    # Normalize QCD to MC prediction
    # -------------------------------

    f_QCD_HT300to500   = TFile("histfiles_80X/hists_QCD_HT300to500_"+muOrEl+"_nom"+append+".root")
    f_QCD_HT500to700   = TFile("histfiles_80X/hists_QCD_HT500to700_"+muOrEl+"_nom"+append+".root")
    f_QCD_HT700to1000  = TFile("histfiles_80X/hists_QCD_HT700to1000_"+muOrEl+"_nom"+append+".root")
    f_QCD_HT1000to1500 = TFile("histfiles_80X/hists_QCD_HT1000to1500_"+muOrEl+"_nom"+append+".root")
    f_QCD_HT1500to2000 = TFile("histfiles_80X/hists_QCD_HT1500to2000_"+muOrEl+"_nom"+append+".root")
    f_QCD_HT2000toInf  = TFile("histfiles_80X/hists_QCD_HT2000toInf_"+muOrEl+"_nom"+append+".root")

    hNorm_QCD_HT300to500   = f_QCD_HT300to500.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hNorm_QCD_HT500to700   = f_QCD_HT500to700.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hNorm_QCD_HT700to1000  = f_QCD_HT700to1000.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hNorm_QCD_HT1000to1500 = f_QCD_HT1000to1500.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hNorm_QCD_HT1500to2000 = f_QCD_HT1500to2000.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")
    hNorm_QCD_HT2000toInf  = f_QCD_HT2000toInf.Get(options.toUnfold+"RecoTop"+options.nbr+"fine")

    hNorm_QCD_HT300to500.Scale(347700. * lum / 37828442.)
    hNorm_QCD_HT500to700.Scale(32100. * lum / 44058594.)
    hNorm_QCD_HT700to1000.Scale(6831. * lum / 29832311.)
    hNorm_QCD_HT1000to1500.Scale(1207. * lum / 4963881.) # was 4980387 -- missing file
    hNorm_QCD_HT1500to2000.Scale(119.9 * lum / 7803965.)
    hNorm_QCD_HT2000toInf.Scale(25.24 * lum / 4047532.)

    QCD_norm = 0.0

    for hist in [hNorm_QCD_HT300to500,hNorm_QCD_HT500to700,hNorm_QCD_HT700to1000,hNorm_QCD_HT1000to1500,hNorm_QCD_HT1500to2000,hNorm_QCD_HT2000toInf]:
        QCD_norm += hist.Integral()
    
    hMeas_QCD.Scale(QCD_norm / hMeas_QCD.Integral())

    # -------------------------------------------------------------------------------------
    # Close histogram files
    # -------------------------------------------------------------------------------------

    f_data.Close()
    f_QCD.Close()
    f_ttbar.Close()
    
    if options.closureTest : 
        f_ttbar_odd.Close()

    else :
        for file in [f_ttbar_nonsemilep, f_T_t, f_Tbar_t, f_T_tW, f_Tbar_tW, f_WJets_HT100to200, f_WJets_HT200to400, f_WJets_HT400to600,
                     f_WJets_HT600to800, f_WJets_HT800to1200, f_WJets_HT1200to2500, f_qcd_ttbar, f_qcd_ttbar_nonsemilep, f_qcd_T_t,
                     f_qcd_Tbar_t, f_qcd_T_tW, f_qcd_Tbar_tW, f_qcd_WJets_HT100to200, f_qcd_WJets_HT200to400, f_qcd_WJets_HT400to600,
                     f_qcd_WJets_HT600to800, f_qcd_WJets_HT800to1200, f_qcd_WJets_HT1200to2500, f_qcd_WJets_HT2500toInf, f_QCD_HT300to500,
                     f_QCD_HT500to700, f_QCD_HT700to1000, f_QCD_HT1000to1500, f_QCD_HT1500to2000, f_QCD_HT2000toInf] :
            file.Close()
        
    # -------------------------------------------------------------------------------------
    # Scale backgrounds if using posterior normalization
    # -------------------------------------------------------------------------------------

    if options.usePost:
        if options.bkgSyst == "Up":
            if options.lepType == "ele" :
                hMeas_tt_nonsemi.Scale(0.83*1.79) 
                hMeas_SingleTop.Scale(1.05*1.79) 
                hMeas_WJets.Scale(1.08*1.79)
                hMeas_QCD.Scale(1.02*1.79)
            else :
                hMeas_tt_nonsemi.Scale(0.83) 
                hMeas_SingleTop.Scale(1.05) 
                hMeas_WJets.Scale(1.08)
                hMeas_QCD.Scale(0.79)
                
        elif options.bkgSyst == "Down":
            if options.lepType == "ele":
                hMeas_tt_nonsemi.Scale(0.71*1.79)
                hMeas_SingleTop.Scale(0.95*1.79) 
                hMeas_WJets.Scale(0.96*1.79)
                hMeas_QCD.Scale(0.32*1.79)
            else :
                hMeas_tt_nonsemi.Scale(0.71)
                hMeas_SingleTop.Scale(0.95) 
                hMeas_WJets.Scale(0.96)
                hMeas_QCD.Scale(0.45)

        else :
            if options.lepType == "ele":
                hMeas_tt_nonsemi.Scale(0.77*1.79)
                hMeas_SingleTop.Scale(1.00*1.79) 
                hMeas_WJets.Scale(1.02*1.79)
                hMeas_QCD.Scale(0.67*1.79)
            else :
                hMeas_tt_nonsemi.Scale(0.77)
                hMeas_SingleTop.Scale(1.00) 
                hMeas_WJets.Scale(1.02)
                hMeas_QCD.Scale(0.62)

    # -------------------------------------------------------------------------------------
    # subtract backgrounds from the data distribution, but not for closure test!!! 
    # -------------------------------------------------------------------------------------

    for hist in [hMeas_SingleTop, hMeas_WJets, hMeas_QCD, hMeas_tt_nonsemi] :
        hMeas.Add(hist, -1.)

    for ibin in xrange( hMeas.GetNbinsX() ) :
        if ( hMeas.GetBinContent( ibin ) < 0.0 ) :
            hMeas.SetBinContent( ibin, 0.0 )

# -------------------------------------------------------------------------------------
# output file for storing histograms to  
# -------------------------------------------------------------------------------------

closureout = ""
if options.closureTest == True: 
    closureout = "_closure"
norm_flag = ""
if options.normalize:
    norm_flag = "_norm"
        
fout = TFile("UnfoldingPlots/unfold_"+options.toUnfold+"_PowhegPythia8_"+options.lepType+"_"+syst_flag+closureout+norm_flag+".root","recreate")
            
# -------------------------------------------------------------------------------------
# do the actual unfolding
# -------------------------------------------------------------------------------------

print "------------ UNFOLDING (syst: " + syst_flag + ") ------------"

n_iter = 3
print " using " + str(n_iter) + " iterations"

unfold = RooUnfoldBayes(response, hMeas, n_iter)
#unfold = RooUnfoldSvd(response, hMeas, 2);

# get the unfolded distribution
hReco = unfold.Hreco()
hReco.Sumw2()

# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------
# TODO: should fix BR

print "WARNING: treatment of branching ratio is wrong!"

hTrue.Scale(1.0/(lum*0.438/3.)) # true @ parton level
hMeas.Scale(1.0/(lum*0.438/3.)) # measured @ reco level
hReco.Scale(1.0/(lum*0.438/3.)) # unfolded to parton level


# -------------------------------------------------------------------------------------
# Adjust for bin width
# -------------------------------------------------------------------------------------

sumReco = 0
sumTrue = 0
sumMeas = 0
sumReco_odd = 0
sumTrue_odd = 0
sumMeas_odd = 0

lowedge = 399.
highedge = 1199.

for ibin in range(1, hTrue.GetXaxis().GetNbins()+1 ) :
        
    # total cross section for pt > 400 GeV
    if hTrue.GetBinLowEdge(ibin) > lowedge :
        sumTrue += hTrue.GetBinContent(ibin)
    if hReco.GetBinLowEdge(ibin) > lowedge :
        sumReco += hReco.GetBinContent(ibin)
    if hMeas.GetBinLowEdge(ibin) > lowedge :
        sumMeas += hMeas.GetBinContent(ibin)

    width = hTrue.GetBinWidth(ibin)

    hTrue.SetBinContent(ibin, hTrue.GetBinContent(ibin) / width )
    hTrue.SetBinError(ibin, hTrue.GetBinError(ibin) / width )

    hMeas.SetBinContent(ibin,  hMeas.GetBinContent(ibin) / width )
    hMeas.SetBinError(ibin,  hMeas.GetBinError(ibin) / width )
        
    hReco.SetBinContent(ibin, hReco.GetBinContent(ibin) / width )
    hReco.SetBinError(ibin, hReco.GetBinError(ibin) / width )

# -------------------------------------------------------------------------------------
# do NORMALIZED differential cross section instead??
# -------------------------------------------------------------------------------------

if options.normalize: 
    hTrue.Scale(1.0/sumTrue)
    hReco.Scale(1.0/sumReco)
    hMeas.Scale(1.0/sumMeas)

# -------------------------------------------------------------------------------------
# print & draw
# -------------------------------------------------------------------------------------

## ratio of unfolded data to generator-level 
hFrac = hReco.Clone()
hFrac.SetName("hFrac")
hFrac.SetTitle(";Top quark p_{T} (GeV);Data/MC")
hFrac.Divide(hTrue)

if options.closureTest:
    print ''
    print '-------------------------------------------------------------------------------------'
    print 'uncertainty from closure test (' + options.lepType + ')'
    print '-------------------------------------------------------------------------------------'
    print 'parton-level'
    for ibin in range(1, hFrac.GetXaxis().GetNbins()+1 ) :
        if hFrac.GetBinLowEdge(ibin) > lowedge and hFrac.GetBinLowEdge(ibin) < highedge:
            print '[' + str(hFrac.GetBinLowEdge(ibin)) + ',' + str(hFrac.GetBinLowEdge(ibin+1)) + '] = ' + str((hFrac.GetBinContent(ibin)-1.0)*100.0) + ' %'            

print ''
print '-------------------------------------------------------------------------------------'
print 'sigma (raw data) = ' + str(int(sumMeas)) + ' pb'
print 'true sigma @ parton-level (pt > 400 GeV)      = ' + str(int(sumTrue)) + ' pb'
print 'measured sigma (unfolded data) @ parton-level = ' + str(int(sumReco)) + ' pb'
print '-------------------------------------------------------------------------------------'
print ''

# -------------------------------------------------------------------------------------
# draw parton-level unfolding
# -------------------------------------------------------------------------------------

c1 = TCanvas("c", "c", 700, 700)
pad1 =  TPad("pad1","pad1",0,0.3,1,1)
pad1.SetBottomMargin(0.05)
pad1.Draw()
pad1.cd()

hReco.SetMarkerStyle(21)
hMeas.SetMarkerStyle(25)

hReco.GetXaxis().SetRangeUser(400.,1200.)
hTrue.GetXaxis().SetRangeUser(400.,1200.)
hMeas.GetXaxis().SetRangeUser(400.,1200.)

xsec_title = ";;d#sigma/dp_{T} [fb/GeV]"
if options.normalize:    
    xsec_title = ";;1/#sigma d#sigma/dp_{T} [1/GeV]"

hReco.SetTitle(xsec_title)
hReco.GetYaxis().SetTitleOffset(1.2)
hReco.SetMinimum(0.0)
max = hTrue.GetMaximum()
max2 = hReco.GetMaximum()
if max2 > max:
	max = max2
hReco.SetAxisRange(0,max*1.15,"Y")
hReco.Draw()
hTrue.Draw('hist same')
hMeas.Draw('same')
hTrue.UseCurrentStyle()
hTrue.SetLineColor(4)
hTrue.GetYaxis().SetTitleSize(25)
hTrue.GetXaxis().SetLabelSize(0)

leg = TLegend(0.5, 0.55, 0.9, 0.75)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)

tt = TLatex()
tt.SetNDC()
tt.SetTextFont(42)

if options.closureTest == False : 
    leg.AddEntry( hReco, 'Unfolded data', 'p')
    leg.AddEntry( hTrue, 'Generated (Powheg)', 'l')
    leg.AddEntry( hMeas, 'Raw data', 'p')
else : 
    leg.AddEntry( hReco, 'Unfolded MC (Powheg)', 'p')
    leg.AddEntry( hTrue, 'Generated (Powheg)', 'l')
    leg.AddEntry( hMeas, 'Reco-level (Powheg)', 'p')
    tt.DrawLatex(0.55,0.45, "MC closure test")
    
leg.Draw()

# write histograms to file
if options.closureTest:
    hReco.SetName("UnfoldedMC")
else:
    hReco.SetName("UnfoldedData")

hReco.Write()
hTrue.Write()
hMeas.Write()

text1 = TLatex()
text1.SetNDC()
text1.SetTextFont(42)
text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 12.4 fb^{-1}, #sqrt{s} = 13 TeV}")

c1.cd()
pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
pad2.SetTopMargin(0.05)
pad2.SetBottomMargin(0.4)
pad2.Draw()
pad2.cd()
pad2.SetGridy()
hFrac.SetMaximum(1.4)
hFrac.SetMinimum(0.6)
if options.normalize == False and options.closureTest == False: 
    hFrac.SetMaximum(1.2)
    hFrac.SetMinimum(0.4)    
hFrac.UseCurrentStyle()
hFrac.GetYaxis().SetTitleSize(25)
hFrac.GetYaxis().SetTitleOffset(2.0)
hFrac.GetXaxis().SetTitleOffset(4.0)
hFrac.GetXaxis().SetLabelSize(25)
hFrac.GetYaxis().SetNdivisions(4,4,0,False)

hFrac.Draw("e")
hFrac.GetXaxis().SetRangeUser(400., 1200.)

c1.Update()

if options.syst=="nom" and options.bkgSyst == "nom":
    c1.Print("UnfoldingPlots/unfolded_ttbar_xs"+closureout+"_"+options.lepType+"_"+options.syst+norm_flag+".pdf", "pdf")

# -------------------------------------------------------------------------------------
# plot response matrices 
# (do this in the end as the normalization otherwise will mess up the unfolding result!)
# -------------------------------------------------------------------------------------

ncontours = 256
stops = [0.00, 1.00]
red   = [0.99, 0.32]
green = [0.99, 0.42]
blue  = [0.99, 0.9]
s = array('d', stops)
r = array('d', red)
g = array('d', green)
b = array('d', blue)
npoints = len(s)
TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
gStyle.SetNumberContours(ncontours)

# -------------------------------------------------------------------------------------
# one-step unfolding
# -------------------------------------------------------------------------------------

gStyle.SetPadRightMargin(0.12)
cr = TCanvas("c_response", "", 800, 600)

hEmpty2D = response.Hresponse().Clone()
hEmpty2D.SetName("empty2D")
hEmpty2D.Reset()
hEmpty2D.GetXaxis().SetTitle("Reconstructed top jet p_{T} (GeV)")
hEmpty2D.GetYaxis().SetTitle("Top quark p_{T} (GeV)")
hEmpty2D.GetXaxis().SetLabelSize(0.045)
hEmpty2D.GetYaxis().SetLabelSize(0.045)
hEmpty2D.Draw()

hResponse2D = response.Hresponse().Clone()
hResponse2D.SetName("plottedResponse")

gStyle.SetPaintTextFormat(".1f")
hResponse2D.Draw("colz,same,text")
hEmpty2D.Draw("axis,same")

#if options.syst=="nom" and options.closureTest==False:
cr.SaveAs("UnfoldingPlots/unfold"+closureout+"_responseMatrix_full_"+options.lepType+"_"+syst_flag+".pdf")

# normalize so that for each bin of true top quark pt(eta), the bins in measured top pt(eta) add up to 100%
nbinsX = hResponse2D.GetNbinsX()
nbinsY = hResponse2D.GetNbinsY()
for iby in range(1,nbinsY+1) :
    rowIntegral = hResponse2D.Integral(1,nbinsX,iby,iby)
    for ibx in range(1,nbinsX+1) :
        binContent = hResponse2D.GetBinContent(ibx,iby)
        newContent = 0
        if rowIntegral > 0:
            newContent = binContent/rowIntegral*100.0
        hResponse2D.SetBinContent(ibx,iby,newContent)

hEmpty2D.Draw()
hResponse2D.Draw("colz,same,text")
hEmpty2D.Draw("axis,same")
#if options.syst=="nom" and options.closureTest==False:
cr.SaveAs("UnfoldingPlots/unfold"+closureout+"_responseMatrix_"+options.lepType+"_"+syst_flag+".pdf")

hEmpty2D.SetAxisRange(410,1150,"X")
hEmpty2D.SetAxisRange(410,1150,"Y")
hResponse2D.SetAxisRange(410,1150,"X")
hResponse2D.SetAxisRange(410,1150,"Y")
hEmpty2D.Draw()
hResponse2D.Draw("colz,same,text")
hEmpty2D.Draw("axis,same")

#if options.syst=="nom" and options.closureTest==False:
cr.SaveAs("UnfoldingPlots/unfold"+closureout+"_responseMatrix_zoom_"+options.lepType+"_"+syst_flag+".pdf")
        
response.Hresponse().SetName("responseMatrix_"+syst_flag)
response.Hresponse().Write()
    

fout.Close()
