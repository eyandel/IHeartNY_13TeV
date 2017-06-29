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

parser.add_option('--lepType', metavar='F', type='string', action='store',
                  default='muon',
                  dest='lepType',
                  help='Lepton type (ele or muon)')

parser.add_option('--type', metavar='F', type='string', action='store',
                  default='full',
                  dest='type',
                  help='')

parser.add_option('--toy', metavar='F', type='string', action='store',
                  default='',
                  dest='toy',
                  help='')

parser.add_option('--regType', metavar='F', type='string', action='store',
                  default='scanTau',
                  dest='regType',
                  help='scanTau or LCurve')

# -------------------------------------------------------------------------------------
# load options & set plot style
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

import sys

from ROOT import gRandom, TH1, TH1D, TH1F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD, TGraph, TUnfold, Double, TSpline, TSpline3, TUnfoldDensity, TAttLine
from array import array

gROOT.Macro("rootlogon.C")
gROOT.SetBatch(True)

gStyle.SetOptStat(000000)
gStyle.SetOptTitle(0);

gStyle.SetTitleFont(43)
gStyle.SetTitleFont(43, "XYZ")
gStyle.SetTitleSize(30, "XYZ")
gStyle.SetLabelFont(43, "XYZ")
gStyle.SetLabelSize(24, "XYZ")

gStyle.SetPadTopMargin(0.07);
gStyle.SetPadRightMargin(0.05);
gStyle.SetPadBottomMargin(0.16);
gStyle.SetPadLeftMargin(0.18);
  
gSystem.Load("RooUnfold/libRooUnfold.so")

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold

print "TUnfold version is " + str(TUnfold.GetTUnfoldVersion())

# -------------------------------------------------------------------------------------
# Define helper functions
# -------------------------------------------------------------------------------------

def makeResponse( response ) :
    Hres = response.HresponseNoOverflow()
    vtru = response.Vtruth()
    etru = response.Etruth()
        
    for j in xrange(1,response.GetNbinsTruth()+1) :
        ntru = 0.0
        for i in xrange(1,response.GetNbinsMeasured()+1) :
            ntru += Hres.GetBinContent(i,j)
        Hres.SetBinContent(0, j, vtru[j-1]-ntru)
        Hres.SetBinError(0, j, etru[j-1])
        
    return Hres

def removeFakes( hMeas, response ):
    if response.FakeEntries() :
        fakes = response.Vfakes()
        fac = response.Vmeasured().Sum() # Measured, from response matrix
        if fac != 0.0 :
            measVec = RooUnfoldResponse.H2V(hMeas,response.GetNbinsMeasured(),0) #Actual measured input
            fac = measVec.Sum() / fac 
            for i in xrange(1,response.GetNbinsMeasured()+1) :
                hMeas.SetBinContent(i,hMeas.GetBinContent(i) - fac * fakes[i-1])

lum = 35867.0

# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"

f_ttbar        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_post.root")
f_ttbar_even   = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_even_post.root")
f_ttbar_odd    = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_odd_post.root")
f_ttbar_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_post.root")
f_ttbar_even_2 = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_even_post.root")
f_ttbar_odd_2  = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_nom_odd_post.root")

# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

if options.type == "full":
    response = f_ttbar.Get("response2_pt_split")
    response2 = f_ttbar_2.Get("response2_pt_split")
elif options.type == "half" :
    response = f_ttbar_odd.Get("response2_pt_split")
    response2 = f_ttbar_odd_2.Get("response2_pt_split")
elif options.type == "each" : #unfold p1 with p2
    response = f_ttbar_2.Get("response2_pt_split")
else :
    response = f_ttbar_odd.Get("response2_pt_split")
    
TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# Get systematic variations
# -------------------------------------------------------------------------------------
Hres_sys = {}
sysnames = ['JECUp','JECDown','JERUp','JERDown','BTagUp','BTagDown','TopTagUp','TopTagDown','lepUp','lepDown','PDFUp','PDFDown','Q2Up','Q2Down','ASUp','ASDown']

for sysname in sysnames:
    f_ttbar_sys        = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+sysname+"_post.root")
    f_ttbar_sys_odd    = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_"+sysname+"_odd_post.root")
    f_ttbar_sys_2      = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_"+sysname+"_post.root")
    f_ttbar_sys_odd_2  = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_p2_"+muOrEl+"_"+sysname+"_odd_post.root")

    if options.type == "full":
        response_sys = f_ttbar_sys.Get("response2_pt_split")
        response2_sys = f_ttbar_sys_2.Get("response2_pt_split")
    elif options.type == "half" :
        response_sys = f_ttbar_sys_odd.Get("response2_pt_split")
        response2_sys = f_ttbar_sys_odd_2.Get("response2_pt_split")
    elif options.type == "each" : #unfold p1 with p2
        response_sys = f_ttbar_sys_2.Get("response2_pt_split")
    else :
        response_sys = f_ttbar_sys_odd.Get("response2_pt_split")

    Hres_tmp = makeResponse(response_sys)
    if options.type == "full" or options.type == "half" :
        Hres_tmp_2 = makeResponse(response2_sys)
        Hres_tmp.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
        Hres_tmp_2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
        Hres_tmp.Add(Hres_tmp_2)
    if options.type == "half" or options.type == "original":
        Hres_tmp.Scale(2.0)
        
    Hres_sys[sysname] = Hres_tmp
        
# -------------------------------------------------------------------------------------
# output file for storing histograms to  
# -------------------------------------------------------------------------------------

fout = TFile("UnfoldingPlots/closureTest_TUnfold_pt_"+muOrEl+"_"+options.type+".root","recreate");

# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

if options.type == "full":
    if options.toy == "up" :
        thisMeas = f_ttbar.Get("ptRecoTopMod2_split").Clone()
        thisTrue = f_ttbar.Get("ptGenTopMod2").Clone()
        thisMeas2 = f_ttbar_2.Get("ptRecoTopMod2_split").Clone()
        thisTrue2 = f_ttbar_2.Get("ptGenTopMod2").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue = f_ttbar.Get("ptGenTopModDown2").Clone() 
        thisMeas2 = f_ttbar_2.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue2 = f_ttbar_2.Get("ptGenTopModDown2").Clone() 
    else : 
        thisMeas = f_ttbar.Get("ptRecoTop2_split").Clone() 
        thisTrue = f_ttbar.Get("ptGenTop2").Clone()
        thisMeas2 = f_ttbar_2.Get("ptRecoTop2_split").Clone() 
        thisTrue2 = f_ttbar_2.Get("ptGenTop2").Clone()
    removeFakes(thisMeas,response)
    removeFakes(thisMeas2,response2)
    thisMeas.Add(thisMeas2)
    thisTrue.Add(thisTrue2)
    thisMeas.Scale(lum * 831.76 / (77229341. + 78006311. * 1191. / 1192.));
    thisTrue.Scale(lum * 831.76 / (77229341. + 78006311. * 1191. / 1192.));

elif options.type == "half":
    if options.toy == "up" :
        thisMeas = f_ttbar_even.Get("ptRecoTopMod2_split").Clone()
        thisTrue = f_ttbar_even.Get("ptGenTopMod2").Clone()
        thisMeas2 = f_ttbar_even_2.Get("ptRecoTopMod2_split").Clone()
        thisTrue2 = f_ttbar_even_2.Get("ptGenTopMod2").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar_even.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue = f_ttbar_even.Get("ptGenTopModDown2").Clone() 
        thisMeas2 = f_ttbar_even_2.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue2 = f_ttbar_even_2.Get("ptGenTopModDown2").Clone() 
    else : 
        thisMeas = f_ttbar_even.Get("ptRecoTop2_split").Clone() 
        thisTrue = f_ttbar_even.Get("ptGenTop2").Clone() 
        thisMeas2 = f_ttbar_even_2.Get("ptRecoTop2_split").Clone() 
        thisTrue2 = f_ttbar_even_2.Get("ptGenTop2").Clone()
    removeFakes(thisMeas,response)
    removeFakes(thisMeas2,response2)
    thisMeas.Add(thisMeas2)
    thisTrue.Add(thisTrue2)
    thisMeas.Scale(lum * 831.76 * 2.0 / (77229341. + 78006311. * 1191. / 1192.));
    thisTrue.Scale(lum * 831.76 * 2.0 / (77229341. + 78006311. * 1191. / 1192.));
    
elif options.type == "each" : #unfold p1 with p2
    if options.toy == "up" :
        thisMeas = f_ttbar.Get("ptRecoTopMod2_split").Clone()
        thisTrue = f_ttbar.Get("ptGenTopMod2").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue = f_ttbar.Get("ptGenTopModDown2").Clone() 
    else : 
        thisMeas = f_ttbar.Get("ptRecoTop2_split").Clone() 
        thisTrue = f_ttbar.Get("ptGenTop2").Clone()
    removeFakes(thisMeas,response)
    thisMeas.Scale( lum * 831.76 / 77229341.)
    thisTrue.Scale( lum * 831.76 / 77229341.)

else :
    if options.toy == "up" :
        thisMeas = f_ttbar_even.Get("ptRecoTopMod2_split").Clone()
        thisTrue = f_ttbar_even.Get("ptGenTopMod2").Clone()
    elif options.toy == "dn" :
        thisMeas = f_ttbar_even.Get("ptRecoTopModDown2_split").Clone() 
        thisTrue = f_ttbar_even.Get("ptGenTopModDown2").Clone() 
    else : 
        thisMeas = f_ttbar_even.Get("ptRecoTop2_split").Clone() 
        thisTrue = f_ttbar_even.Get("ptGenTop2").Clone()
    removeFakes(thisMeas,response)
    thisMeas.Scale( lum * 831.76 * 2.0 / 77229341.)
    thisTrue.Scale( lum * 831.76 * 2.0 / 77229341.)

thisMeas.SetName("recolevel") 
thisTrue.SetName("truthlevel") 


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
#
# RUN PSEUDO EXPERIMENTS !!!
#
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Convert for TUnfold
# -------------------------------------------------------------------------------------

Hres = makeResponse(response)
if options.type == "half" or options.type == "full" :
    Hres_2 = makeResponse(response2)
    Hres.Scale(77229341. / (77229341. + 78006311. * 1191. / 1192.))
    Hres_2.Scale(78006311. * 1191. / 1192. / (77229341. + 78006311. * 1191. / 1192.))
    Hres.Add(Hres_2)
if options.type == "half" or options.type == "original":
    Hres.Scale(2.0)

# -------------------------------------------------------------------------------------
# generate toys
# -------------------------------------------------------------------------------------

hToy = thisMeas.Clone() 
hToy.SetName("toys")
hToy.Reset() 

hToy_i = [] 

ntoys = 1000

ntot = thisMeas.Integral()
nentries = int(thisMeas.GetEntries())
nbinsMeas = thisMeas.GetNbinsX()
nbinsTrue = thisTrue.GetNbinsX()

for itoy in xrange(0,ntoys) :
    hToy_tmp = hToy.Clone() 
    hToy_tmp.SetName("toys"+str(itoy))
     
    hToy_tmp.FillRandom(thisMeas,nentries)
    hToy_tmp.Scale(ntot/hToy_tmp.Integral())
    for ibin in xrange(1,thisMeas.GetNbinsX()+1):
        hToy_tmp.SetBinError(ibin,thisMeas.GetBinError(ibin))

    hToy_i.append(hToy_tmp)


# -------------------------------------------------------------------------------------
# UNFOLDING FOR TOYS
# -------------------------------------------------------------------------------------

hDiff_bin = []
h_tau = TH1F("tau_toys",";log(Tau);Toys",100,-4.0,-1.0)

for ibin in xrange(0,nbinsTrue) :
    hDiff_tmp = TH1F("diff_bin"+str(ibin),";(truth-unfolded)/truth; number of events",100,-1,1)
    hDiff_bin.append(hDiff_tmp)
        
for itoy in xrange(0,ntoys) :
    unfold_tmp = TUnfoldDensity(Hres,TUnfold.kHistMapOutputVert, TUnfold.kRegModeDerivative, TUnfold.kEConstraintNone, TUnfoldDensity.kDensityModeBinWidth)
    unfold_tmp.SetInput(hToy_i[itoy])
    for sysname in sysnames:
        unfold_tmp.AddSysError(Hres_sys[sysname],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)
    
    if options.regType == "LCurve" :
        logTauX_tmp = TSpline3()
        logTauY_tmp = TSpline3()
        lCurve_tmp = TGraph()
        iBest_tmp = unfold_tmp.ScanLcurve(30,0.,0.,lCurve_tmp,logTauX_tmp,logTauY_tmp)
        logTau_tmp = Double(0)
        x_tmp = Double(0)
        logTauX_tmp.GetKnot(iBest_tmp,logTau_tmp,x_tmp)
        h_tau.Fill(logTau_tmp)
    elif options.regType == "ScanTau" :
        scanResult_tmp = TSpline3()
        iBest_tmp = unfold_tmp.ScanTau(100,0.0001,0.1,scanResult_tmp,TUnfoldDensity.kEScanTauRhoAvgSys) #CAN CHANGE
        Tau_tmp = Double(0)
        rho_tmp = Double(0)
        scanResult_tmp.GetKnot(iBest_tmp,Tau_tmp,rho_tmp)
        h_tau.Fill(Tau_tmp)
    else :
        unfold_tmp.DoUnfold(0)
    
    hReco_tmp = unfold_tmp.GetOutput("tmp_output")

    for ibin in range(0, nbinsTrue):
        hDiff_bin[ibin].Fill( (thisTrue.GetBinContent(ibin+1) - hReco_tmp.GetBinContent(ibin+1)) / thisTrue.GetBinContent(ibin+1))

for ibin in xrange(0,nbinsTrue) :
    c = TCanvas()
    hDiff_bin[ibin].Draw("same")
    c.SaveAs("UnfoldingPlots/pull"+options.toy+"_"+options.regType+"_"+options.lepType+"_bin"+str(ibin)+"_"+options.type+".pdf")

if options.regType is not "None" :
    c2 = TCanvas()
    h_tau.Draw()
    c2.SaveAs("UnfoldingPlots/TauFromToys_"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type+".pdf")

hBias_pt = thisTrue.Clone() 
hBias_pt.SetName("biasVsPt_"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type)
hBias_pt.Reset()

for ibin in xrange(0,nbinsTrue) :
    mean = hDiff_bin[ibin].GetMean()
    err = hDiff_bin[ibin].GetRMS()

    hBias_pt.SetBinContent(ibin+1, mean)
    hBias_pt.SetBinError(ibin+1, err)
            
ccc = TCanvas()
gPad.SetGridy()
hBias_pt.GetYaxis().SetTitle("Bias")
hBias_pt.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
hBias_pt.SetAxisRange(400,1199,"X")
hBias_pt.SetAxisRange(-0.5,0.5,"Y")
hBias_pt.Draw()
ccc.SaveAs("UnfoldingPlots/bias_vspt"+options.toy+"_"+options.regType+"_"+options.lepType+"_"+options.type+".pdf")

hBias_pt.Write()

# -------------------------------------------------------------------------------------
# Done with toys, doing actual unfolding
# -------------------------------------------------------------------------------------

unfold = TUnfoldDensity(Hres,TUnfold.kHistMapOutputVert, TUnfold.kRegModeDerivative, TUnfold.kEConstraintNone, TUnfoldDensity.kDensityModeBinWidth)
unfold.SetInput(thisMeas)
for sysname in sysnames:
    unfold.AddSysError(Hres_sys[sysname],sysname,TUnfold.kHistMapOutputVert,TUnfoldDensity.kSysErrModeMatrix)

if options.regType == "LCurve" :
    logTauX = TSpline3()
    logTauY = TSpline3()
    lCurve = TGraph()
    bestLCurve = TGraph(1)
    iBest = unfold.ScanLcurve(30,0.,0.,lCurve,logTauX,logTauY)
    Tau = Double(0)
    x = Double(0)
    y = Double(0)
    logTauX.GetKnot(iBest,Tau,x)
    logTauY.GetKnot(iBest,Tau,y)
    bestLCurve.SetPoint(1,x,y)
    bestLCurve.SetMarkerColor(2)

elif options.regType == "ScanTau" :
    bestTau = TGraph(1)
    scanResult = TSpline3()
    iBest = unfold.ScanTau(100,0.0001,0.1,scanResult,TUnfoldDensity.kEScanTauRhoAvgSys)
    Tau = Double(0)
    rho = Double(0)
    scanResult.GetKnot(iBest,Tau,rho)
    bestTau.SetPoint(1,Tau,rho)
    bestTau.SetMarkerColor(2)

else :
    unfold.DoUnfold(0)

print "chi**2=" + str(unfold.GetChi2A()) + "+" + str(unfold.GetChi2L()) + " / " + str(unfold.GetNdf())

# unfolded distribution (histogram)
thisReco = unfold.GetOutput("reco")
        
# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------
# TODO: should fix BR

thisTrue.Scale(1.0/(lum*0.438/3.)) # true @ parton level
thisMeas.Scale(1.0/(lum*0.438/3.)) # measured @ reco level
thisReco.Scale(1.0/(lum*0.438/3.)) # unfolded to parton level

# -------------------------------------------------------------------------------------
# Adjust for bin width
# -------------------------------------------------------------------------------------

for ibin in range(1, nbinsTrue+1 ) :

    width = thisTrue.GetBinWidth(ibin)
    
    thisTrue.SetBinContent(ibin, thisTrue.GetBinContent(ibin) / width )
    thisTrue.SetBinError(ibin, thisTrue.GetBinError(ibin) / width )
    
    thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
    thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )
    
for ibin in range(1, nbinsMeas+1) :
        
    width = thisMeas.GetBinWidth(ibin)
    
    thisMeas.SetBinContent(ibin,  thisMeas.GetBinContent(ibin) / width )
    thisMeas.SetBinError(ibin,  thisMeas.GetBinError(ibin) / width )

# -------------------------------------------------------------------------------------
# draw parton-level unfolding
# -------------------------------------------------------------------------------------

## ratio of unfolded data to generator-level

hFrac = thisReco.Clone()
hFrac.SetName("hFrac")
hFrac.SetTitle(";Top quark p_{T} (GeV);Data/MC")
hFrac.Divide(thisTrue)

c1 = TCanvas("c", "c", 700, 700)
pad1 =  TPad("pad1","pad1",0,0.3,1,1)
pad1.SetBottomMargin(0.05);
pad1.Draw();
pad1.cd();

thisReco.SetMarkerStyle(21)
thisMeas.SetMarkerStyle(25);

thisReco.GetXaxis().SetRangeUser(400.,1199.)
thisTrue.GetXaxis().SetRangeUser(400.,1199.)
thisMeas.GetXaxis().SetRangeUser(400.,1199.)

xsec_title = ";;d#sigma/dp_{T} [fb/GeV]"

thisReco.SetTitle(xsec_title)
thisReco.GetYaxis().SetTitleOffset(1.2)
thisReco.SetMinimum(0.0)
max = thisTrue.GetMaximum()
max2 = thisReco.GetMaximum()
if max2 > max:
	max = max2
thisReco.SetAxisRange(0,max*1.15,"Y")
thisReco.Draw()
thisTrue.Draw('hist same')
thisMeas.Draw('same')
thisTrue.UseCurrentStyle()
thisTrue.SetLineColor(4);
thisTrue.GetYaxis().SetTitleSize(25)
thisTrue.GetXaxis().SetLabelSize(0)

leg = TLegend(0.5, 0.55, 0.9, 0.75)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)

tt = TLatex()
tt.SetNDC()
tt.SetTextFont(42)
leg.AddEntry( thisReco, 'Unfolded MC (Powheg)', 'p')
leg.AddEntry( thisTrue, 'Generated (Powheg)', 'l')
leg.AddEntry( thisMeas, 'Reco-level (Powheg)', 'p')
tt.DrawLatex(0.55,0.45, "MC closure test")
leg.Draw()

# write histograms to file
thisReco.SetName("UnfoldedMC")

thisReco.Write()
thisTrue.Write()
thisMeas.Write()

text1 = TLatex()
text1.SetNDC()
text1.SetTextFont(42)
text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 35.9 fb^{-1}, #sqrt{s} = 13 TeV}")

c1.cd();
pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
pad2.SetTopMargin(0.05);
pad2.SetBottomMargin(0.4);
pad2.Draw();
pad2.cd();
pad2.SetGridy()
hFrac.SetMaximum(1.8)
hFrac.SetMinimum(0.2)
hFrac.UseCurrentStyle()
hFrac.GetYaxis().SetTitleSize(25)
hFrac.GetYaxis().SetTitleOffset(2.0)
hFrac.GetXaxis().SetTitleOffset(4.0)
hFrac.GetXaxis().SetLabelSize(25)
hFrac.GetYaxis().SetNdivisions(4,4,0,False)

hFrac.Draw("e")
hFrac.GetXaxis().SetRangeUser(400., 1199.)

c1.Update()

c1.SaveAs("UnfoldingPlots/closure_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+"_result.pdf")

# -------------------------------------------------------------------------------------
# Plot L-curve scan
# -------------------------------------------------------------------------------------

if not options.regType == "None" :
    c2 = TCanvas("c2", "c2", 700, 700)
    c2.cd()
    if options.regType == "LCurve" :
        lCurve.GetXaxis().SetTitle("log L_{1}")
        lCurve.GetYaxis().SetTitle("log L_{2} / #tau^{2}")
        lCurve.Draw()
        bestLCurve.Draw("*")
    else :
        #dummy = TGraph()
        #dummy.GetXaxis.SetRangeUser(-4.0,-1.0)
        #dummy.GetXaxis.SetTitle("log(#tau)")
        #dummy.GetYaxis.SetRangeUser(0.5,1.0)
        #dummy.GetYaxis.SetTitle("Correlation")
        #dummy.Draw()
        scanResult.Draw("P")
        bestTau.Draw("*")
        
    tl2 = TLatex()
    tl2.SetNDC()
    tl2.SetTextFont(42)
    legend = "log(#tau) = %.3e" % Tau
    tl2.DrawLatex(0.55,0.8,legend)

    c2.SaveAs("UnfoldingPlots/TauScan_"+options.regType+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

fout.Close()
