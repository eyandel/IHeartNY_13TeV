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

parser.add_option('--nbr', metavar='F', type='string', action='store',
                  default='5',
                  dest='nbr',
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

from ROOT import gRandom, TH1, TH1D, TH1F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD, TGraph, TUnfold, Double, TSpline, TSpline3, TUnfoldDensity
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
# cross sections, efficiencies, total number of events
# -------------------------------------------------------------------------------------

lum = 12337.98 # was 12358.75 (pb-1), but missing some data
if options.lepType == "ele":
    lum = 12267.67 # was 12295.65 (pb-1), but missing some data

PowhegPythia8_norm = 831.76 * lum / 182123200.                                                                                    

eff_closure = 2.0
if options.type == "full":
    eff_closure = 1.0

    
# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"


# In the below, file named f_..._odd will be the one from which response matrix is extracted from

if options.type == "full":
    f_ttbar     = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_post.root")
else :
    f_ttbar     = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_even_post.root")
    f_ttbar_odd = TFile("histfiles_80X/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_odd_post.root")

    
# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

if options.type == "full":
    response = f_ttbar.Get("response_pt"+options.nbr+"fine")
else :
    response = f_ttbar_odd.Get("response_pt"+options.nbr+"fine")
TH1.AddDirectory(0)


# -------------------------------------------------------------------------------------
# output file for storing histograms to  
# -------------------------------------------------------------------------------------


fout = TFile("UnfoldingPlots/closureTest"+options.nbr+"_"+options.type+".root","recreate");


# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

hTrue = f_ttbar.Get("ptGenTop"+options.nbr)
hTrueUp = f_ttbar.Get("ptGenTopMod"+options.nbr)
hTrueDn = f_ttbar.Get("ptGenTopModDown"+options.nbr)

hTrue.Scale(PowhegPythia8_norm * eff_closure)
hTrue.SetName("ptGenTop"+options.nbr+"_true")

hTrueUp.Scale(PowhegPythia8_norm * eff_closure)
hTrueDn.Scale(PowhegPythia8_norm * eff_closure)
hTrueUp.SetName("ptGenTopUp"+options.nbr+"_true")
hTrueDn.SetName("ptGenTopDn"+options.nbr+"_true")


hMeas = f_ttbar.Get("ptRecoTop"+options.nbr+"fine")
hMeasUp = f_ttbar.Get("ptRecoTopMod"+options.nbr+"fine")
hMeasDn = f_ttbar.Get("ptRecoTopModDown"+options.nbr+"fine")

hMeas.Scale(PowhegPythia8_norm * eff_closure)
hMeas.SetName("ptRecoTop"+options.nbr+"fine"+"_measured")

hMeasUp.Scale(PowhegPythia8_norm * eff_closure)
hMeasDn.Scale(PowhegPythia8_norm * eff_closure)
hMeasUp.SetName("ptRecoTopUp"+options.nbr+"fine"+"_measured")
hMeasDn.SetName("ptRecoTopDn"+options.nbr+"fine"+"_measured")

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
#
# RUN PSEUDO EXPERIMENTS !!!
#
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# do this for nominal / UP / DN 

if options.toy == "up" :
    thisMeas = hMeasUp.Clone() 
    thisTrue = hTrueUp.Clone() 
elif options.toy == "dn" :
    thisMeas = hMeasDn.Clone() 
    thisTrue = hTrueDn.Clone() 
else : 
    thisMeas = hMeas.Clone() 
    thisTrue = hTrue.Clone() 

thisMeas.SetName("recolevel") 
thisTrue.SetName("truthlevel") 

# -------------------------------------------------------------------------------------
# Convert for TUnfold
# -------------------------------------------------------------------------------------

Hres = response.HresponseNoOverflow()

tru = response.Vtruth()
for j in xrange(1,response.GetNbinsTruth()) :
    ntru = 0.0
    for i in xrange(1,response.GetNbinsMeasured()) :
        ntru += Hres.GetBinContent(i,j)
    Hres.SetBinContent(response.GetNbinsMeasured()+1, j, tru[j-1]-ntru)
        
if response.FakeEntries() :
    fakes = response.Vfakes()
    fac = response.Vmeasured().Sum()
    if fac != 0.0 :
        measVec = RooUnfoldRespone.H2V(thisMeas,response.GetNbinsMeasured(),0)
        fac = measVec.Sum() / fac
    for i in xrange(1,response.GetNbinsMeasured()) :
        thisMeas.SetBinContent(i,thisMeas.GetBinContent(i) - fac * fakes[i-1])

# -------------------------------------------------------------------------------------
# generate toys
# -------------------------------------------------------------------------------------

hToy = thisMeas.Clone() 
hToy.SetName("toys")
hToy.Reset() 

hToy_i = [] 

ntoys = 10000

ntot = thisMeas.GetSum()
nbinsMeas = thisMeas.GetNbinsX()
nbinsTrue = thisTrue.GetNbinsX()

for itoy in xrange(0,ntoys) :
    hToy_tmp = hToy.Clone() 
    hToy_tmp.SetName("toys"+str(itoy))
     
    hToy_tmp.FillRandom(thisMeas,10000)
    hToy_tmp.Scale(ntot/hToy_tmp.GetSum())

    hToy_i.append(hToy_tmp)


# -------------------------------------------------------------------------------------
# UNFOLDING FOR TOYS
# -------------------------------------------------------------------------------------

hDiff_bin = [] # for each iterations, list of bins (i.e. [1st it, bin0], [1st it, bin 1], ... , [1st it, bin n], [2nd it, bin0] , ...  

for ibin in xrange(0,nbinsTrue) :
        hDiff_tmp = TH1F("diff_bin"+str(ibin),"; unfolded-truth; number of events",100,-1,1)
        hDiff_bin.append(hDiff_tmp)
        
lowedge = 399.
highedge = 1199.

for itoy in xrange(0,ntoys) :

    unfold_tmp = TUnfoldDensity(Hres,TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintNone, TUnfoldDensity.kDensityModeBinWidth)
    unfold_tmp.SetInput(hToy_i[itoy])
    if options.regType == "LCurve" :
        logTauX_tmp = TSpline3()
        logTauY_tmp = TSpline3()
        lCurve_tmp = TGraph()
        unfold_tmp.ScanLcurve(30,0.,0.,lCurve_tmp,logTauX_tmp,logTauY_tmp)
    else :
        scanResult = TSpline3()
        unfold_tmp.ScanTau(100,0.0001,0.1,scanResult,TUnfoldDensity.kEScanTauRhoAvg)
    
    hReco_tmp = unfold_tmp.GetOutput("tmp_output")

    for ibin in range(0, nbinsTrue):
        if thisTrue.GetBinLowEdge(ibin+1) > lowedge and thisTrue.GetBinLowEdge(ibin+1) < highedge:
            hDiff_bin[ibin].Fill( (thisTrue.GetBinContent(ibin+1) - hReco_tmp.GetBinContent(ibin+1)) / thisTrue.GetBinContent(ibin+1))

for ibin in xrange(0,nbinsTrue) :
    c = TCanvas()
    hDiff_bin[ibin].Draw("same")
    c.SaveAs("pull"+options.toy+"_"+options.lepType+"_bin"+str(ibin)+".pdf")


# sum relative bias squared and standard deviation squared (sqrt(bias^2)_{sum over bins} / N_bins) and plot vs # iterations 
hBias_pt = thisTrue.Clone() 
hBias_pt.SetName("bias_vs_pt")
hBias_pt.Reset()

for ibin in xrange(0,nbinsTrue) :
    mean = hDiff_bin[ibin].GetMean()
    err = hDiff_bin[ibin].GetRMS()

    hBias_pt.SetBinContent(ibin+1, mean)
    hBias_pt.SetBinError(ibin+1, err)
            
ccc = TCanvas()
gPad.SetGridy()
for bin in xrange(0,nbinsTrue):
    print 'bin ', bin, hBias_pt.GetBinContent(bin+1)
hBias_pt.GetYaxis().SetTitle("Bias")
hBias_pt.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
hBias_pt.SetAxisRange(410,1100,"X")
hBias_pt.SetAxisRange(-0.5,0.5,"Y")
hBias_pt.Draw()
ccc.SaveAs("bias_vspt"+options.toy+"_"+options.lepType+".pdf")

# -------------------------------------------------------------------------------------
# do the unfolding for different nbr of iterations 
# -------------------------------------------------------------------------------------

unfold = TUnfoldDensity(Hres,TUnfold.kHistMapOutputVert, TUnfold.kRegModeCurvature, TUnfold.kEConstraintNone, TUnfoldDensity.kDensityModeBinWidth)
unfold.SetInput(thisMeas)

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

else :
    bestTau = TGraph(1)
    scanResult = TSpline3()
    iBest = unfold.ScanTau(100,0.0001,0.1,scanResult,TUnfoldDensity.kEScanTauRhoAvg)
    Tau = Double(0)
    rho = Double(0)
    scanResult.GetKnot(iBest,Tau,rho)
    bestTau.SetPoint(1,Tau,rho)
    bestTau.SetMarkerColor(2)

print "chi**2=" + str(unfold.GetChi2A()) + "+" + str(unfold.GetChi2L()) + " / " + str(unfold.GetNdf())

# unfolded distribution (histogram)
thisReco = unfold.GetOutput("reco")

lowedge = 399.
highedge = 1199.

for ibin in range(0, nbinsTrue):

    if thisTrue.GetBinLowEdge(ibin+1) < lowedge or thisTrue.GetBinLowEdge(ibin+1) > highedge:
        
        thisTrue.SetBinContent(ibin+1,0)
        thisTrue.SetBinError(ibin+1,0)
        
        thisReco.SetBinContent(ibin+1,0)
        thisReco.SetBinError(ibin+1,0)
        
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

    if thisTrue.GetBinLowEdge(ibin) < lowedge or thisTrue.GetBinLowEdge(ibin) > highedge :
        continue
    
    width = thisTrue.GetBinWidth(ibin)

    print thisTrue.GetBinLowEdge(ibin)
    
    thisTrue.SetBinContent(ibin, thisTrue.GetBinContent(ibin) / width )
    thisTrue.SetBinError(ibin, thisTrue.GetBinError(ibin) / width )

    thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
    thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )

for ibin in range(1, nbinsMeas+1) :
    if thisMeas.GetBinLowEdge(ibin) < lowedge or thisMeas.GetBinLowEdge(ibin) > highedge :
        continue
    
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

thisReco.GetXaxis().SetRangeUser(400.,1200.)
thisTrue.GetXaxis().SetRangeUser(400.,1200.)
thisMeas.GetXaxis().SetRangeUser(400.,1200.)

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
text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 12.4 fb^{-1}, #sqrt{s} = 13 TeV}")

c1.cd();
pad2 =  TPad("pad2","pad2",0,0.0,1,0.28)
pad2.SetTopMargin(0.05);
pad2.SetBottomMargin(0.4);
pad2.Draw();
pad2.cd();
pad2.SetGridy()
hFrac.SetMaximum(1.4)
hFrac.SetMinimum(0.6)
hFrac.UseCurrentStyle()
hFrac.GetYaxis().SetTitleSize(25)
hFrac.GetYaxis().SetTitleOffset(2.0)
hFrac.GetXaxis().SetTitleOffset(4.0)
hFrac.GetXaxis().SetLabelSize(25)
hFrac.GetYaxis().SetNdivisions(4,4,0,False)

hFrac.Draw("e")
hFrac.GetXaxis().SetRangeUser(400., 1200.)

c1.Update()

c1.SaveAs("UnfoldingPlots/closure"+options.nbr+"_"+options.lepType+"_"+options.toy+"_"+options.type+"_result.pdf")

# -------------------------------------------------------------------------------------
# Plot L-curve scan
# -------------------------------------------------------------------------------------

c2 = TCanvas("c2", "c2", 700, 700)
c2.cd()
if options.regType == "LCurve" :
    lCurve.Draw()
    bestLCurve.Draw("*")
else :
    scanResult.Draw("P")
    bestTau.Draw("*")

tl2 = TLatex()
tl2.SetNDC()
tl2.SetTextFont(42)
legend = "#tau = %.3e" % Tau
tl2.DrawLatex(0.55,0.8,legend)

c2.SaveAs("TauScan"+options.nbr+"_"+options.lepType+"_"+options.toy+"_"+options.type+".pdf")

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

gStyle.SetPadRightMargin(0.12);
cr = TCanvas("c_response", "", 800, 600)

hEmpty2D = Hres.Clone()
hEmpty2D.SetName("empty2D")
hEmpty2D.Reset()
hEmpty2D.GetXaxis().SetTitle("Reconstructed top jet p_{T} (GeV)")
hEmpty2D.GetYaxis().SetTitle("Top quark p_{T} (GeV)")
hEmpty2D.GetXaxis().SetLabelSize(0.045)
hEmpty2D.GetYaxis().SetLabelSize(0.045)
hEmpty2D.Draw()

hResponse2D = Hres.Clone()
hResponse2D.SetName("plottedResponse")

gStyle.SetPaintTextFormat(".1f")
hResponse2D.Draw("colz,same,text")
hEmpty2D.Draw("axis,same")


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

cr.SaveAs("UnfoldingPlots/closure"+options.nbr+"_"+options.lepType+"_"+options.type+"_response.pdf")

hEmpty2D.SetAxisRange(410,1150,"X")
hEmpty2D.SetAxisRange(410,1150,"Y")
hResponse2D.SetAxisRange(410,1150,"X")
hResponse2D.SetAxisRange(410,1150,"Y")
hEmpty2D.Draw()
hResponse2D.Draw("colz,same,text")
hEmpty2D.Draw("axis,same")

cr.SaveAs("UnfoldingPlots/closure"+options.nbr+"_"+options.lepType+"_"+options.type+"_response_zoom.pdf")
        
Hres.SetName("responseMatrix"+options.nbr+"_"+options.lepType+"_"+options.type)
Hres.Write()
    

fout.Close()
