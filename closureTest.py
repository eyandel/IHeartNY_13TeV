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


# -------------------------------------------------------------------------------------
# load options & set plot style
# -------------------------------------------------------------------------------------

(options, args) = parser.parse_args()
argv = []

import sys

from ROOT import gRandom, TH1, TH1D, TH1F, cout, TFile, gSystem, TCanvas, TPad, gPad, gROOT, gStyle, THStack, TLegend, TLatex, TColor, TMath, TVectorD
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

# -------------------------------------------------------------------------------------
# cross sections, efficiencies, total number of events
# -------------------------------------------------------------------------------------

lum = 35867.0
PowhegPythia8_norm = 831.76 * lum / 77229341.

eff_closure = 2.0
if options.type == "full":
    eff_closure = 1.0

best_iter = 10;
    
# -------------------------------------------------------------------------------------
#  read histogram files
# -------------------------------------------------------------------------------------

muOrEl = "mu"
if options.lepType=="ele":
    muOrEl = "el"

# In the below, file named f_..._odd will be the one from which response matrix is extracted from

if options.type == "full":
    f_ttbar     = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_post.root")
else :
    f_ttbar     = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_even_post.root")
    f_ttbar_odd = TFile("histfiles_full2016/hists_PowhegPythia8_fullTruth_"+muOrEl+"_nom_odd_post.root")
    
# -------------------------------------------------------------------------------------
# Get response matrix
# -------------------------------------------------------------------------------------

if options.type == "full":
    response = f_ttbar.Get("response_pt")
else :
    response = f_ttbar_odd.Get("response_pt")
response.UseOverflow()
TH1.AddDirectory(0)

# -------------------------------------------------------------------------------------
# output file for storing histograms to  
# -------------------------------------------------------------------------------------

fout = TFile("UnfoldingPlots/closureTest_dAgostini_pt_"+muOrEl+"_"+options.type+".root","recreate");

# -------------------------------------------------------------------------------------
# read & normalize histograms
# -------------------------------------------------------------------------------------

hTrue = f_ttbar.Get("ptGenTop")
hTrueUp = f_ttbar.Get("ptGenTopMod")
hTrueDn = f_ttbar.Get("ptGenTopModDown")

hTrue.Scale(PowhegPythia8_norm * eff_closure)
hTrueUp.Scale(PowhegPythia8_norm * eff_closure)
hTrueDn.Scale(PowhegPythia8_norm * eff_closure)

hMeas = f_ttbar.Get("ptRecoTop")
hMeasUp = f_ttbar.Get("ptRecoTopMod")
hMeasDn = f_ttbar.Get("ptRecoTopModDown")

hMeas.Scale(PowhegPythia8_norm * eff_closure)
hMeasUp.Scale(PowhegPythia8_norm * eff_closure)
hMeasDn.Scale(PowhegPythia8_norm * eff_closure)

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
#
# RUN PSEUDO EXPERIMENTS !!!
#
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# do this for nominal / UP / DN 
if options.toy == "up" :
    myReco = hMeasUp.Clone() 
    myTrue = hTrueUp.Clone() 
elif options.toy == "dn" :
    myReco = hMeasDn.Clone() 
    myTrue = hTrueDn.Clone() 
else : 
    myReco = hMeas.Clone() 
    myTrue = hTrue.Clone() 

myReco.SetName("recolevel") 
myTrue.SetName("truthlevel") 

# -------------------------------------------------------------------------------------
# generate toys
# -------------------------------------------------------------------------------------

hToy = myReco.Clone() 
hToy.SetName("toys")
hToy.Reset() 

hToy_i = [] 

ntoys = 10000

ntot = myReco.GetSum()
nbins = myReco.GetNbinsX()

for itoy in xrange(0,ntoys) :
    hToy_tmp = hToy.Clone() 
    hToy_tmp.SetName("toys"+str(itoy))
     
    hToy_tmp.FillRandom(myReco,10000)
    hToy_tmp.Scale(ntot/hToy_tmp.GetSum())

    hToy_i.append(hToy_tmp)

# -------------------------------------------------------------------------------------
# UNFOLDING FOR TOYS
# -------------------------------------------------------------------------------------

tot_iter = 15
n_iter = array('i', [])

hDiff_bin = [] #Array of pull hists for each bin and iteration
               #Order is [1st it, bin0], [1st it, bin 1], ... , [1st it, bin n], [2nd it, bin0] , ...  

for iiter in xrange(1,tot_iter+1) :
    for ibin in xrange(0,nbins) :
        hDiff_tmp = TH1F("diff_bin"+str(ibin)+"_"+str(iiter),";(truth-unfolded)/truth; number of events",100,-0.5,0.5)
        hDiff_bin.append(hDiff_tmp)
        
for iiter in xrange(1,tot_iter+1) :
    for itoy in xrange(0,ntoys) :

        unfold = RooUnfoldBayes(response, hToy_i[itoy], iiter)
        hReco_tmp = unfold.Hreco()

        for ibin in range(0, nbins):
            hDiff_bin[(iiter-1)*nbins + ibin].Fill((myTrue.GetBinContent(ibin+1) - hReco_tmp.GetBinContent(ibin+1)) / myTrue.GetBinContent(ibin+1))                

# Plotting
color = [1,2,3,4,5,6,7,8,9,14]

for ibin in xrange(0,nbins) :
    c = TCanvas()
    if hDiff_bin[ibin].GetSum() == 0: 
        continue
    ll = TLegend(0.75, 0.5, 0.9, 0.9)
    ll.SetFillStyle(0)
    ll.SetTextFont(42)
    ll.SetTextSize(0.035)
    ll.SetBorderSize(0)

    for iiter in xrange(0,tot_iter) :
        if iiter > 9: 
            continue
        hDiff_bin[iiter*nbins + ibin].SetLineColor(color[iiter])
        hDiff_bin[iiter*nbins + ibin].Draw("same")
        ll.AddEntry(hDiff_bin[iiter*nbins + ibin],'Iter '+str(iiter+1), 'l')

    ll.Draw()    
    c.SaveAs("UnfoldingPlots/pull"+options.toy+"_"+options.lepType+"_bin"+str(ibin)+"_"+options.type+".pdf")


# sum relative bias squared and standard deviation squared (sqrt(bias^2)_{sum over bins} / N_bins) and plot vs # iterations 

bias = TH1F("bias",";# iterations; ", tot_iter,0,tot_iter)
std  = TH1F("std", ";# iterations; ", tot_iter,0,tot_iter)
biasstd = TH1F("biasstd",";# iterations; ", tot_iter,0,tot_iter)

hBias_pt = myReco.Clone() 
hBias_pt.SetName("bias_vs_pt")
hBias_pt.Reset()
print 'hBias_pt has ' + str(hBias_pt.GetNbinsX()) + ' bins'

for iiter in xrange(0,tot_iter) :

    sum_bias = 0
    sum_std = 0
    sum_biasstd = 0

    filledbins = 0

    for ibin in xrange(0,nbins-1) : #Skipping last bin
        if hDiff_bin[iiter*nbins + ibin].GetSum() == 0:
            continue
        filledbins += 1
        mean = hDiff_bin[iiter*nbins + ibin].GetMean()
        err = hDiff_bin[iiter*nbins + ibin].GetRMS()
        sum_bias += mean*mean
        sum_std += err*err
        sum_biasstd += (mean*mean + err*err)

        if iiter+1 == best_iter:
            hBias_pt.SetBinContent(ibin+1, mean)
            hBias_pt.SetBinError(ibin+1, err)
            
    sum_bias = math.sqrt(sum_bias/filledbins)
    sum_std = math.sqrt(sum_std/filledbins)
    sum_biasstd = math.sqrt(sum_biasstd/filledbins)
    
    bias.SetBinContent(iiter+1, sum_bias)
    std.SetBinContent(iiter+1, sum_std)
    biasstd.SetBinContent(iiter+1, sum_biasstd)    

cc = TCanvas()
biasstd.SetMinimum(0)
biasstd.SetMaximum(biasstd.GetMaximum()*1.5)
biasstd.Draw()
bias.SetLineColor(4)
bias.Draw("same")
std.SetLineColor(2)
std.Draw("same")

ll = TLegend(0.3, 0.7, 0.6, 0.9)
ll.SetFillStyle(0)
ll.SetTextFont(42)
ll.SetTextSize(0.035)
ll.SetBorderSize(0)
tt = TLatex()
tt.SetNDC()
tt.SetTextFont(42)
ll.AddEntry( std, '#sqrt{#sum std^{2}/N_{bin}}', 'l')
ll.AddEntry( bias, '#sqrt{#sum bias^{2}/N_{bin}}', 'l')
ll.AddEntry( biasstd, '#sqrt{#sum (std^{2} + bias^{2})/N_{bin}}', 'l')
ll.Draw()
cc.SaveAs("UnfoldingPlots/biasstd"+options.toy+"_"+options.lepType+"_"+options.type+".pdf")

ccc = TCanvas()
gPad.SetGridy()
hBias_pt.GetYaxis().SetTitle("Bias")
hBias_pt.GetXaxis().SetTitle("Top quark p_{T} (GeV)")
hBias_pt.SetAxisRange(400.,1199.,"X")
hBias_pt.SetAxisRange(-0.5,0.5,"Y")
hBias_pt.Draw()
ccc.SaveAs("UnfoldingPlots/bias_vspt"+options.toy+"_"+options.lepType+"_"+options.type+".pdf")

# -------------------------------------------------------------------------------------
# Done with toys, now doing single unfolding
# -------------------------------------------------------------------------------------

n_iter = array('i', [])
pvalues = array('f', [])
pvaluesUp = array('f', [])
pvaluesDn = array('f', [])

hErrorRatio = []
hErrorRatioUp = []
hErrorRatioDn = []

hReco_iter = []
hRecoUp_iter = []
hRecoDn_iter = []

for iiter in xrange(0,tot_iter) :
    hError = hMeas.Clone() 
    hError.SetName("ErrorRatio_iter"+str(iiter))
    hError.Reset()
    hErrorRatio.append(hError)

    hErrorUp = hMeas.Clone() 
    hErrorUp.SetName("ErrorRatioUp_iter"+str(iiter))
    hErrorUp.Reset()
    hErrorRatioUp.append(hErrorUp)

    hErrorDn = hMeas.Clone() 
    hErrorDn.SetName("ErrorRatioDn_iter"+str(iiter))
    hErrorDn.Reset()
    hErrorRatioDn.append(hErrorDn)

nbins = hMeas.GetNbinsX()
print 'hMeas has ' + str(nbins) + ' bins'

for i_iter in xrange(1,tot_iter+1) :

    unfold   = RooUnfoldBayes(response, hMeas, i_iter)
    unfoldUp = RooUnfoldBayes(response, hMeasUp, i_iter)
    unfoldDn = RooUnfoldBayes(response, hMeasDn, i_iter)

    # unfolded distribution (histogram)
    hReco = unfold.Hreco()
    hRecoUp = unfoldUp.Hreco()
    hRecoDn = unfoldDn.Hreco()

    hReco_iter.append(hReco)
    hRecoUp_iter.append(hRecoUp)
    hRecoDn_iter.append(hRecoDn)
    
    # error ratios
    for ibin in range(0, nbins):

        rel_err_reco = 0
        rel_err_unfold = 0
        rel_err_recoUp = 0
        rel_err_unfoldUp = 0
        rel_err_recoDn = 0
        rel_err_unfoldDn = 0
        
        if hMeas.GetBinContent(ibin+1) > 0:
            rel_err_reco = hMeas.GetBinError(ibin+1) / hMeas.GetBinContent(ibin+1)
        if hReco.GetBinContent(ibin+1) > 0:
            rel_err_unfold = hReco.GetBinError(ibin+1) / hReco.GetBinContent(ibin+1)

        if hMeasUp.GetBinContent(ibin+1) > 0:
            rel_err_recoUp = hMeasUp.GetBinError(ibin+1) / hMeasUp.GetBinContent(ibin+1)
        if hRecoUp.GetBinContent(ibin+1) > 0:
            rel_err_unfoldUp = hRecoUp.GetBinError(ibin+1) / hRecoUp.GetBinContent(ibin+1)

        if hMeasDn.GetBinContent(ibin+1) > 0:
            rel_err_recoDn = hMeasDn.GetBinError(ibin+1) / hMeasDn.GetBinContent(ibin+1)
        if hRecoDn.GetBinContent(ibin+1) > 0:
            rel_err_unfoldDn = hRecoDn.GetBinError(ibin+1) / hRecoDn.GetBinContent(ibin+1)
        
        err_ratio = 0
        err_ratioUp = 0
        err_ratioDn = 0
        if rel_err_reco > 0:
            err_ratio = rel_err_unfold / rel_err_reco
        if rel_err_recoUp > 0:
            err_ratioUp = rel_err_unfoldUp / rel_err_recoUp
        if rel_err_recoDn > 0:
            err_ratioDn = rel_err_unfoldDn / rel_err_recoDn
            
        hErrorRatio[i_iter-1].SetBinContent(ibin+1, err_ratio)
        hErrorRatioUp[i_iter-1].SetBinContent(ibin+1, err_ratioUp)
        hErrorRatioDn[i_iter-1].SetBinContent(ibin+1, err_ratioDn)
    
    # unfolded & measured (vectors)
    vReco = unfold.Vreco()
    vMeas = unfold.Vmeasured()
    vFake = response.Vfakes()

    vRecoUp = unfoldUp.Vreco()
    vRecoDn = unfoldDn.Vreco()
    vMeasUp = unfoldUp.Vmeasured()
    vMeasDn = unfoldDn.Vmeasured()

    # response matrix as matrix instead of 2D histogram: (row,column)=(measured,truth)
    mResponse = response.Mresponse()

    # refolded (vector)
    vRefold = TVectorD(vReco)
    vRefold *= mResponse
    print 'vRefold has length ' + str(vRefold.GetNrows())

    vRefoldUp = TVectorD(vRecoUp)
    vRefoldDn = TVectorD(vRecoDn)
    vRefoldUp *= mResponse
    vRefoldDn *= mResponse

    hRefold = TH1D("hRefold", "refolded", nbins, 0, nbins)
    hMeasCheck = TH1D("hMeasCheck", "measured check", nbins, 0, nbins)

    hRefoldUp = TH1D("hRefoldUp", "refoldedUp", nbins, 0, nbins)
    hRefoldDn = TH1D("hRefoldDn", "refoldedDn", nbins, 0, nbins)
    hMeasCheckUp = TH1D("hMeasCheckUp", "measured check up", nbins, 0, nbins)
    hMeasCheckDn = TH1D("hMeasCheckDn", "measured check down", nbins, 0, nbins)

    for ibin in range(0, nbins):

        print 'vRefold['+str(ibin)+'] = ' + str(vRefold[ibin])
        hRefold.SetBinContent(ibin+1, vRefold[ibin])
        hRefoldUp.SetBinContent(ibin+1, vRefoldUp[ibin])
        hRefoldDn.SetBinContent(ibin+1, vRefoldDn[ibin])

        print 'vMeas['+str(ibin)+'] = ' + str(vMeas[ibin])
        print 'vMeas-Fake['+str(ibin)+'] = ' + str(vMeas[ibin]-vFake[ibin])
        hMeasCheck.SetBinContent(ibin+1, vMeas[ibin]-vFake[ibin])
        hMeasCheckUp.SetBinContent(ibin+1, vMeasUp[ibin]-vFake[ibin])
        hMeasCheckDn.SetBinContent(ibin+1, vMeasDn[ibin]-vFake[ibin])


    pvalue = hMeasCheck.Chi2Test(hRefold, "WW")
    pvalueUp = hMeasCheckUp.Chi2Test(hRefoldUp, "WW")
    pvalueDn = hMeasCheckDn.Chi2Test(hRefoldDn, "WW")

    n_iter.append(i_iter)
    pvalues.append(pvalue)
    pvaluesUp.append(pvalueUp)
    pvaluesDn.append(pvalueDn)

print ""
for i in xrange(0,tot_iter-1):
    print 'iter # ', n_iter[i], ' pvalue = ', pvalues[i], ' pvalueUP = ', pvaluesUp[i], ' pvalueDN = ', pvaluesDn[i]
print ""

cerr = TCanvas()
cleg = TLegend(0.5, 0.5, 0.9, 0.9)
cleg.SetFillStyle(0)
cleg.SetTextFont(42)
cleg.SetTextSize(0.045)
cleg.SetBorderSize(0)

for iiter in xrange(0,tot_iter) :
    if iiter > 9:
        continue
    hErrorRatio[iiter].SetAxisRange(0.5,2.0,"Y")
    hErrorRatio[iiter].SetAxisRange(400,1199,"X")
    hErrorRatio[iiter].GetYaxis().SetTitle("Error ratio")
    hErrorRatio[iiter].GetXaxis().SetTitle("Top quark p_{T} (GeV)")
    hErrorRatio[iiter].SetLineColor(color[iiter])
    hErrorRatio[iiter].Draw("same")
    cleg.AddEntry(hErrorRatio[iiter], 'Iteration '+str(iiter), 'l')
cleg.Draw()
cerr.SaveAs("UnfoldingPlots/errorRatio_"+options.lepType+"_"+options.type+".pdf")

cerrUp = TCanvas()
for iiter in xrange(0,tot_iter) :
    if iiter > 9:
        continue
    hErrorRatioUp[iiter].SetAxisRange(0.5,2.0,"Y")
    hErrorRatioUp[iiter].SetAxisRange(400,1199,"X")
    hErrorRatioUp[iiter].GetYaxis().SetTitle("Error ratio")
    hErrorRatioUp[iiter].GetXaxis().SetTitle("Top quark p_{T} (GeV)")
    hErrorRatioUp[iiter].SetLineColor(color[iiter])
    hErrorRatioUp[iiter].Draw("same")
cleg.Draw()
cerrUp.SaveAs("UnfoldingPlots/errorRatioUp_"+options.lepType+"_"+options.type+".pdf")

cerrDn = TCanvas()
for iiter in xrange(0,tot_iter) :
    if iiter > 9:
        continue
    hErrorRatioDn[iiter].SetAxisRange(0.5,2.0,"Y")
    hErrorRatioDn[iiter].SetAxisRange(400,1199,"X")
    hErrorRatioDn[iiter].GetYaxis().SetTitle("Error ratio")
    hErrorRatioDn[iiter].GetXaxis().SetTitle("Top quark p_{T} (GeV)")
    hErrorRatioDn[iiter].SetLineColor(color[iiter])
    hErrorRatioDn[iiter].Draw("same")
cleg.Draw()
cerrDn.SaveAs("UnfoldingPlots/errorRatioDn_"+options.lepType+"_"+options.type+".pdf")

# -------------------------------------------------------------------------------------
# Translate to cross section (not events) in bins of pt N/L/BR)
# -------------------------------------------------------------------------------------
# TODO: should fix BR

if options.toy == "up" :
    thisReco = hRecoUp_iter[best_iter-1]
    thisTrue = hTrueUp
    thisMeas = hMeasUp
elif options.toy == "dn" :
    thisReco = hRecoDn_iter[best_iter-1]
    thisTrue = hTrueDn
    thisMeas = hMeasDn
else :
    thisReco = hReco_iter[best_iter-1]
    thisTrue = hTrue
    thisMeas = hMeas

thisTrue.Scale(1.0/(lum*0.438/3.)) # true @ parton level
thisMeas.Scale(1.0/(lum*0.438/3.)) # measured @ reco level
thisReco.Scale(1.0/(lum*0.438/3.)) # unfolded to parton level

# -------------------------------------------------------------------------------------
# Adjust for bin width
# -------------------------------------------------------------------------------------

for ibin in range(1, thisTrue.GetXaxis().GetNbins()+1 ) :

    width = thisTrue.GetBinWidth(ibin)

    thisTrue.SetBinContent(ibin, thisTrue.GetBinContent(ibin) / width )
    thisTrue.SetBinError(ibin, thisTrue.GetBinError(ibin) / width )

    thisMeas.SetBinContent(ibin,  thisMeas.GetBinContent(ibin) / width )
    thisMeas.SetBinError(ibin,  thisMeas.GetBinError(ibin) / width )

    thisReco.SetBinContent(ibin, thisReco.GetBinContent(ibin) / width )
    thisReco.SetBinError(ibin, thisReco.GetBinError(ibin) / width )

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
text1.DrawLatex(0.55,0.8, "#scale[1.0]{L = 35.9 fb^{-1}, #sqrt{s} = 13 TeV}")

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
hFrac.GetXaxis().SetRangeUser(400.,1200.)

c1.Update()

c1.SaveAs("UnfoldingPlots/closure_"+options.lepType+"_"+options.toy+"_"+options.type+"_result.pdf")

fout.Close()
