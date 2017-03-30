#! /usr/bin/env python
import os
import glob
import time
import math


# -------------------------------------------------------------------------------------
# define selection values
# -------------------------------------------------------------------------------------

# muons
MIN_MU_PT  = 50.0
MAX_MU_ETA = 2.1

# electrons
MIN_EL_PT  = 50.0
MAX_EL_ETA = 2.5

# jets
MIN_JET_PT  = 50.0
MAX_JET_ETA = 2.4

TOP_PT_CUT = 350.0

MIN_TOP_MASS = 105.0
MAX_TOP_MASS = 220.0
# -------------------------------------------------------------------------------------
# define input options
# -------------------------------------------------------------------------------------

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--files', metavar='F', type='string', action='store',
                  dest='files',
                  help='Input files')

parser.add_option('--outname', metavar='F', type='string', action='store',
                  default='test_iheartNY',
                  dest='outname',
                  help='Name of output file')

parser.add_option('--mttGenMax', metavar='J', type='float', action='store',
                  default=None,
                  dest='mttGenMax',
                  help='Maximum generator-level m_ttbar [GeV] to stitch together the ttbar samples.')

parser.add_option('--semilep', metavar='J', type='float',action='store',
                  default=None,
                  dest='semilep',
                  help='Select only semileptonic ttbar decays (1) or only non-semileptonic ttbar decays (-1) or no such cut (None)')

parser.add_option('--usePuppi', metavar='M', action='store_true',
                  default=False,
                  dest='usePuppi',
                  help='Use Puppi jets instead of CHS')

parser.add_option('--isMC', metavar='M', action='store_true',
                  default=False,
                  dest='isMC',
                  help='Running on Monte Carlo')

parser.add_option('--jecIOV', metavar='F', type='string', action='store',
                  default=None,
                  dest='jecIOV',
                  help='Run period for JEC/JER (BCD,EF,G,H)')

parser.add_option('--debug', metavar='M', action='store_true',
                  default=False,
                  dest='debug',
                  help='Print out debug statements')

parser.add_option('--fullTruth', metavar='M', action='store_true',
                  default=False,
                  dest='fullTruth',
                  help='Save additional tree with truth info for fail-reco events')

parser.add_option('--puFile', metavar='F', type='string', action='store',
                  default=None,
                  dest='puFile',
                  help='Pileup reweighting file')

parser.add_option('--muOrEl', metavar='F', type='string', action='store',
                  default=None,
                  dest='muOrEl',
                  help='Muon or electron sample')

(options, args) = parser.parse_args()
argv = []

import ROOT
ROOT.gROOT.Macro("rootlogon.C")

from array import *

import sys
from DataFormats.FWLite import Events, Handle, Runs

if options.fullTruth and not (options.isMC and options.semilep == 1):
    sys.exit("ERROR: can only use fullTruth option with semilep TTbar MC!\n")
    
# -------------------------------------------------------------------------------------
# jet energy corrections
# -------------------------------------------------------------------------------------

ROOT.gSystem.Load('libCondFormatsJetMETObjects')

jetname = "chs"
if options.usePuppi:
    jetname = "Puppi"

if options.isMC : 
    L3JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L3Absolute_AK4PF'+jetname+'.txt')
    L3JetParAK4 = ROOT.JetCorrectorParameters(L3JecStrAK4)
    L2JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L2Relative_AK4PF'+jetname+'.txt')
    L2JetParAK4 = ROOT.JetCorrectorParameters(L2JecStrAK4)
    L1JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L1FastJet_AK4PF'+jetname+'.txt')
    L1JetParAK4 = ROOT.JetCorrectorParameters(L1JecStrAK4)
    UncJecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_Uncertainty_AK4PF'+jetname+'.txt')
    UncertJetAK4 = ROOT.JetCorrectionUncertainty(UncJecStrAK4)
    L3JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L3Absolute_AK8PF'+jetname+'.txt')
    L3JetParAK8 = ROOT.JetCorrectorParameters(L3JecStrAK8)
    L2JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L2Relative_AK8PF'+jetname+'.txt')
    L2JetParAK8 = ROOT.JetCorrectorParameters(L2JecStrAK8)
    L1JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L1FastJet_AK8PF'+jetname+'.txt')
    L1JetParAK8 = ROOT.JetCorrectorParameters(L1JecStrAK8)
    UncJecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_Uncertainty_AK8PF'+jetname+'.txt')
    UncertJetAK8 = ROOT.JetCorrectionUncertainty(UncJecStrAK8)
    
else :
    L3JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L3Absolute_AK4PF'+jetname+'.txt')
    L3JetParAK4 = ROOT.JetCorrectorParameters(L3JecStrAK4)
    L2JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L2Relative_AK4PF'+jetname+'.txt')
    L2JetParAK4 = ROOT.JetCorrectorParameters(L2JecStrAK4)
    L1JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L1FastJet_AK4PF'+jetname+'.txt')
    L1JetParAK4 = ROOT.JetCorrectorParameters(L1JecStrAK4)
    ResJecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L2L3Residual_AK4PF'+jetname+'.txt')
    ResJetParAK4 = ROOT.JetCorrectorParameters(ResJecStrAK4)
    UncJecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_Uncertainty_AK4PF'+jetname+'.txt')
    UncertJetAK4 = ROOT.JetCorrectionUncertainty(UncJecStrAK4)
    L3JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L3Absolute_AK8PF'+jetname+'.txt')
    L3JetParAK8 = ROOT.JetCorrectorParameters(L3JecStrAK8)
    L2JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L2Relative_AK8PF'+jetname+'.txt')
    L2JetParAK8 = ROOT.JetCorrectorParameters(L2JecStrAK8)
    L1JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L1FastJet_AK8PF'+jetname+'.txt')
    L1JetParAK8 = ROOT.JetCorrectorParameters(L1JecStrAK8)
    ResJecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_L2L3Residual_AK8PF'+jetname+'.txt')
    ResJetParAK8 = ROOT.JetCorrectorParameters(ResJecStrAK8)
    UncJecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016'+options.jecIOV+'V3_DATA_Uncertainty_AK8PF'+jetname+'.txt')
    UncertJetAK8 = ROOT.JetCorrectionUncertainty(UncJecStrAK8)

#  load JetCorrectorParameter objects into vector (order matters!)
vParJecAK4 = ROOT.std.vector(ROOT.JetCorrectorParameters)()
vParJecAK4.push_back(L1JetParAK4)
vParJecAK4.push_back(L2JetParAK4)
vParJecAK4.push_back(L3JetParAK4)
if not options.isMC : 
    vParJecAK4.push_back(ResJetParAK4)

ak4JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK4)

vParJecAK8 = ROOT.std.vector(ROOT.JetCorrectorParameters)()
vParJecAK8.push_back(L1JetParAK8)
vParJecAK8.push_back(L2JetParAK8)
vParJecAK8.push_back(L3JetParAK8)
if not options.isMC : 
    vParJecAK8.push_back(ResJetParAK8)

ak8JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK8)

    
# -------------------------------------------------------------------------------------
# define helper classes that use ROOT
# -------------------------------------------------------------------------------------

def findClosestInList( p41, p4list ) :
    minDR = 9999.
    ret = None
    for j in range(0,len(p4list) ):
        dR = p4list[j].DeltaR(p41)
        if dR < minDR :
            minDR = dR
            ret = p4list[j]
    return ret


class GenTopQuark :
    pdgId = 6                    # 6 = top, -6 = antitop
    p4 = ROOT.TLorentzVector()
    decay = 0                    # 0 = hadronic, 1 = leptonic
    def __init__( self, pdgId, p4, decay ) :
        self.pdgId = pdgId
        self.p4 = p4
        self.decay = decay
    def match( self, jets ) :
        return findClosestInList( self.p4, jets )



# -------------------------------------------------------------------------------------
# input and output files
# -------------------------------------------------------------------------------------

filelist = file( options.files )
filesraw = filelist.readlines()
files = []

for ifile in filesraw : 
    if len( ifile ) > 2 :
        s = ifile.rstrip()
        files.append( s )

print 'Running on files: '
print files

# @@@ Create output root file

f = ROOT.TFile(options.outname+".root", "RECREATE")
f.cd()

#~ Tree initializations

recoTree = ROOT.TTree("recoTree", "recoTree")

# parton level
genTopPt        = ROOT.vector('float')()
genTopEta       = ROOT.vector('float')()
genTopPhi       = ROOT.vector('float')()
genMuPt         = ROOT.vector('float')()
genMuEta        = ROOT.vector('float')()
genMuPhi        = ROOT.vector('float')()
genElPt         = ROOT.vector('float')()
genElEta        = ROOT.vector('float')()
genElPhi        = ROOT.vector('float')()
genTTbarMass    = ROOT.vector('float')()

# particle level
#genAK4jetPt     = ROOT.vector('float')()
#genAK4jetEta    = ROOT.vector('float')()
#genAK4jetPhi    = ROOT.vector('float')()
#genAK4jetMass   = ROOT.vector('float')()
genAK8jetPt     = ROOT.vector('float')()
genAK8jetEta    = ROOT.vector('float')()
genAK8jetPhi    = ROOT.vector('float')()
genAK8jetMass   = ROOT.vector('float')()
partMuPt        = ROOT.vector('float')()
partMuEta       = ROOT.vector('float')()
partMuPhi       = ROOT.vector('float')()
partElPt        = ROOT.vector('float')()
partElEta       = ROOT.vector('float')()
partElPhi       = ROOT.vector('float')()

# reco level
metPt                  = ROOT.vector('float')()
metPhi                 = ROOT.vector('float')()
ht                     = ROOT.vector('float')()
muPt                   = ROOT.vector('float')()
muEta                  = ROOT.vector('float')()
muPhi                  = ROOT.vector('float')()
muMiniIso              = ROOT.vector('float')()
muPtRelPt15            = ROOT.vector('float')()
muPtRelPt20            = ROOT.vector('float')()
muPtRelPt25            = ROOT.vector('float')()
muPtRelPt30            = ROOT.vector('float')()
muPtRelPt35            = ROOT.vector('float')()
muPtRelPt40            = ROOT.vector('float')()
muPtRelPt45            = ROOT.vector('float')()
mudRPt15               = ROOT.vector('float')()
mudRPt20               = ROOT.vector('float')()
mudRPt25               = ROOT.vector('float')()
mudRPt30               = ROOT.vector('float')()
mudRPt35               = ROOT.vector('float')()
mudRPt40               = ROOT.vector('float')()
mudRPt45               = ROOT.vector('float')()
muTight                = ROOT.vector('int')()
muCharge               = ROOT.vector('float')()
elPt                   = ROOT.vector('float')()
elEta                  = ROOT.vector('float')()
elPhi                  = ROOT.vector('float')()
elMiniIso              = ROOT.vector('float')()
elPtRelPt15            = ROOT.vector('float')()
elPtRelPt20            = ROOT.vector('float')()
elPtRelPt25            = ROOT.vector('float')()
elPtRelPt30            = ROOT.vector('float')()
elPtRelPt35            = ROOT.vector('float')()
elPtRelPt40            = ROOT.vector('float')()
elPtRelPt45            = ROOT.vector('float')()
eldRPt15               = ROOT.vector('float')()
eldRPt20               = ROOT.vector('float')()
eldRPt25               = ROOT.vector('float')()
eldRPt30               = ROOT.vector('float')()
eldRPt35               = ROOT.vector('float')()
eldRPt40               = ROOT.vector('float')()
eldRPt45               = ROOT.vector('float')()
elTight                = ROOT.vector('int')()
elPassConVeto          = ROOT.vector('int')()
elCharge               = ROOT.vector('float')()
ak4jetPt               = ROOT.vector('float')()
ak4jetEta              = ROOT.vector('float')()
ak4jetPhi              = ROOT.vector('float')()
ak4jetMass             = ROOT.vector('float')()
ak4jetCSV              = ROOT.vector('float')()
ak4jetVtxMass          = ROOT.vector('float')()
ak4jetJECunc           = ROOT.vector('float')()
ak4jetPtJERup          = ROOT.vector('float')()
ak4jetPtJERdown        = ROOT.vector('float')()
ak4jetEtaJERup         = ROOT.vector('float')()
ak4jetEtaJERdown       = ROOT.vector('float')()
ak4jetPhiJERup         = ROOT.vector('float')()
ak4jetPhiJERdown       = ROOT.vector('float')()
ak4jetMassJERup        = ROOT.vector('float')()
ak4jetMassJERdown      = ROOT.vector('float')()
ak4jetHadronFlavour    = ROOT.vector('float')()
ak8jetPt               = ROOT.vector('float')()
ak8jetEta              = ROOT.vector('float')()
ak8jetPhi              = ROOT.vector('float')()
ak8jetMass             = ROOT.vector('float')()
ak8jetMassPruned       = ROOT.vector('float')()
ak8jetMassFiltered     = ROOT.vector('float')()
ak8jetMassTrimmed      = ROOT.vector('float')()   
ak8jetTau1             = ROOT.vector('float')()
ak8jetTau2             = ROOT.vector('float')()
ak8jetTau3             = ROOT.vector('float')()
ak8jetCSV              = ROOT.vector('float')()
ak8jetSDmass           = ROOT.vector('float')()
ak8jetSDsubjet0pt      = ROOT.vector('float')()
ak8jetSDsubjet0eta     = ROOT.vector('float')()
ak8jetSDsubjet0phi     = ROOT.vector('float')()
ak8jetSDsubjet0mass    = ROOT.vector('float')()
ak8jetSDsubjet0CSV     = ROOT.vector('float')()
ak8jetSDsubjet1pt      = ROOT.vector('float')()
ak8jetSDsubjet1eta     = ROOT.vector('float')()
ak8jetSDsubjet1phi     = ROOT.vector('float')()
ak8jetSDsubjet1mass    = ROOT.vector('float')()
ak8jetSDsubjet1CSV     = ROOT.vector('float')()
ak8jetJECunc           = ROOT.vector('float')()
ak8jetPtJERup          = ROOT.vector('float')()
ak8jetPtJERdown        = ROOT.vector('float')()
ak8jetEtaJERup         = ROOT.vector('float')()
ak8jetEtaJERdown       = ROOT.vector('float')()
ak8jetPhiJERup         = ROOT.vector('float')()
ak8jetPhiJERdown       = ROOT.vector('float')()
ak8jetMassJERup        = ROOT.vector('float')()
ak8jetMassJERdown      = ROOT.vector('float')()
eventWeight_nom        = ROOT.vector('float')()
eventWeight_puUp       = ROOT.vector('float')()
eventWeight_puDown     = ROOT.vector('float')()
eventWeight_Q2Up       = ROOT.vector('float')()
eventWeight_Q2Down     = ROOT.vector('float')()
eventWeight_PDF        = ROOT.vector('float')()
eventWeight_alphaUp    = ROOT.vector('float')()
eventWeight_alphaDown  = ROOT.vector('float')()
truthChannel       = ROOT.vector('int')()

if options.semilep == 1:
    recoTree.Branch('genTopPt'               , genTopPt               )
    recoTree.Branch('genTopEta'              , genTopEta              )
    recoTree.Branch('genTopPhi'              , genTopPhi              )
    recoTree.Branch('genMuPt'                , genMuPt                )
    recoTree.Branch('genMuEta'               , genMuEta               )
    recoTree.Branch('genMuPhi'               , genMuPhi               )
    recoTree.Branch('genElPt'                , genElPt                )
    recoTree.Branch('genElEta'               , genElEta               )
    recoTree.Branch('genElPhi'               , genElPhi               )
    recoTree.Branch('genTTbarMass'           , genTTbarMass           )
    
    #recoTree.Branch('genAK4jetPt'            , genAK4jetPt            )
    #recoTree.Branch('genAK4jetEta'           , genAK4jetEta           )
    #recoTree.Branch('genAK4jetPhi'           , genAK4jetPhi           )
    #recoTree.Branch('genAK4jetMass'          , genAK4jetMass          )
    recoTree.Branch('genAK8jetPt'            , genAK8jetPt            )
    recoTree.Branch('genAK8jetEta'           , genAK8jetEta           )
    recoTree.Branch('genAK8jetPhi'           , genAK8jetPhi           )
    recoTree.Branch('genAK8jetMass'          , genAK8jetMass          )
    recoTree.Branch('partMuPt'               , partMuPt               )
    recoTree.Branch('partMuEta'              , partMuEta              )
    recoTree.Branch('partMuPhi'              , partMuPhi              )
    recoTree.Branch('partElPt'               , partElPt               )
    recoTree.Branch('partElEta'              , partElEta              )
    recoTree.Branch('partElPhi'              , partElPhi              )

recoTree.Branch('metPt'                  , metPt                  )
recoTree.Branch('metPhi'                 , metPhi                 )
recoTree.Branch('ht'                     , ht                     )
recoTree.Branch('muPt'                   , muPt                   )
recoTree.Branch('muEta'                  , muEta                  )
recoTree.Branch('muPhi'                  , muPhi                  )
recoTree.Branch('muMiniIso'              , muMiniIso              )
recoTree.Branch('muPtRelPt15'            , muPtRelPt15            )
recoTree.Branch('muPtRelPt20'            , muPtRelPt20            )
recoTree.Branch('muPtRelPt25'            , muPtRelPt25            )
recoTree.Branch('muPtRelPt30'            , muPtRelPt30            )
recoTree.Branch('muPtRelPt35'            , muPtRelPt35            )
recoTree.Branch('muPtRelPt40'            , muPtRelPt40            )
recoTree.Branch('muPtRelPt45'            , muPtRelPt45            )
recoTree.Branch('mudRPt15'               , mudRPt15               )
recoTree.Branch('mudRPt20'               , mudRPt20               )
recoTree.Branch('mudRPt25'               , mudRPt25               )
recoTree.Branch('mudRPt30'               , mudRPt30               )
recoTree.Branch('mudRPt35'               , mudRPt35               )
recoTree.Branch('mudRPt40'               , mudRPt40               )
recoTree.Branch('mudRPt45'               , mudRPt45               )
recoTree.Branch('muTight'                , muTight                )
recoTree.Branch('muCharge'               , muCharge               )
recoTree.Branch('elPt'                   , elPt                   )
recoTree.Branch('elEta'                  , elEta                  )
recoTree.Branch('elPhi'                  , elPhi                  )
recoTree.Branch('elMiniIso'              , elMiniIso              )
recoTree.Branch('elPtRelPt15'            , elPtRelPt15            )
recoTree.Branch('elPtRelPt20'            , elPtRelPt20            )
recoTree.Branch('elPtRelPt25'            , elPtRelPt25            )
recoTree.Branch('elPtRelPt30'            , elPtRelPt30            )
recoTree.Branch('elPtRelPt35'            , elPtRelPt35            )
recoTree.Branch('elPtRelPt40'            , elPtRelPt40            )
recoTree.Branch('elPtRelPt45'            , elPtRelPt45            )
recoTree.Branch('eldRPt15'               , eldRPt15               )
recoTree.Branch('eldRPt20'               , eldRPt20               )
recoTree.Branch('eldRPt25'               , eldRPt25               )
recoTree.Branch('eldRPt30'               , eldRPt30               )
recoTree.Branch('eldRPt35'               , eldRPt35               )
recoTree.Branch('eldRPt40'               , eldRPt40               )
recoTree.Branch('eldRPt45'               , eldRPt45               )
recoTree.Branch('elTight'                , elTight                )
recoTree.Branch('elPassConVeto'          , elPassConVeto          )
recoTree.Branch('elCharge'               , elCharge               )
recoTree.Branch('ak4jetPt'               , ak4jetPt               )
recoTree.Branch('ak4jetEta'              , ak4jetEta              )
recoTree.Branch('ak4jetPhi'              , ak4jetPhi              )
recoTree.Branch('ak4jetMass'             , ak4jetMass             )
recoTree.Branch('ak4jetCSV'              , ak4jetCSV              )
recoTree.Branch('ak4jetVtxMass'          , ak4jetVtxMass          )
recoTree.Branch('ak4jetJECunc'           , ak4jetJECunc           )
if options.isMC:
    recoTree.Branch('ak4jetHadronFlavour'    , ak4jetHadronFlavour    )
    recoTree.Branch('ak4jetPtJERup'          , ak4jetPtJERup          )
    recoTree.Branch('ak4jetPtJERdown'        , ak4jetPtJERdown        )
    recoTree.Branch('ak4jetEtaJERup'         , ak4jetEtaJERup         )
    recoTree.Branch('ak4jetEtaJERdown'       , ak4jetEtaJERdown       )
    recoTree.Branch('ak4jetPhiJERup'         , ak4jetPhiJERup         )
    recoTree.Branch('ak4jetPhiJERdown'       , ak4jetPhiJERdown       )
    recoTree.Branch('ak4jetMassJERup'        , ak4jetMassJERup        )
    recoTree.Branch('ak4jetMassJERdown'      , ak4jetMassJERdown      )
recoTree.Branch('ak8jetPt'               , ak8jetPt               )
recoTree.Branch('ak8jetEta'              , ak8jetEta              )
recoTree.Branch('ak8jetPhi'              , ak8jetPhi              )
recoTree.Branch('ak8jetMass'             , ak8jetMass          )
recoTree.Branch('ak8jetMassPruned'       , ak8jetMassPruned       )
recoTree.Branch('ak8jetMassFiltered'     , ak8jetMassFiltered     )
recoTree.Branch('ak8jetMassTrimmed'      , ak8jetMassTrimmed      )
recoTree.Branch('ak8jetTau1'             , ak8jetTau1             )
recoTree.Branch('ak8jetTau2'             , ak8jetTau2             )
recoTree.Branch('ak8jetTau3'             , ak8jetTau3             )
recoTree.Branch('ak8jetCSV'              , ak8jetCSV              )
recoTree.Branch('ak8jetSDmass'           , ak8jetSDmass           )
recoTree.Branch('ak8jetSDsubjet0pt'      , ak8jetSDsubjet0pt      )
recoTree.Branch('ak8jetSDsubjet0eta'     , ak8jetSDsubjet0eta     )
recoTree.Branch('ak8jetSDsubjet0phi'     , ak8jetSDsubjet0phi     )
recoTree.Branch('ak8jetSDsubjet0mass'    , ak8jetSDsubjet0mass    )
recoTree.Branch('ak8jetSDsubjet0CSV'     , ak8jetSDsubjet0CSV     )
recoTree.Branch('ak8jetSDsubjet1pt'      , ak8jetSDsubjet1pt      ) 
recoTree.Branch('ak8jetSDsubjet1eta'     , ak8jetSDsubjet1eta     )
recoTree.Branch('ak8jetSDsubjet1phi'     , ak8jetSDsubjet1phi     )
recoTree.Branch('ak8jetSDsubjet1mass'    , ak8jetSDsubjet1mass    )
recoTree.Branch('ak8jetSDsubjet1CSV'     , ak8jetSDsubjet1CSV     )
recoTree.Branch('ak8jetJECunc'           , ak8jetJECunc           )
if options.isMC:
    recoTree.Branch('ak8jetPtJERup'          , ak8jetPtJERup          )
    recoTree.Branch('ak8jetPtJERdown'        , ak8jetPtJERdown        )
    recoTree.Branch('ak8jetEtaJERup'         , ak8jetEtaJERup         )
    recoTree.Branch('ak8jetEtaJERdown'       , ak8jetEtaJERdown       )
    recoTree.Branch('ak8jetPhiJERup'         , ak8jetPhiJERup         )
    recoTree.Branch('ak8jetPhiJERdown'       , ak8jetPhiJERdown       )
    recoTree.Branch('ak8jetMassJERup'        , ak8jetMassJERup        )
    recoTree.Branch('ak8jetMassJERdown'      , ak8jetMassJERdown      )

recoTree.Branch('eventWeight_nom'           , eventWeight_nom           )

if options.isMC:
    recoTree.Branch('eventWeight_puUp'          , eventWeight_puUp          )
    recoTree.Branch('eventWeight_puDown'        , eventWeight_puDown        )
    if options.semilep is not None:
        recoTree.Branch('eventWeight_Q2Up'   , eventWeight_Q2Up )
        recoTree.Branch('eventWeight_Q2Down' , eventWeight_Q2Down )
        recoTree.Branch('eventWeight_PDF'    , eventWeight_PDF )
        recoTree.Branch('eventWeight_alphaUp', eventWeight_alphaUp)
        recoTree.Branch('eventWeight_alphaDown', eventWeight_alphaDown)
    if options.semilep == 1:
        recoTree.Branch('truthChannel'       , truthChannel           )

if options.fullTruth:
    trueTree = ROOT.TTree("trueTree", "trueTree")
    trueTree.Branch('genTopPt'               , genTopPt               )
    trueTree.Branch('genTopEta'              , genTopEta              )
    trueTree.Branch('genTopPhi'              , genTopPhi              )
    trueTree.Branch('genTTbarMass'           , genTTbarMass           )
    trueTree.Branch('truthChannel'           , truthChannel           )
    trueTree.Branch('eventWeight_nom'        , eventWeight_nom        )
    trueTree.Branch('eventWeight_puUp'       , eventWeight_puUp       )
    trueTree.Branch('eventWeight_puDown'     , eventWeight_puDown     )
    trueTree.Branch('eventWeight_Q2Up'       , eventWeight_Q2Up       )
    trueTree.Branch('eventWeight_Q2Down'     , eventWeight_Q2Down     )
    trueTree.Branch('eventWeight_PDF'        , eventWeight_PDF        )
    trueTree.Branch('eventWeight_alphaUp'    , eventWeight_alphaUp    )
    trueTree.Branch('eventWeight_alphaDown'  , eventWeight_alphaDown  )

# -------------------------------------------------------------------------------------
# define all variables to be read from input files
# -------------------------------------------------------------------------------------

events = Events (files)
runs = Runs ( files)

jetname = "CHS"
massType = "CHS"
if options.usePuppi:
    jetname = "Puppi"
    massType = ""

pvChiHandle  = Handle("std::vector<float>")
pvChiLabel   = ( "vertexInfo", "chi" )
pvRhoHandle  = Handle("std::vector<float>")
pvRhoLabel   = ( "vertexInfo", "rho" )
pvZHandle    = Handle("std::vector<float>")
pvZLabel     = ( "vertexInfo", "z" )
pvNdofHandle = Handle("std::vector<std::int>")
pvNdofLabel  = ( "vertexInfo", "ndof" )

trigNameHandle = Handle( "std::vector<std::string>")
trigNameLabel = ("TriggerUserData", "triggerNameTree")
trigBitsHandle = Handle( "std::vector<float>")
trigBitsLabel = ("TriggerUserData", "triggerBitTree")
trigPrescalesHandle = Handle( "std::vector<int>")
trigPrescalesLabel = ("TriggerUserData", "triggerPrescaleTree")

metFilterNameHandle = Handle( "std::vector<std::string>")
metFilterNameLabel  = ("METUserData", "triggerNameTree")
metFilterBitsHandle = Handle( "std::vector<float>")
metFilterBitsLabel  = ("METUserData", "triggerBitTree")
BadChargedCandidateFilterHandle    = Handle("bool")
BadChargedCandidateFilterLabel     = ("BadChargedCandidateFilter", "")
BadPFMuonFilterHandle = Handle("bool")
BadPFMuonFilterLabel = ("BadPFMuonFilter","")

LHERunHandle = Handle("LHERunInfoProduct")
LHERunLabel  = ("externalLHEProducer","")
LHEHandle = Handle("LHEEventProduct")
LHELabel = ("externalLHEProducer","")

# genParticles
if options.isMC and (options.semilep is not None):
    genParticlesHandle     = Handle("std::vector<reco::GenParticle>")
    genParticlesLabel      = ("filteredPrunedGenParticles", "")
    
    ## only have AK8 gen jets currently (in addition to gen jets matched to AK4/AK8 reco jets)
    ak8GenJetPtHandle   = Handle("std::vector<float>")
    ak8GenJetPtLabel    = ("genJetsAK8", "genJetsAK8Pt")
    ak8GenJetEtaHandle  = Handle("std::vector<float>")
    ak8GenJetEtaLabel   = ("genJetsAK8", "genJetsAK8Eta")
    ak8GenJetPhiHandle  = Handle("std::vector<float>")
    ak8GenJetPhiLabel   = ("genJetsAK8", "genJetsAK8Phi")
    ak8GenJetEHandle    = Handle("std::vector<float>")
    ak8GenJetELabel     = ("genJetsAK8", "genJetsAK8E")
    
metHandle    = Handle("std::vector<float>")
metLabel     = ("metFull", "metFullPt")
metphiHandle = Handle("std::vector<float>")
metphiLabel  = ("metFull", "metFullPhi")

# lepton variables
muonPtHandle      = Handle("std::vector<float>")
muonPtLabel       = ("muons", "muPt")
muonEtaHandle     = Handle( "std::vector<float>")
muonEtaLabel      = ("muons", "muEta")
muonPhiHandle     = Handle( "std::vector<float>")
muonPhiLabel      = ("muons", "muPhi")
muonChargeHandle  = Handle( "std::vector<float>")
muonChargeLabel   = ("muons", "muCharge")
muTightHandle     = Handle("std::vector<float>")
muTightLabel      = ("muons" , "muIsTightMuon" )
muMediumHandle    = Handle("std::vector<float>")
muMediumLabel     = ("muons" , "muIsMediumMuon" )
muLooseHandle     = Handle("std::vector<float>")
muLooseLabel      = ("muons" , "muIsLooseMuon" )
muPFHandle        = Handle( "std::vector<float>" )
muPFLabel         = ("muons" , "muIsPFMuon" )
muGlobalHandle    = Handle( "std::vector<float>" )
muGlobalLabel     = ("muons" , "muIsGlobalMuon" )
muTrackerHandle   = Handle( "std::vector<float>" )
muTrackerLabel    = ("muons" , "muIsTrackerMuon" )
muTrkChi2Handle   = Handle( "std::vector<float>" )
muTrkChi2Label    = ("muons" , "muGlbTrkNormChi2" )
muHitHandle       = Handle( "std::vector<float>" )
muHitLabel        = ("muons" , "muNumberValidMuonHits")
muMatchStatHandle = Handle( "std::vector<float>" )
muMatchStatLabel  = ("muons" , "muNumberMatchedStations" )
muPixHitHandle    = Handle( "std::vector<float>" )
muPixHitLabel     = ("muons" , "muNumberValidPixelHits" )
muTrkLayerHandle  = Handle( "std::vector<float>" ) 
muTrkLayerLabel   = ("muons" , "muNumberTrackerLayers")
muDxyHandle       = Handle( "std::vector<float>" )
muDxyLabel        = ("muons" , "muDxy")
muDzHandle        = Handle( "std::vector<float>" )
muDzLabel         = ("muons" , "muDz" )
muMiniIsoHandle   = Handle( "std::vector<float>" )
muMiniIsoLabel    = ("muons" , "muMiniIso" )

muKeyHandle = Handle("std::vector<std::vector<int> >")
muKeyLabel = ("muonKeys")

electronPtHandle          = Handle( "std::vector<float>")
electronPtLabel           = ("electrons", "elPt")
electronEtaHandle         = Handle( "std::vector<float>")
electronEtaLabel          = ("electrons", "elEta")
electronPhiHandle         = Handle( "std::vector<float>")
electronPhiLabel          = ("electrons", "elPhi")
electronChargeHandle      = Handle( "std::vector<float>")
electronChargeLabel       = ("electrons", "elCharge")
electronTightHandle       = Handle("std::vector<float>")
electronTightLabel        = ("electrons" , "elvidTight" )
electronMediumHandle      = Handle("std::vector<float>")
electronMediumLabel       = ("electrons" , "elvidMedium" )
electronLooseHandle       = Handle("std::vector<float>")
electronLooseLabel        = ("electrons" , "elvidLoose" )
electronfull5x5sieeHandle = Handle("std::vector<float>")
electronfull5x5sieeLabel  = ( "electrons" , "elfull5x5siee")
electrondEtaInHandle      = Handle("std::vector<float>")
electrondEtaInLabel       = ( "electrons" , "eldEtaInSeed" )
electrondPhiInHandle      = Handle("std::vector<float>")
electrondPhiInLabel       = ( "electrons" , "eldPhiIn" )
electronHoEHandle         = Handle("std::vector<float>")
electronHoELabel          = ( "electrons" , "elHoE" )
electronooEmooPHandle     = Handle("std::vector<float>")
electronooEmooPLabel      = ( "electrons" , "elooEmooP")
electronMiniIsoHandle     = Handle("std::vector<float>")
electronMiniIsoLabel      = ( "electrons" , "elMiniIso" )
electronSCEtaHandle       = Handle("std::vector<float>")
electronSCEtaLabel        = ( "electrons" , "elSCEta")
electronMissHitsHandle    = Handle("std::vector<float>")
electronMissHitsLabel     = ( "electrons" , "elmissHits")
electronConVetoHandle     = Handle("std::vector<float>")
electronConVetoLabel      = ( "electrons" , "elhasMatchedConVeto")
electronMediumNoIsoHandle = Handle("std::vector<float>")
electronMediumNoIsoLabel  = ( "electrons" , "elvidMediumnoiso")
electronTightNoIsoHandle  = Handle("std::vector<float>")
electronTightNoIsoLabel   = ( "electrons" , "elvidTightnoiso")
electronDxyHandle         = Handle("std::vector<float>")
electronDxyLabel          = ( "electrons" , "elDxy")
electronDzHandle          = Handle("std::vector<float>")
electronDzLabel           = ( "electrons" , "elDz")

elKeyHandle = Handle("std::vector<std::vector<int> >")
elKeyLabel = ( "electronKeys" )

# AK4 jet collection
ak4JetPtHandle   = Handle( "std::vector<float>" )
ak4JetPtLabel    = ("jetsAK4"+jetname, "jetAK4"+jetname+"Pt")
ak4JetEtaHandle  = Handle( "std::vector<float>" )
ak4JetEtaLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"Eta")
ak4JetPhiHandle  = Handle( "std::vector<float>" )
ak4JetPhiLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"Phi")
ak4JetEHandle = Handle( "std::vector<float>" )
ak4JetELabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"E")
ak4JetCSVHandle  = Handle( "std::vector<float>" )
ak4JetCSVLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"CSVv2")
ak4JetVtxMassHandle = Handle( "std::vector<float>" )
ak4JetVtxMassLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"SV0mass")    
ak4JetAreaHandle = Handle( "std::vector<float>" )
ak4JetAreaLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"jetArea")

ak4JetNeuHadEnergyFracHandle = Handle("std::vector<float>")
ak4JetNeuHadEnergyFracLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"neutralHadronEnergyFrac")
ak4JetNeuEmEnergyFracHandle = Handle("std::vector<float>")
ak4JetNeuEmEnergyFracLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"neutralEmEnergyFrac")
ak4JetChHadEnergyFracHandle = Handle("std::vector<float>")
ak4JetChHadEnergyFracLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"chargedHadronEnergyFrac")
ak4JetChEmEnergyFracHandle = Handle("std::vector<float>")
ak4JetChEmEnergyFracLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"chargedEmEnergyFrac")
ak4JetNeuMultiHandle = Handle("std::vector<float>")
ak4JetNeuMultiLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"neutralMultiplicity")
ak4JetChMultiHandle = Handle("std::vector<float>")
ak4JetChMultiLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"chargedMultiplicity")

ak4JetJECHandle = Handle("std::vector<float>")
ak4JetJECLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"jecFactor0")

# JER variables
ak4MatchedGenJetPtHandle   = Handle("std::vector<float>")
ak4MatchedGenJetPtLabel    = ("jetsAK4"+jetname, "jetAK4"+jetname+"GenJetPt")
ak4MatchedGenJetEtaHandle  = Handle("std::vector<float>")
ak4MatchedGenJetEtaLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"GenJetEta")
ak4MatchedGenJetPhiHandle  = Handle("std::vector<float>")
ak4MatchedGenJetPhiLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"GenJetPhi")
ak4MatchedGenJetEnergyHandle = Handle("std::vector<float>")
ak4MatchedGenJetEnergyLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"GenJetE")
ak4JERSFnomHandle  = Handle("std::vector<float>")
ak4JERSFnomLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"JERSF")
ak4JERSFupHandle   = Handle("std::vector<float>")
ak4JERSFupLabel    = ("jetsAK4"+jetname, "jetAK4"+jetname+"JERSFUp")
ak4JERSFdownHandle = Handle("std::vector<float>")
ak4JERSFdownLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"JERSFDown")
ak4JERHandle       = Handle("std::vector<float>")
ak4JERLabel        = ("jetsAK4"+jetname, "jetAK4"+jetname+"PtResolution")
ak4JetHadronFlavourHandle = Handle("std::vector<float>")
ak4JetHadronFlavourLabel = ("jetsAK4"+jetname, "jetAK4"+jetname+"HadronFlavour")

ak4JetKeysHandle = Handle("std::vector<std::vector<int> >")
ak4JetKeysLabel = ( "jetKeysAK4"+jetname , "" )

# AK8 jet collection
ak8JetPtHandle   = Handle( "std::vector<float>" )
ak8JetPtLabel    = ("jetsAK8"+jetname, "jetAK8"+jetname+"Pt")
ak8JetEtaHandle  = Handle( "std::vector<float>" )
ak8JetEtaLabel   = ("jetsAK8"+jetname, "jetAK8"+jetname+"Eta")
ak8JetPhiHandle  = Handle( "std::vector<float>" )
ak8JetPhiLabel   = ("jetsAK8"+jetname, "jetAK8"+jetname+"Phi")
ak8JetEHandle = Handle( "std::vector<float>" )
ak8JetELabel  = ("jetsAK8"+jetname, "jetAK8"+jetname+"E")
ak8JetTrimMassHandle = Handle("std::vector<float>")
ak8JetTrimMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"trimmedMass"+massType )
ak8JetPrunMassHandle = Handle("std::vector<float>")
ak8JetPrunMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"prunedMass"+massType )
ak8JetFiltMassHandle = Handle("std::vector<float>")
ak8JetFiltMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"filteredMass"+massType )
ak8JetTau1Handle = Handle("std::vector<float>")
ak8JetTau1Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau1"+massType )
ak8JetTau2Handle = Handle("std::vector<float>")
ak8JetTau2Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau2"+massType )
ak8JetTau3Handle = Handle("std::vector<float>")
ak8JetTau3Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau3"+massType )
ak8JetCSVHandle = Handle("std::vector<float>")               
ak8JetCSVLabel = ( "jetsAK8"+jetname , "jetAK8"+jetname+"CSVv2" )

ak8JetSoftDropMassHandle = Handle("std::vector<float>")
ak8JetSoftDropMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"softDropMass"+massType )
ak8JetSoftDropSubjet0Handle    = Handle("std::vector<float>")
ak8JetSoftDropSubjet0Label     = ("jetsAK8"+jetname, "jetAK8"+jetname+"vSubjetIndex0")
ak8JetSoftDropSubjet1Handle    = Handle("std::vector<float>")
ak8JetSoftDropSubjet1Label     = ("jetsAK8"+jetname, "jetAK8"+jetname+"vSubjetIndex1")

sjSoftDropPtHandle              = Handle( "std::vector<float>")
sjSoftDropPtLabel               = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"Pt")
sjSoftDropEtaHandle             = Handle( "std::vector<float>")
sjSoftDropEtaLabel              = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"Eta")
sjSoftDropPhiHandle             = Handle( "std::vector<float>")
sjSoftDropPhiLabel              = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"Phi")
sjSoftDropEHandle            = Handle( "std::vector<float>")
sjSoftDropELabel             = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"E")
sjSoftDropCSVHandle           = Handle( "std::vector<float>")
sjSoftDropCSVLabel            = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"CSVv2")

# JETID variables
ak8JetNeuHadEnergyFracHandle = Handle("std::vector<float>")
ak8JetNeuHadEnergyFracLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"neutralHadronEnergyFrac")
ak8JetNeuEmEnergyFracHandle = Handle("std::vector<float>")
ak8JetNeuEmEnergyFracLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"neutralEmEnergyFrac")
ak8JetChHadEnergyFracHandle = Handle("std::vector<float>")
ak8JetChHadEnergyFracLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"chargedHadronEnergyFrac")
ak8JetChEmEnergyFracHandle = Handle("std::vector<float>")
ak8JetChEmEnergyFracLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"chargedEmEnergyFrac")
ak8JetNeuMultiHandle = Handle("std::vector<float>")
ak8JetNeuMultiLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"neutralMultiplicity")
ak8JetChMultiHandle = Handle("std::vector<float>")
ak8JetChMultiLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"chargedMultiplicity")

#JEC
ak8JetJECHandle = Handle("std::vector<float>")
ak8JetJECLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"jecFactor0")
ak8JetAreaHandle = Handle( "std::vector<float>" )
ak8JetAreaLabel  = ("jetsAK8"+jetname, "jetAK8"+jetname+"jetArea")

#JER
ak8MatchedGenJetPtHandle   = Handle("std::vector<float>")
ak8MatchedGenJetPtLabel    = ("jetsAK8"+jetname, "jetAK8"+jetname+"GenJetPt")
ak8MatchedGenJetEtaHandle  = Handle("std::vector<float>")
ak8MatchedGenJetEtaLabel   = ("jetsAK8"+jetname, "jetAK8"+jetname+"GenJetEta")
ak8MatchedGenJetPhiHandle  = Handle("std::vector<float>")
ak8MatchedGenJetPhiLabel   = ("jetsAK8"+jetname, "jetAK8"+jetname+"GenJetPhi")
ak8MatchedGenJetEnergyHandle = Handle("std::vector<float>")
ak8MatchedGenJetEnergyLabel  = ("jetsAK8"+jetname, "jetAK8"+jetname+"GenJetE")
ak8JERSFnomHandle  = Handle("std::vector<float>")
ak8JERSFnomLabel   = ("jetsAK8"+jetname, "jetAK8"+jetname+"JERSF")
ak8JERSFupHandle   = Handle("std::vector<float>")
ak8JERSFupLabel    = ("jetsAK8"+jetname, "jetAK8"+jetname+"JERSFUp")
ak8JERSFdownHandle = Handle("std::vector<float>")
ak8JERSFdownLabel  = ("jetsAK8"+jetname, "jetAK8"+jetname+"JERSFDown")
ak8JERHandle       = Handle("std::vector<float>")
ak8JERLabel        = ("jetsAK8"+jetname, "jetAK8"+jetname+"PtResolution")

ak8JetKeysHandle = Handle("std::vector<std::vector<int> >")
ak8JetKeysLabel = ( "jetKeysAK8"+jetname , "" )

rhoHandle = Handle("double")
rhoLabel = ("fixedGridRhoFastjetAll", "")

puNtrueIntHandle = Handle("std::int")
puNtrueIntLabel = ( "eventUserData" , "puNtrueInt" )

runnumHandle = Handle("uint")
runnumLabel = ("eventInfo", "evtInfoRunNumber")

# -------------------------------------------------------------------------------------
# Histograms to save
# -------------------------------------------------------------------------------------

h_NtrueIntPU     = ROOT.TH1D("h_NtrueIntPU"    , "", 50,0,50 )
h_NPV_noweight   = ROOT.TH1F("h_NPV_noweight"  , "", 50,0,50 )
h_NPV            = ROOT.TH1F("h_NPV"           , "", 50,0,50 )
h_muPrescale     = ROOT.TH1F("h_muPrescale"    , "", 50,0,50 )
h_elPrescale     = ROOT.TH1F("h_elPrescale"    , "", 50,0,50 )
h_cutflow        = ROOT.TH1F("h_cutflow"       , "", 8,0.5,8.5 )

# -------------------------------------------------------------------------------------
# Get pileup weights
# -------------------------------------------------------------------------------------

if options.isMC and (options.puFile is not None):
    fPUweight      = ROOT.TFile(options.puFile)
    hPUweight_nom  = fPUweight.Get("PUweight_true")
    hPUweight_up   = fPUweight.Get("PUweight_up")
    hPUweight_down = fPUweight.Get("PUweight_down")

# -------------------------------------------------------------------------------------
# reset various counters
# -------------------------------------------------------------------------------------

ntotal = 0       # total number of events

nInIOV = 0

# Parton level checks
nMuChParton = 0
nElChParton = 0
nDilepParton = 0
nHadParton = 0
nTauOnlyParton = 0
nWrongHemisphere = 0

# Event quality
nPassNPV = 0
nPassLHE = 0
nPassMetFilter = 0
nPassRho = 0

nPassParton = 0

nPassParticle = 0

nPassLep = 0
nPassAK4jet = 0
nPassAK8jet = 0
nEventsPass = 0

nLooseNotManualEl = 0.0
nManualNotLooseEl = 0.0
nMediumNotManualEl = 0.0
nManualNotMediumEl = 0.0
nTightNotManualEl = 0.0
nManualNotTightEl = 0.0
nMediumNoIsoNotManualEl = 0.0
nManualNotMediumNoIsoEl = 0.0
nTightNoIsoNotManualEl = 0.0
nManualNotTightNoIsoEl = 0.0
nEleRaw = 0.0
nLooseNotManualMu = 0.0
nManualNotLooseMu = 0.0
nTightNotManualMu = 0.0
nManualNotTightMu = 0.0
nMuRaw = 0.0

nMatchAK4Jet = 0
nGoodMatchAK4Jet = 0
nMatchAK8Jet = 0
nGoodMatchAK8Jet = 0

smearfunc = ROOT.TRandom3()

# -------------------------------                                                                                                               
# Get MET filter names, indices of mu and el triggers   
# -------------------------------
muTrigIndices = {}
elTrigIndices = {}
metFiltIndices = {}

for run in runs:

    #Also get LHE information
    #if options.isMC :
    #    run.getByLabel( LHERunLabel, LHERunHandle)
    #    LHERunInfo = LHERunHandle.product()
    #    iter = LHERunInfo.headers_begin()
    #    while iter!=LHERunInfo.headers_end() :
            #print iter.tag()
    #        for line in iter.lines():
    #            if ('<weight' in line) or ('</weightgroup' in line):
    #                print line.rstrip()
    #        iter.next()

    runnumber = run.runAuxiliary().run()

    gotName = run.getByLabel( metFilterNameLabel, metFilterNameHandle )
    if not gotName:
        print 'Error! no MET filter in run info'
    filterNameStrings = metFilterNameHandle.product()
    metFilts = []
    if len(filterNameStrings) != 0:
        for ifilt in xrange(0, len(filterNameStrings) ) :
            if any(s in filterNameStrings[ifilt] for s in ("HBHENoiseFilter","HBHENoiseIsoFilter","globalTightHalo2016Filter","goodVertices","eeBadScFilter","EcalDeadCellTriggerPrimitiveFilter")):
                metFilts.append(ifilt)


    metFiltIndices[runnumber] = metFilts

    gotTrigName = run.getByLabel( trigNameLabel, trigNameHandle )
    if not gotTrigName:
        print 'Error! no trigger names in run info'
    triggerNameStrings = trigNameHandle.product()
    itrig = 0
    for name in triggerNameStrings:
        if "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50" in name :
            elTrigIndices[runnumber] = itrig
        if "HLT_Mu40_eta2p1_PFJet200_PFJet50" in name :
            muTrigIndices[runnumber] = itrig
        if options.debug :
            if "HLT_Ele" in name or "HLT_Mu" in name:
                print name
        itrig += 1

if options.debug :
    print muTrigIndices
    print elTrigIndices
    print metFiltIndices

# -------------------------------------------------------------------------------------
# start looping over events
# -------------------------------------------------------------------------------------

print "Start looping over events!"

for event in events :

    if options.semilep == 1:
        genTopPt.clear()
        genTopEta.clear()
        genTopPhi.clear()
        genMuPt.clear()
        genMuEta.clear()
        genMuPhi.clear()
        genElPt.clear()
        genElEta.clear()
        genElPhi.clear()
        genTTbarMass.clear()
        genAK8jetPt.clear()
        genAK8jetEta.clear()
        genAK8jetPhi.clear()
        genAK8jetMass.clear()
        partMuPt.clear()
        partMuEta.clear()
        partMuPhi.clear()
        partElPt.clear()
        partElEta.clear()
        partElPhi.clear()
    metPt.clear()
    metPhi.clear()
    ht.clear()
    muPt.clear()
    muEta.clear()
    muPhi.clear()
    muMiniIso.clear()
    muPtRelPt15.clear()
    muPtRelPt20.clear()
    muPtRelPt25.clear()
    muPtRelPt30.clear()
    muPtRelPt35.clear()
    muPtRelPt40.clear()
    muPtRelPt45.clear()
    mudRPt15.clear()
    mudRPt20.clear()
    mudRPt25.clear()
    mudRPt30.clear()
    mudRPt35.clear()
    mudRPt40.clear()
    mudRPt45.clear()
    muTight.clear()
    muCharge.clear()
    elPt.clear()
    elEta.clear()
    elPhi.clear()
    elMiniIso.clear()
    elPtRelPt15.clear()
    elPtRelPt20.clear()
    elPtRelPt25.clear()
    elPtRelPt30.clear()
    elPtRelPt35.clear()
    elPtRelPt40.clear()
    elPtRelPt45.clear()
    eldRPt15.clear()
    eldRPt20.clear()
    eldRPt25.clear()
    eldRPt30.clear()
    eldRPt35.clear()
    eldRPt40.clear()
    eldRPt45.clear()
    elTight.clear()
    elPassConVeto.clear()
    elCharge.clear()
    ak4jetPt.clear()
    ak4jetEta.clear()
    ak4jetPhi.clear()
    ak4jetMass.clear()
    ak4jetCSV.clear()
    ak4jetVtxMass.clear()
    ak4jetJECunc.clear()
    if options.isMC:
        ak4jetHadronFlavour.clear()
        ak4jetPtJERup.clear()
        ak4jetPtJERdown.clear()
        ak4jetEtaJERup.clear()
        ak4jetEtaJERdown.clear()
        ak4jetPhiJERup.clear()
        ak4jetPhiJERdown.clear()
        ak4jetMassJERup.clear()
        ak4jetMassJERdown.clear()
    ak8jetPt.clear()
    ak8jetEta.clear()
    ak8jetPhi.clear()
    ak8jetMass.clear()
    ak8jetMassPruned.clear()
    ak8jetMassFiltered.clear()
    ak8jetMassTrimmed.clear()   
    ak8jetTau1.clear()
    ak8jetTau2.clear()
    ak8jetTau3.clear()
    ak8jetCSV.clear()
    ak8jetSDmass.clear()
    ak8jetSDsubjet0pt.clear()
    ak8jetSDsubjet0eta.clear()
    ak8jetSDsubjet0phi.clear()
    ak8jetSDsubjet0mass.clear()
    ak8jetSDsubjet0CSV.clear()
    ak8jetSDsubjet1pt.clear()
    ak8jetSDsubjet1eta.clear()
    ak8jetSDsubjet1phi.clear()
    ak8jetSDsubjet1mass.clear()
    ak8jetSDsubjet1CSV.clear()
    ak8jetJECunc.clear()
    eventWeight_nom.clear()
    if options.isMC:
        ak8jetPtJERup.clear()
        ak8jetPtJERdown.clear()
        ak8jetEtaJERup.clear()
        ak8jetEtaJERdown.clear()
        ak8jetPhiJERup.clear()
        ak8jetPhiJERdown.clear()
        ak8jetMassJERup.clear()
        ak8jetMassJERdown.clear()
        eventWeight_puUp.clear()
        eventWeight_puDown.clear()
        if options.semilep is not None:
            eventWeight_Q2Up.clear()
            eventWeight_Q2Down.clear()
            eventWeight_PDF.clear()
            eventWeight_alphaUp.clear()
            eventWeight_alphaDown.clear()
        if options.semilep == 1:
            truthChannel.clear()

    weight_nom = 1.0 #event weight
    weight_puUp = 1.0
    weight_puDown = 1.0
    weight_Q2up = 1.0
    weight_Q2down = 1.0
    weight_PDF = 1.0
    weight_alphaUp = 1.0
    weight_alphaDown = 1.0

    mttbarGen = -1.0
    
    if ntotal % 1000 == 0 :
      print  '--------- Processing Event ' + str(ntotal)
    ntotal += 1

    if options.debug :
        if ntotal > 1000:
            continue
        print 'Event ' + str(ntotal)

    event.getByLabel(runnumLabel,runnumHandle)
    runnumber = runnumHandle.product()[0]

    # Discard events in wrong IOV
    if not options.isMC:
        if options.debug :
            print 'Run number is ' + str(runnumber) + ', JEC IOV is ' + options.jecIOV
            if "BCD" in options.jecIOV:
                print options.jecIOV + ' = BCD'
            if runnumber == 273425 :
                print str(runnumber) + ' = 273425'

        if not (("BCD" in options.jecIOV and 1 <= runnumber <= 276811) or 
                ("EF" in options.jecIOV and 276831 <= runnumber <= 278801) or
                ("G" in options.jecIOV and 278802 <= runnumber <= 280385) or
                ("H" in options.jecIOV and 280919 <= runnumber)):
            continue

    nInIOV += 1

    # -------------------------------------------------------------------------------------
    # Get parton-level info
    # -------------------------------------------------------------------------------------

    if options.isMC and options.semilep is not None :

        genMuons = []
        genElectrons = []
        isSemiLeptonicGen = True
        isMuon = False
        isElectron = False
        leptonicTop = 0 #start by assuming hadronic decay
        leptonicAntitop = 0
        
        event.getByLabel( genParticlesLabel, genParticlesHandle )
                
        if genParticlesHandle.isValid() == False :
            continue
        if options.debug:
            print 'Got gen particles!'

        genParticles  = genParticlesHandle.product()
        
        p4Top = ROOT.TLorentzVector()
        p4Antitop = ROOT.TLorentzVector()
        p4Lepton = ROOT.TLorentzVector()
        p4HadTop = ROOT.TLorentzVector()
        p4LepTop = ROOT.TLorentzVector()
        nPromptLepton = 0

        if options.debug :
            print 'There are ' + str(len(genParticles)) + ' particles.'

        # loop over gen particles
        for genParticle in genParticles :
            
            if options.debug and genParticle.isHardProcess():
                print 'Gen particle with ID ' + str(genParticle.pdgId()) + ' and status ' + str(genParticle.status())
            # Find tops -- |pdgID| = 6, isLastCopy
            if genParticle.pdgId() == 6 and genParticle.isLastCopy() :
                p4Top.SetPtEtaPhiM( genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.mass() )                    
            if genParticle.pdgId() == -6 and genParticle.isLastCopy():
                p4Antitop.SetPtEtaPhiM( genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.mass() )
            
            # Find prompt leptons
            if (abs(genParticle.pdgId()) == 11 or abs(genParticle.pdgId()) == 13) and genParticle.isPromptFinalState():
                nPromptLepton += 1
                if options.debug:
                    print 'Found prompt mu/el'
                if genParticle.pdgId() > 0:
                    leptonicAntitop = 1 #Prompt lepton comes from W-, aka antitop
                else :
                    leptonicTop = 1 #Prompt antilepton comes from W+, aka top
            if abs(genParticle.pdgId()) == 15 and genParticle.isPromptDecayed() :
                nPromptLepton += 1
                if options.debug:
                    print 'Found prompt tau'
                if genParticle.pdgId() > 0:
                    leptonicAntitop = 1
                else :
                    leptonicTop = 1

            # Find good muons / electrons
            if abs(genParticle.pdgId()) == 11 : 
                if genParticle.isPromptFinalState() or genParticle.isDirectPromptTauDecayProductFinalState():
                    if options.debug and genParticle.isDirectPromptTauDecayProductFinalState():
                        print 'Found el from tau'
                    isElectron = True
                    p4Electron = ROOT.TLorentzVector()
                    p4Electron.SetPtEtaPhiM( genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.mass() )
                    genElectrons.append(p4Electron)

            if abs(genParticle.pdgId()) == 13 :
                if genParticle.isPromptFinalState() or genParticle.isDirectPromptTauDecayProductFinalState():
                    if options.debug and genParticle.isDirectPromptTauDecayProductFinalState() :
                        print 'Found mu from tau'
                    isMuon = True
                    p4Muon = ROOT.TLorentzVector()
                    p4Muon.SetPtEtaPhiM( genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.mass() )
                    genMuons.append(p4Muon)
                
        if (nPromptLepton == 1) and (isMuon == True) and (isElectron == False) :
            p4Lepton = genMuons[0]
            nMuChParton += 1
            if options.debug :
                print 'Event in muon channel'
            isSemiLeptonicGen = True
        elif (nPromptLepton == 1) and (isMuon == False) and (isElectron == True) :
            p4Lepton = genElectrons[0]
            nElChParton += 1
            if options.debug :
                print 'Event in electron channel'
            isSemiLeptonicGen = True
        else :
            if nPromptLepton == 0:
                nHadParton += 1
            if nPromptLepton == 2:
                nDilepParton += 1
            if nPromptLepton == 1 and (isMuon == False) and (isElectron == False) :
                nTauOnlyParton += 1
            if options.debug :
                print 'Event not semileptonic'
                print 'nPromptLepton = ' + str(nPromptLepton)
                print 'isMuon = ' + str(isMuon)
                print 'isElectron = ' + str(isElectron)
            isSemiLeptonicGen = False

        if options.semilep > 0 and isSemiLeptonicGen == False:
            continue
        if options.semilep < 0 and isSemiLeptonicGen == True:
            continue

        # Store gen information for signal
        if options.semilep == 1:

            # Get truth channel
            if isMuon:
                truthChannel.push_back(0)
            elif isElectron:
                truthChannel.push_back(1)
            else:
                truthChannel.push_back(2)

            ttbarGen = p4Top + p4Antitop

            # Get hadronic top
            if not leptonicTop + leptonicAntitop == 1 :
                print 'Error: decay chain says event is not semileptonic'
            if leptonicTop :
                if options.debug :
                    print 'Top decays leptonically!'
                p4HadTop = p4Antitop
                p4LepTop = p4Top
            else :
                if options.debug :
                    print 'Antitop decays leptonically!'
                p4HadTop = p4Top
                p4LepTop = p4Antitop

            #Check assignment is correct by comparing to lepton
            if p4HadTop.DeltaR(p4Lepton) < p4LepTop.DeltaR(p4Lepton) :
                nWrongHemisphere += 1
                if options.debug :
                    print 'Error: hadronic top is closer to lepton! dR(t_had,l) = ' + str(p4HadTop.DeltaR(p4Lepton)) + ' < dR(t_lep,l) = ' + str(p4LepTop.DeltaR(p4Lepton))

            # ------------------------------------------------------------
            # Store parton-level information 
            # ------------------------------------------------------------
        
            nPassParton += 1
            genTopPt.push_back(p4HadTop.Perp())
            genTopEta.push_back(p4HadTop.Eta())
            genTopPhi.push_back(p4HadTop.Phi())
            genTTbarMass.push_back(ttbarGen.M())
        
            if len(genMuons) != 0:
                genMuPt.push_back(genMuons[0].Perp())
                genMuEta.push_back(genMuons[0].Eta())
                genMuPhi.push_back(genMuons[0].Phi())
            
            if len(genElectrons) != 0:
                genElPt.push_back(genElectrons[0].Perp())
                genElEta.push_back(genElectrons[0].Eta())
                genElPhi.push_back(genElectrons[0].Phi())   

            # -------------------------------------------------------------------------------------
            # read gen jets
            # -------------------------------------------------------------------------------------
        
            # Get AK8 gen jet information
            genTops = []
            event.getByLabel( ak8GenJetPtLabel, ak8GenJetPtHandle )
            if ak8GenJetPtHandle.isValid() :
                event.getByLabel( ak8GenJetEtaLabel, ak8GenJetEtaHandle )
                event.getByLabel( ak8GenJetPhiLabel, ak8GenJetPhiHandle )
                event.getByLabel( ak8GenJetELabel, ak8GenJetEHandle )
            
                ak8GenJetPt   = ak8GenJetPtHandle.product()
                ak8GenJetEta  = ak8GenJetEtaHandle.product()
                ak8GenJetPhi  = ak8GenJetPhiHandle.product()
                ak8GenJetE    = ak8GenJetEHandle.product()
                
                # loop over AK8 gen jets
                if len(ak8GenJetPt) != 0:
                    for iak8 in xrange( len(ak8GenJetPt) ) :
                        p4 = ROOT.TLorentzVector()
                        p4.SetPtEtaPhiE( ak8GenJetPt[iak8], ak8GenJetEta[iak8], ak8GenJetPhi[iak8], ak8GenJetE[iak8] )
                        if abs(ak8GenJetEta[iak8]) < MAX_JET_ETA and p4.M() > MIN_TOP_MASS and p4.M() < MAX_TOP_MASS: # no pt cut                              
                            genTops.append(p4)

            # -------------------------------------------------------------------------------------
            # Get particle-level leptons
            # -------------------------------------------------------------------------------------

            partMu = []
            for iMuon in genMuons:
                if iMuon.Perp() > MIN_MU_PT and abs(iMuon.Eta()) < MAX_MU_ETA:  ## pt>50, |eta|<2.1
                    partMu.append(iMuon)
                
            partEl = []
            for iEle in genElectrons:
                if iEle.Perp() > MIN_MU_PT and abs(iEle.Eta()) < MAX_MU_ETA:  ## pt>50, |eta|<2.1  (same selection as for muons here!)
                    partEl.append(iEle)

            # -------------------------------------------------------------------------------
            # Store particle-level info if event good at particle level
            # -------------------------------------------------------------------------------

            if len(genTops) >= 1 and (len(partMu) + len(partEl)) >= 1:
                nPassParticle += 1

                for itjet in genTops :
                    genAK8jetPt.push_back(itjet.Perp())
                    genAK8jetEta.push_back(itjet.Eta())
                    genAK8jetPhi.push_back(itjet.Phi())
                    genAK8jetMass.push_back(itjet.M())
                
                for imu in partMu :
                    partMuPt.push_back(imu.Perp())
                    partMuEta.push_back(imu.Eta())
                    partMuPhi.push_back(imu.Phi())
                
                for iel in partEl :
                    partElPt.push_back(iel.Perp())
                    partElEta.push_back(iel.Eta())
                    partElPhi.push_back(iel.Phi())

    # -------------------
    #
    # R E C O   L E V E L
    #
    # -------------------

    passReco = True #This will be set to false later if cuts fail

    # ------------------------------
    # Require a good primary vertex
    # ------------------------------

    event.getByLabel( pvChiLabel, pvChiHandle )
    event.getByLabel( pvRhoLabel, pvRhoHandle )
    event.getByLabel( pvZLabel, pvZHandle )
    event.getByLabel( pvNdofLabel, pvNdofHandle )
    
    pv_chi  = pvChiHandle.product()
    pv_rho  = pvRhoHandle.product()
    pv_z    = pvZHandle.product()
    pv_ndof = pvNdofHandle.product()
    NPV = 0
    
    for ivtx in xrange( len(pv_chi) ) :
        if abs(pv_z[ivtx]) < 24. and pv_ndof[ivtx] > 4 and abs(pv_rho[ivtx]) < 2.0 :
            NPV += 1

    # -------------------------------
    # Do pileup reweighting if MC
    # -------------------------------
    
    if options.isMC :
        event.getByLabel(puNtrueIntLabel, puNtrueIntHandle)
        puNTrueInt = puNtrueIntHandle.product()[0] 
        h_NtrueIntPU.Fill( puNTrueInt )

        if options.puFile is not None :
            weight_nom  *= hPUweight_nom.GetBinContent( hPUweight_nom.GetXaxis().FindBin( puNTrueInt ) )
            weight_puUp   *= hPUweight_up.GetBinContent( hPUweight_up.GetXaxis().FindBin( puNTrueInt ) )
            weight_puDown *= hPUweight_down.GetBinContent( hPUweight_down.GetXaxis().FindBin( puNTrueInt ) )
            if options.debug :
                print 'Event weight after pileup reweigting is : ' + str(weight_nom)

    #Fill NPV hists before and after PU reweighting
    h_NPV_noweight.Fill(NPV)
    h_NPV.Fill(NPV,weight_nom)

    if NPV == 0 :
        passReco = False
        if not options.fullTruth:
            continue
    else :
        nPassNPV += 1

    # -------------------------------
    # Get Q2 and NNPDF weights
    # -------------------------------
    if options.isMC and options.semilep is not None:
        gotLHE = event.getByLabel( LHELabel, LHEHandle )
        if not gotLHE:
            passReco = False
            if not options.fullTruth:
                continue
        else :
            LHEEventInfo = LHEHandle.product()
           
            if len(LHEEventInfo.weights()) < 111: #Some events have an incomplete LHE weight vector...
                passReco = False
                if not options.fullTruth:
                    continue
            else :
                nPassLHE += 1

                #Q^2 up and down
                for i_lhe in range(0,9):
                    if i_lhe != 5 and i_lhe != 7:
                        Q2wgt = LHEEventInfo.weights()[i_lhe].wgt
                        Q2wgt_frac = Q2wgt/(LHEEventInfo.weights()[0].wgt)
                        weight_Q2up = max(weight_Q2up, Q2wgt_frac)
                        weight_Q2down = min(weight_Q2down, Q2wgt_frac)

                #NNPDF3 PDF (symmetric)
                NNPDF3wgtAvg = 0.0
                NNPDF3wgtRMS = 0.0
                PDFcentral = 1.0
                PDFstart = 9
                PDFend = 109

                for i_lhePDF in range(PDFstart,PDFend):
                    NNPDF3wgt = LHEEventInfo.weights()[i_lhePDF].wgt
                    NNPDF3wgt_frac = NNPDF3wgt/(PDFcentral)
                    NNPDF3wgtAvg += NNPDF3wgt_frac
        
                NNPDF3wgtAvg = NNPDF3wgtAvg/100.0
        
                for i_lhePDF in range(PDFstart,PDFend):
                    NNPDF3wgt = LHEEventInfo.weights()[i_lhePDF].wgt
                    NNPDF3wgt_frac = NNPDF3wgt/(PDFcentral)
                    NNPDF3wgtRMS += (NNPDF3wgt_frac - NNPDF3wgtAvg)*(NNPDF3wgt_frac - NNPDF3wgtAvg)
 
                weight_PDF = math.sqrt(NNPDF3wgtRMS/99.0)

                weight_alphaUp = LHEEventInfo.weights()[109].wgt
                weight_alphaDown = LHEEventInfo.weights()[110].wgt

    # -------------------------------
    # Require met filter
    # -------------------------------
    
    metFilt = True
    gotBits = event.getByLabel( metFilterBitsLabel, metFilterBitsHandle )

    if gotBits == False  :
        print 'No MET filter bits!'
        metFilt = False
        passReco = False
        if not options.fullTruth:
            continue
    else :
        filterBits = metFilterBitsHandle.product()

        metFilterIndices = metFiltIndices[runnumber]
        for ifilt in metFilterIndices :
            if filterBits[ifilt] == 0:
                if options.debug :
                    print 'MET filter fails'
                metFilt = False

    # Now check additional MET filters
    event.getByLabel(BadChargedCandidateFilterLabel, BadChargedCandidateFilterHandle)
    badChargedCandidateFilter = BadChargedCandidateFilterHandle.product()
    event.getByLabel(BadPFMuonFilterLabel, BadPFMuonFilterHandle)
    badPFMuonFilter = BadPFMuonFilterHandle.product()
    
    if not (badChargedCandidateFilter and badPFMuonFilter):
        metFilt = False

    if metFilt == False :
        passReco = False
        if not options.fullTruth:
            continue
    else :
        nPassMetFilter += 1
                
    # -------------------------------
    # Require a trigger
    # -------------------------------

    passMuTrig = False
    passElTrig = False
    prescale = 1.0
    
    event.getByLabel( trigBitsLabel, trigBitsHandle )
    event.getByLabel( trigPrescalesLabel, trigPrescalesHandle )
    
    triggerBits = trigBitsHandle.product()
    triggerPrescales = trigPrescalesHandle.product()
    
    if triggerBits[muTrigIndices[runnumber]] == 1 :
        passMuTrig = True
        prescale = prescale * triggerPrescales[muTrigIndices[runnumber]]
        h_muPrescale.Fill(triggerPrescales[muTrigIndices[runnumber]])

    if triggerBits[elTrigIndices[runnumber]] == 1 :
        passElTrig = True
        prescale = prescale * triggerPrescales[elTrigIndices[runnumber]]
        h_elPrescale.Fill(triggerPrescales[elTrigIndices[runnumber]])
            
    weight_nom = weight_nom * prescale #Currently prescale has both mu and el prescales if both triggers fired
    weight_puUp = weight_puUp * prescale 
    weight_puDown = weight_puDown * prescale
        
    if not options.isMC :
        if "mu" in options.muOrEl and not passMuTrig:
            continue
        if "el" in options.muOrEl and not passElTrig:
            continue
    
    # -------------------------------------------------------------------------------------
    # read event rho value
    # -------------------------------------------------------------------------------------

    event.getByLabel( rhoLabel, rhoHandle )
    if len(rhoHandle.product()) == 0 :
        passReco = False
        print "Event has no rho values."
        if not options.fullTruth :
            continue
    else :
        rho = rhoHandle.product()[0]
        nPassRho += 1
    
    # -------------------------------------------------------------------------------------
    # get electrons
    # -------------------------------------------------------------------------------------

    elCand = []
    elCandKey = []
    lepCand = []
    lepCandKey = []

    event.getByLabel (electronPtLabel, electronPtHandle)
    if electronPtHandle.isValid() :
        electronPts = electronPtHandle.product()
        event.getByLabel (electronEtaLabel, electronEtaHandle)
        electronEtas = electronEtaHandle.product()
        event.getByLabel (electronPhiLabel, electronPhiHandle)
        electronPhis = electronPhiHandle.product()
        event.getByLabel (electronChargeLabel, electronChargeHandle)
        electronCharges = electronChargeHandle.product()
        event.getByLabel (electronTightLabel, electronTightHandle)
        isElectronTight = electronTightHandle.product()
        event.getByLabel (electronMediumLabel, electronMediumHandle)
        isElectronMedium = electronMediumHandle.product()
        event.getByLabel (electronLooseLabel, electronLooseHandle)
        isElectronLoose = electronLooseHandle.product()
        event.getByLabel (electronfull5x5sieeLabel, electronfull5x5sieeHandle)
        electronFull5x5siees = electronfull5x5sieeHandle.product()
        event.getByLabel (electrondEtaInLabel, electrondEtaInHandle)
        electronDEtaIns = electrondEtaInHandle.product()
        event.getByLabel (electrondPhiInLabel, electrondPhiInHandle)
        electronDPhiIns = electrondPhiInHandle.product()
        event.getByLabel (electronHoELabel, electronHoEHandle)
        electronHoEs = electronHoEHandle.product()
        event.getByLabel (electronooEmooPLabel, electronooEmooPHandle)
        electronooEmooPs = electronooEmooPHandle.product()
        event.getByLabel (electronMiniIsoLabel, electronMiniIsoHandle)
        electronMiniIsos = electronMiniIsoHandle.product()
        event.getByLabel (electronSCEtaLabel, electronSCEtaHandle)
        electronSCEtas = electronSCEtaHandle.product()
        event.getByLabel (electronMissHitsLabel, electronMissHitsHandle)
        electronMissingHits = electronMissHitsHandle.product()
        event.getByLabel (electronConVetoLabel, electronConVetoHandle)
        isElectronConVeto = electronConVetoHandle.product()
        event.getByLabel (electronMediumNoIsoLabel, electronMediumNoIsoHandle)
        isElectronMediumNoIso = electronMediumNoIsoHandle.product()
        event.getByLabel (electronTightNoIsoLabel, electronTightNoIsoHandle)
        isElectronTightNoIso = electronTightNoIsoHandle.product()
        event.getByLabel (electronDxyLabel, electronDxyHandle)
        electronDxys = electronDxyHandle.product()
        event.getByLabel (electronDzLabel, electronDzHandle)
        electronDzs = electronDzHandle.product()

        event.getByLabel (elKeyLabel, elKeyHandle)
        elKeys = elKeyHandle.product()

        nTrueMedium = 0
        nFakeMedium = 0
        nTrueTight = 0
        nFakeTight = 0
        for ielectronPt in range(0,len(electronPts)) :
            electronPt = electronPts[ielectronPt]
            electronEta = electronEtas[ielectronPt]
            electronPhi = electronPhis[ielectronPt]
            electronMass = 0.0
            if (electronPt < MIN_EL_PT or abs(electronEta) > MAX_EL_ETA ) :
                continue
            nEleRaw += 1.0
            manualEisMedium = False
            if abs( electronSCEtas[ielectronPt] ) <= 1.479 :
                if abs(electronDEtaIns[ielectronPt]) < 0.00311 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.103 :
                        if electronFull5x5siees[ielectronPt] < 0.00998 :
                            if electronHoEs[ielectronPt] <  0.253 :
                                if electronooEmooPs[ielectronPt] <  0.134 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        #if not isElectronConVeto[ielectronPt] :
                                        manualEisMedium = True
            else :
                if abs(electronDEtaIns[ielectronPt]) < 0.00609 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.045 :
                        if electronFull5x5siees[ielectronPt] < 0.0298 :
                            if electronHoEs[ielectronPt] <  0.0878 :
                                if electronooEmooPs[ielectronPt] <  0.13 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        #if not isElectronConVeto[ielectronPt] :
                                        manualEisMedium = True
                                            
            manualEisTight = False
            if abs( electronSCEtas[ielectronPt] ) <= 1.479 :
                if abs(electronDEtaIns[ielectronPt]) < 0.00308 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.0816 :
                        if electronFull5x5siees[ielectronPt] < 0.00998 :
                            if electronHoEs[ielectronPt] <  0.0414 :
                                if electronooEmooPs[ielectronPt] <  0.0129 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        #if not isElectronConVeto[ielectronPt]:
                                        manualEisTight = True
            else :
                if abs(electronDEtaIns[ielectronPt]) < 0.00605 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.0394 :
                        if electronFull5x5siees[ielectronPt] < 0.0292 :
                            if electronHoEs[ielectronPt] <  0.0641 :
                                if electronooEmooPs[ielectronPt] <  0.0129 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        #if not isElectronConVeto[ielectronPt] :
                                        manualEisTight = True
                                            
            #eIsTight = isElectronTight[ielectronPt]
            #eIsMedium = isElectronMedium[ielectronPt]
            #eIsTightNoIso = isElectronTightNoIso[ielectronPt]
            #eIsMediumNoIso = isElectronMediumNoIso[ielectronPt]
            
            #if eIsMedium and not manualEisMedium :
            #    nMediumNotManualEl += 1.0
            #if manualEisMedium and not eIsMedium :
            #    nManualNotMediumEl += 1.0
            #if eIsTight and not manualEisTight :
            #    nTightNotManualEl += 1.0
            #if manualEisTight and not eIsTight :
            #    nManualNotTightEl += 1.0

            #if manualEisMedium and not eIsMediumNoIso :
            #    print 'Medium manual cuts do not match no iso'
            #    print 'electron abs(eta): ' + str(abs(electronSCEtas[ielectronPt]))
            #    print 'full5x5_sigmaIetaIeta: ' + str(electronFull5x5siees[ielectronPt])
            #    print 'abs(dEtaInSeed): ' + str(abs(electronDEtaIns[ielectronPt]))
            #    print 'abs(dPhiIn): ' + str(abs(electronDPhiIns[ielectronPt]))
            #    print 'H/E: '+ str(electronHoEs[ielectronPt])
            #    print 'abs(1/E-1/p): ' + str(electronooEmooPs[ielectronPt])
            #    print 'expected missing inner hits: ' + str(electronMissingHits[ielectronPt])
            #    nManualNotMediumNoIsoEl += 1.0
            #if manualEisTight and not eIsTightNoIso :
            #    print 'Tight manual cuts do not match no iso'
            #    print 'electron abs(eta): ' + str(abs(electronSCEtas[ielectronPt]))
            #    print 'full5x5_sigmaIetaIeta: ' + str(electronFull5x5siees[ielectronPt])
            #    print 'abs(dEtaInSeed): ' + str(abs(electronDEtaIns[ielectronPt]))
            #    print 'abs(dPhiIn): ' + str(abs(electronDPhiIns[ielectronPt]))
            #    print 'H/E: '+ str(electronHoEs[ielectronPt])
            #    print 'abs(1/E-1/p): ' + str(electronooEmooPs[ielectronPt])
            #    print 'expected missing inner hits: ' + str(electronMissingHits[ielectronPt])
            #    nManualNotTightNoIsoEl += 1.0
            #if eIsMediumNoIso and not manualEisMedium :
            #    nMediumNoIsoNotManual += 1.0
            #if eIsTightNoIso and not manualEisTight :
            #    nTightNoIsoNotManual += 1.0

            # Additional d0 / dz cuts
            passD0 = False
            passDz = False 
            if abs( electronSCEtas[ielectronPt] ) <= 1.479 :
                if electronDxys[ielectronPt] < 0.05 :
                    passD0 = True
                if electronDzs[ielectronPt] < 0.10 :
                    passDz = True
            else :
                if electronDxys[ielectronPt] < 0.10 :
                    passD0 = True
                if electronDzs[ielectronPt] < 0.20 :
                    passDz = True

            if manualEisMedium and passD0 and passDz and isElectronConVeto[ielectronPt] :
                nTrueMedium += 1.0
            if manualEisMedium and passD0 and passDz and not isElectronConVeto[ielectronPt] :
                nFakeMedium += 1.0
            if manualEisTight and passD0 and passDz and isElectronConVeto[ielectronPt] :
                nTrueTight += 1.0
            if manualEisTight and passD0 and passDz and not isElectronConVeto[ielectronPt] :
                nFakeTight += 1.0

            if not (manualEisMedium and passD0 and passDz) :
                continue

            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( electronPt, electronEta, electronPhi, electronMass )

            elPt.push_back(electronPt)
            elEta.push_back(electronEta)
            elPhi.push_back(electronPhi)
            elCharge.push_back(electronCharges[ielectronPt])
            elMiniIso.push_back(electronMiniIsos[ielectronPt])

            if isElectronConVeto[ielectronPt] : #Electrons for cleaning pass full Medium ID
                elCand.append(p4)
                lepCand.append(p4)
                elCandKey.append(elKeys[ielectronPt])
                lepCandKey.append(elKeys[ielectronPt])
                
            if manualEisTight :
                elTight.push_back(1)
            else :
                elTight.push_back(0)
                
            if isElectronConVeto[ielectronPt] :
                elPassConVeto.push_back(1)
            else :
                elPassConVeto.push_back(0)
                
    # --------------------------
    # get muons
    # --------------------------

    muCand = []
    muCandKey = []

    event.getByLabel (muonPtLabel, muonPtHandle)
    if muonPtHandle.isValid() : 
        muonPts = muonPtHandle.product()
        event.getByLabel (muonEtaLabel, muonEtaHandle)
        muonEtas = muonEtaHandle.product()
        event.getByLabel (muonPhiLabel, muonPhiHandle)
        muonPhis = muonPhiHandle.product()
        event.getByLabel (muonChargeLabel, muonChargeHandle)
        muonCharges = muonChargeHandle.product()
        event.getByLabel (muTightLabel, muTightHandle)
        isTightMuon = muTightHandle.product()
        event.getByLabel (muMediumLabel, muMediumHandle)
        isMediumMuon = muMediumHandle.product()
        event.getByLabel (muLooseLabel, muLooseHandle)
        isLooseMuon = muLooseHandle.product()
        event.getByLabel (muPFLabel, muPFHandle)
        isPFMuon = muPFHandle.product()
        event.getByLabel (muGlobalLabel, muGlobalHandle)
        isGlobalMuon = muGlobalHandle.product()
        event.getByLabel (muTrackerLabel, muTrackerHandle)
        isTrackerMuon = muTrackerHandle.product()
        event.getByLabel (muTrkChi2Label, muTrkChi2Handle)
        muTrkChi2s = muTrkChi2Handle.product()
        event.getByLabel (muHitLabel, muHitHandle)
        muHits = muHitHandle.product()
        event.getByLabel (muMatchStatLabel, muMatchStatHandle)
        muMatchStats = muMatchStatHandle.product()
        event.getByLabel (muPixHitLabel, muPixHitHandle)
        muPixHits = muPixHitHandle.product()
        event.getByLabel (muTrkLayerLabel, muTrkLayerHandle)
        muTrkLayers = muTrkLayerHandle.product()
        event.getByLabel (muDxyLabel, muDxyHandle)
        muDxys = muDxyHandle.product()
        event.getByLabel (muDzLabel, muDzHandle)
        muDzs = muDzHandle.product()
        event.getByLabel (muMiniIsoLabel, muMiniIsoHandle)
        muMiniIsos = muMiniIsoHandle.product()

        event.getByLabel (muKeyLabel, muKeyHandle)
        muKeys = muKeyHandle.product()
        
        for imuonPt in range(0,len(muonPts)) :
            muonPt = muonPts[imuonPt]
            muonEta = muonEtas[imuonPt]
            muonPhi = muonPhis[imuonPt]
            muonMass = 0.105658
            if muonPt < MIN_MU_PT or abs(muonEta) > MAX_MU_ETA :
                continue
            nMuRaw += 1.0
            manualMuIsTight = False
            if isPFMuon[imuonPt] :
                if isGlobalMuon[imuonPt] :
                    if muTrkChi2s[imuonPt] < 10 :
                        if muHits[imuonPt] > 0. :
                            if muMatchStats[imuonPt] > 1.0 :
                                if muDxys[imuonPt] < 0.2 : #2 mm
                                    if muDzs[imuonPt] < 0.5 : #5 mm
                                        if muPixHits[imuonPt] > 0. :
                                            if muTrkLayers[imuonPt] > 5.0 :
                                                manualMuIsTight = True
            muIsTight = isTightMuon[imuonPt]
            muIsMedium = isMediumMuon[imuonPt]
            if muIsTight and not manualMuIsTight :
                nTightNotManualMu += 1.0
            if manualMuIsTight and not muIsTight :
                nManualNotTightMu += 1.0

            if not muIsMedium :
                continue
            
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( muonPt, muonEta, muonPhi, muonMass )
            
            muPt.push_back(muonPt)
            muEta.push_back(muonEta)
            muPhi.push_back(muonPhi)
            muCharge.push_back(muonCharges[imuonPt])
            muMiniIso.push_back(muMiniIsos[imuonPt])

            muCand.append(p4)
            lepCand.append(p4)
            muCandKey.append(muKeys[imuonPt])
            lepCandKey.append(muKeys[imuonPt])

            if muIsTight :
                muTight.push_back(1)
            else :
                muTight.push_back(0)
                
    # -------------------------------------------------------------------------------------
    # check that we have exactly one lepton candidate
    # -------------------------------------------------------------------------------------

    if not ((passMuTrig and len(muCand) == 1 and (nTrueMedium + nFakeMedium == 0)) or (passElTrig and (nTrueMedium == 1 or (nTrueMedium == 0 and nFakeMedium >= 1)) and len(muCand) == 0)):
        passReco = False
        if not options.fullTruth :
            continue
    else :
        nPassLep += 1

    # -------------------------------------------------------------------------------------
    # read MET 
    # -------------------------------------------------------------------------------------

    event.getByLabel (metLabel, metHandle)
    mets = metHandle.product()
    metRaw = mets[0]
    event.getByLabel (metphiLabel, metphiHandle)
    metphis = metphiHandle.product()
    metphi = metphis[0]
    met_px = metRaw * math.cos( metphi )
    met_py = metRaw * math.sin( metphi )            
    met = math.sqrt(met_px*met_px + met_py*met_py)
    metv = ROOT.TLorentzVector()
    metv.SetPtEtaPhiM( met, 0.0, metphi, 0.0)
            
    # -------------------------------------------------------------------------------------
    # read AK4 jet information
    # -------------------------------------------------------------------------------------

    event.getByLabel (ak4JetPtLabel, ak4JetPtHandle)
    if ak4JetPtHandle.isValid(): 
        ak4JetPts = ak4JetPtHandle.product()
        event.getByLabel (ak4JetEtaLabel, ak4JetEtaHandle)
        ak4JetEtas = ak4JetEtaHandle.product()
        event.getByLabel (ak4JetPhiLabel, ak4JetPhiHandle)
        ak4JetPhis = ak4JetPhiHandle.product()
        event.getByLabel (ak4JetELabel, ak4JetEHandle)
        ak4JetEs = ak4JetEHandle.product()
        event.getByLabel (ak4JetCSVLabel, ak4JetCSVHandle)
        ak4JetCSVs = ak4JetCSVHandle.product()
        event.getByLabel (ak4JetVtxMassLabel, ak4JetVtxMassHandle)
        ak4JetVtxMasses = ak4JetVtxMassHandle.product()
        event.getByLabel (ak4JetAreaLabel, ak4JetAreaHandle)
        ak4JetAreas = ak4JetAreaHandle.product()

        # variables needed for jet ID 
        event.getByLabel( ak4JetNeuHadEnergyFracLabel, ak4JetNeuHadEnergyFracHandle )
        event.getByLabel( ak4JetNeuEmEnergyFracLabel, ak4JetNeuEmEnergyFracHandle )
        event.getByLabel( ak4JetChHadEnergyFracLabel, ak4JetChHadEnergyFracHandle )
        event.getByLabel( ak4JetChEmEnergyFracLabel, ak4JetChEmEnergyFracHandle )
        event.getByLabel( ak4JetNeuMultiLabel, ak4JetNeuMultiHandle )
        event.getByLabel( ak4JetChMultiLabel, ak4JetChMultiHandle )

        ak4JetNeuHadEnergyFracs = ak4JetNeuHadEnergyFracHandle.product()
        ak4JetNeuEmEnergyFracs = ak4JetNeuEmEnergyFracHandle.product()
        ak4JetChHadEnergyFracs = ak4JetChHadEnergyFracHandle.product()
        ak4JetChEmEnergyFracs = ak4JetChEmEnergyFracHandle.product()
        ak4JetNeuMultis = ak4JetNeuMultiHandle.product()
        ak4JetChMultis = ak4JetChMultiHandle.product()

        # applied JEC
        event.getByLabel( ak4JetJECLabel, ak4JetJECHandle )
        ak4JetJECs = ak4JetJECHandle.product()

        #JER
        event.getByLabel(ak4MatchedGenJetPtLabel, ak4MatchedGenJetPtHandle)
        ak4MatchedGenJetPts = ak4MatchedGenJetPtHandle.product()
        event.getByLabel(ak4MatchedGenJetEtaLabel, ak4MatchedGenJetEtaHandle)
        ak4MatchedGenJetEtas = ak4MatchedGenJetEtaHandle.product()
        event.getByLabel(ak4MatchedGenJetPhiLabel, ak4MatchedGenJetPhiHandle)
        ak4MatchedGenJetPhis = ak4MatchedGenJetPhiHandle.product()
        event.getByLabel(ak4MatchedGenJetEnergyLabel, ak4MatchedGenJetEnergyHandle)
        ak4MatchedGenJetEnergys = ak4MatchedGenJetEnergyHandle.product()
        event.getByLabel(ak4JERSFnomLabel, ak4JERSFnomHandle)
        ak4JERSFnoms = ak4JERSFnomHandle.product()
        event.getByLabel(ak4JERSFupLabel, ak4JERSFupHandle)
        ak4JERSFups = ak4JERSFupHandle.product()
        event.getByLabel(ak4JERSFdownLabel, ak4JERSFdownHandle)
        ak4JERSFdowns = ak4JERSFdownHandle.product()
        event.getByLabel(ak4JERLabel, ak4JERHandle)
        ak4JERs = ak4JERHandle.product()
        event.getByLabel(ak4JetHadronFlavourLabel, ak4JetHadronFlavourHandle)
        ak4jetHadronFlavours = ak4JetHadronFlavourHandle.product()

        # jet constituents 
        event.getByLabel( ak4JetKeysLabel, ak4JetKeysHandle )
        ak4JetKeys = ak4JetKeysHandle.product()
        
        jetsFor2DPt15 = []
        jetsFor2DPt20 = []
        jetsFor2DPt25 = []
        jetsFor2DPt30 = []
        jetsFor2DPt35 = []
        jetsFor2DPt40 = []
        jetsFor2DPt45 = []

        HT = 0.0

        # -------------------------------------------------------------------------------------
        # loop over AK4 jets
        # -------------------------------------------------------------------------------------

        for ijet in xrange( len(ak4JetPts) ) :

            # jet IDs must be calculated prior to JECs, so remove them 
            jetP4Pre = ROOT.TLorentzVector()
            jetP4Pre.SetPtEtaPhiE( ak4JetPts[ijet], ak4JetEtas[ijet], ak4JetPhis[ijet], ak4JetEs[ijet] )

            # get the raw jet energy
            jetP4Raw = jetP4Pre * ak4JetJECs[ijet]

            # calculate the neutral/charged hadron/em energy fractions
            nhf = ak4JetNeuHadEnergyFracs[ijet]
            nef = ak4JetNeuEmEnergyFracs[ijet]
            chf = ak4JetChHadEnergyFracs[ijet]
            cef = ak4JetChEmEnergyFracs[ijet]
            nconstituents = ak4JetNeuMultis[ijet] + ak4JetChMultis[ijet]
            nchmult = ak4JetChMultis[ijet] 

            # require loose jet ID
            if (nhf >= 0.99 or nef >= 0.99 or chf <= 0.00 or cef >= 0.99 or nconstituents <= 1 or nchmult <= 0) :
                continue

            # clean leptons from jets 
            zeroedEnergy = False

            cleanedLepton = 0
            for ilep in xrange( len(lepCand) ) :
                if lepCand[ilep].DeltaR(jetP4Raw) < 0.4 :
                    # check jet daughters close to the lepton
                    #pfcands = int(nconstituents)
                    #for ipf in range(0, pfcands) :
                    ak4daughters = ak4JetKeys[ijet]
                    for daughterKey in ak4daughters :
                        # if jet daughter matches lepton, remove lepton p4 from (raw) jet p4
                        #if ak4JetKeys[ijet][ipf] in lepCandKey[ilep]: 
                        if daughterKey in lepCandKey[ilep]:
                            if options.debug:
                                print 'Event {0:d}, removing lepton with pt/eta/phi = {1:6.2f},{2:6.2f},{3:6.2f} from AK4 jet with pt/eta/phi = {4:6.2f},{5:6.2f},{6:6.2f}'.format(ntotal,lepCand[ilep].Perp(),lepCand[ilep].Eta(),lepCand[ilep].Phi(),jetP4Raw.Perp(), jetP4Raw.Eta(), jetP4Raw.Phi())
                            if lepCand[ilep].E() > jetP4Raw.E() :
                                zeroedEnergy = True
                            cleanedLepton = lepCand[ilep]
                            jetP4Raw -= lepCand[ilep]
                            break

            if zeroedEnergy: 
                continue                

            # apply back jet energy corrections
            ak4JetCorrector.setJetEta( jetP4Raw.Eta() )
            ak4JetCorrector.setJetPt( jetP4Raw.Perp() )
            ak4JetCorrector.setJetE( jetP4Raw.E() )
            ak4JetCorrector.setJetA( ak4JetAreas[ijet] )
            ak4JetCorrector.setRho( rho )
            ak4JetCorrector.setNPV( NPV )

            newJEC = ak4JetCorrector.getCorrection()
            jetP4 = jetP4Raw*newJEC

            if cleanedLepton != 0: #Correct MET if JEC changed due to cleaning
                metv += jetP4Pre
                metv -= cleanedLepton
                metv -= jetP4
            else:                  #Correct MET due to adjusting JEC
                metv += jetP4Pre
                metv -= jetP4

            UncertJetAK4.setJetEta(jetP4Raw.Eta())
            UncertJetAK4.setJetPhi(jetP4Raw.Phi())
            UncertJetAK4.setJetPt(jetP4.Perp())
            unc = UncertJetAK4.getUncertainty(True)

            # JER
            jetP4_JERup = ROOT.TLorentzVector
            jetP4_JERdown = ROOT.TLorentzVector

            metv += jetP4 #Now correct MET for JER
            
            # Scale jet pt if there is a matched gen jet
            if options.isMC:
                if ak4MatchedGenJetPts[ijet] > 0:
                    genJetP4 = ROOT.TLorentzVector()
                    genJetP4.SetPtEtaPhiE(ak4MatchedGenJetPts[ijet],ak4MatchedGenJetEtas[ijet],ak4MatchedGenJetPhis[ijet],ak4MatchedGenJetEnergys[ijet])
                    if jetP4.DeltaR(genJetP4) < 0.2 and abs(jetP4.Perp() - genJetP4.Perp()) < (3 * ak4JERs[ijet] * jetP4.Perp()) : # Do matching requirement
                        jetP4 -= genJetP4
                        jetP4 *= ak4JERSFnoms[ijet]
                        jetP4 += genJetP4
                        jetP4_JERup = jetP4 - genJetP4
                        jetP4_JERup *= ak4JERSFups[ijet]
                        jetP4_JERup += genJetP4
                        jetP4_JERdown = jetP4 - genJetP4
                        jetP4_JERdown *= ak4JERSFdowns[ijet]
                        jetP4_JERdown += genJetP4
                    else:
                        jetP4 *= (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)))
                        jetP4_JERup   = jetP4 * (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFups[ijet]*ak4JERSFups[ijet]-1)))
                        jetP4_JERdown = jetP4 * (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFdowns[ijet]*ak4JERSFdowns[ijet]-1)))
                else:
                    jetP4 *= (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)))
                    jetP4_JERup   = jetP4 * (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFups[ijet]*ak4JERSFups[ijet]-1)))
                    jetP4_JERdown = jetP4 * (1.0 + smearfunc.Gaus(0.0,ak4JERs[ijet])*math.sqrt(max(0.0,ak4JERSFdowns[ijet]*ak4JERSFdowns[ijet]-1)))

            metv -= jetP4

            if jetP4.Perp() > 15. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt15.append(jetP4)
            if jetP4.Perp() > 20. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt20.append(jetP4)
            if jetP4.Perp() > 25. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt25.append(jetP4)
            if jetP4.Perp() > 30. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt30.append(jetP4)
            if jetP4.Perp() > 35. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt35.append(jetP4)
            if jetP4.Perp() > 40. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt40.append(jetP4)
            if jetP4.Perp() > 45. and abs(jetP4.Eta()) < 3.0:
                jetsFor2DPt45.append(jetP4)

            if jetP4.Perp() < MIN_JET_PT or abs(jetP4.Eta()) > MAX_JET_ETA:
                continue

            # calculate HT
            HT += jetP4.Perp()

            ak4jetPt.push_back(jetP4.Perp())
            ak4jetEta.push_back(jetP4.Eta())
            ak4jetPhi.push_back(jetP4.Phi())
            ak4jetMass.push_back(jetP4.M())
            ak4jetCSV.push_back(ak4JetCSVs[ijet])
            if (ak4JetVtxMasses[ijet] >= 0.0) :
                ak4jetVtxMass.push_back(ak4JetVtxMasses[ijet])
            else :
                ak4jetVtxMass.push_back(-1.0)
            ak4jetJECunc.push_back(unc)
            if options.isMC:
                ak4jetHadronFlavour.push_back(ak4jetHadronFlavours[ijet])
                ak4jetPtJERup.push_back(jetP4_JERup.Perp())
                ak4jetPtJERdown.push_back(jetP4_JERdown.Perp())
                ak4jetEtaJERup.push_back(jetP4_JERup.Eta())
                ak4jetEtaJERdown.push_back(jetP4_JERdown.Eta())
                ak4jetPhiJERup.push_back(jetP4_JERup.Phi())
                ak4jetPhiJERdown.push_back(jetP4_JERdown.Phi())
                ak4jetMassJERup.push_back(jetP4_JERup.M())
                ak4jetMassJERdown.push_back(jetP4_JERdown.M())

    # -------------------------------------------------------
    # If not storing full truth information, require AK4 jet
    # -------------------------------------------------------

    if len(ak4jetPt) == 0 :
        passReco = False
        if not options.fullTruth :
            continue
    else :
        nPassAK4jet += 1
        
    # -------------------------------------------------------------------------------------
    # now we can calculate electron / muon 2D cuts
    # -------------------------------------------------------------------------------------

    if len(elCand) > 0:
        for electron in elCand :

            if len(jetsFor2DPt15) >= 1:             
                eleJetPt15 = findClosestInList(electron, jetsFor2DPt15)
                elPtRelPt15.push_back(electron.Perp(eleJetPt15.Vect()))
                eldRPt15.push_back(electron.DeltaR(eleJetPt15))

            if len(jetsFor2DPt20) >= 1:             
                eleJetPt20 = findClosestInList(electron, jetsFor2DPt20)
                elPtRelPt20.push_back(electron.Perp(eleJetPt20.Vect()))
                eldRPt20.push_back(electron.DeltaR(eleJetPt20))
                    
            if len(jetsFor2DPt25) >= 1:             
                eleJetPt25 = findClosestInList(electron, jetsFor2DPt25)
                elPtRelPt25.push_back(electron.Perp(eleJetPt25.Vect()))
                eldRPt25.push_back(electron.DeltaR(eleJetPt25))

            if len(jetsFor2DPt30) >= 1:             
                eleJetPt30 = findClosestInList(electron, jetsFor2DPt30)
                elPtRelPt30.push_back(electron.Perp(eleJetPt30.Vect()))
                eldRPt30.push_back(electron.DeltaR(eleJetPt30))

            if len(jetsFor2DPt35) >= 1:             
                eleJetPt35 = findClosestInList(electron, jetsFor2DPt35)
                elPtRelPt35.push_back(electron.Perp(eleJetPt35.Vect()))
                eldRPt35.push_back(electron.DeltaR(eleJetPt35))

            if len(jetsFor2DPt40) >= 1:             
                eleJetPt40 = findClosestInList(electron, jetsFor2DPt40)
                elPtRelPt40.push_back(electron.Perp(eleJetPt40.Vect()))
                eldRPt40.push_back(electron.DeltaR(eleJetPt40))

            if len(jetsFor2DPt45) >= 1:             
                eleJetPt45 = findClosestInList(electron, jetsFor2DPt45)
                elPtRelPt45.push_back(electron.Perp(eleJetPt45.Vect()))
                eldRPt45.push_back(electron.DeltaR(eleJetPt45))

    if len(muCand) > 0:
        for muon in muCand :

            if len(jetsFor2DPt15) >= 1: 
                muJetPt15 = findClosestInList(muon, jetsFor2DPt15)
                muPtRelPt15.push_back(muon.Perp(muJetPt15.Vect()))
                mudRPt15.push_back(muon.DeltaR(muJetPt15))

            if len(jetsFor2DPt20) >= 1: 
                muJetPt20 = findClosestInList(muon, jetsFor2DPt20)
                muPtRelPt20.push_back(muon.Perp(muJetPt20.Vect()))
                mudRPt20.push_back(muon.DeltaR(muJetPt20))

            if len(jetsFor2DPt25) >= 1: 
                muJetPt25 = findClosestInList(muon, jetsFor2DPt25)
                muPtRelPt25.push_back(muon.Perp(muJetPt25.Vect()))
                mudRPt25.push_back(muon.DeltaR(muJetPt25))

            if len(jetsFor2DPt30) >= 1: 
                muJetPt30 = findClosestInList(muon, jetsFor2DPt30)
                muPtRelPt30.push_back(muon.Perp(muJetPt30.Vect()))
                mudRPt30.push_back(muon.DeltaR(muJetPt30))

            if len(jetsFor2DPt35) >= 1: 
                muJetPt35 = findClosestInList(muon, jetsFor2DPt35)
                muPtRelPt35.push_back(muon.Perp(muJetPt35.Vect()))
                mudRPt35.push_back(muon.DeltaR(muJetPt35))

            if len(jetsFor2DPt40) >= 1: 
                muJetPt40 = findClosestInList(muon, jetsFor2DPt40)
                muPtRelPt40.push_back(muon.Perp(muJetPt40.Vect()))
                mudRPt40.push_back(muon.DeltaR(muJetPt40))

            if len(jetsFor2DPt45) >= 1: 
                muJetPt45 = findClosestInList(muon, jetsFor2DPt45)
                muPtRelPt45.push_back(muon.Perp(muJetPt45.Vect()))
                mudRPt45.push_back(muon.DeltaR(muJetPt45))


    # -------------------------------------------------------------------------------------
    # read variables for AK8 jets
    # -------------------------------------------------------------------------------------

    event.getByLabel (ak8JetPtLabel, ak8JetPtHandle)
    if ak8JetPtHandle.isValid() :
        ak8JetPt = ak8JetPtHandle.product()    

        event.getByLabel (ak8JetPhiLabel, ak8JetPhiHandle)
        ak8JetPhi = ak8JetPhiHandle.product()
        event.getByLabel (ak8JetEtaLabel, ak8JetEtaHandle)
        ak8JetEta = ak8JetEtaHandle.product()
        event.getByLabel (ak8JetELabel, ak8JetEHandle)
        ak8JetE = ak8JetEHandle.product()
        event.getByLabel (ak8JetTrimMassLabel, ak8JetTrimMassHandle)
        ak8JetTrimMass = ak8JetTrimMassHandle.product()
        event.getByLabel (ak8JetPrunMassLabel, ak8JetPrunMassHandle)
        ak8JetPrunMass = ak8JetPrunMassHandle.product()
        event.getByLabel (ak8JetFiltMassLabel, ak8JetFiltMassHandle)
        ak8JetFiltMass = ak8JetFiltMassHandle.product()
        event.getByLabel (ak8JetTau1Label, ak8JetTau1Handle)
        ak8JetTau1 = ak8JetTau1Handle.product()
        event.getByLabel (ak8JetTau2Label, ak8JetTau2Handle)
        ak8JetTau2 = ak8JetTau2Handle.product()
        event.getByLabel (ak8JetTau3Label, ak8JetTau3Handle)
        ak8JetTau3 = ak8JetTau3Handle.product()
        event.getByLabel (ak8JetCSVLabel, ak8JetCSVHandle)
        ak8JetCSV = ak8JetCSVHandle.product()
        event.getByLabel (ak8JetSoftDropMassLabel, ak8JetSoftDropMassHandle)
        ak8JetSoftDropMass = ak8JetSoftDropMassHandle.product()
        event.getByLabel (ak8JetSoftDropSubjet0Label, ak8JetSoftDropSubjet0Handle)
        ak8JetSDSj0 = ak8JetSoftDropSubjet0Handle.product()
        event.getByLabel (ak8JetSoftDropSubjet1Label, ak8JetSoftDropSubjet1Handle)
        ak8JetSDSj1 = ak8JetSoftDropSubjet1Handle.product()

        # Get subjet info
        event.getByLabel(sjSoftDropPtLabel, sjSoftDropPtHandle)
        sjSoftDropPt   = sjSoftDropPtHandle.product()
        event.getByLabel(sjSoftDropEtaLabel, sjSoftDropEtaHandle)
        sjSoftDropEta  = sjSoftDropEtaHandle.product()
        event.getByLabel(sjSoftDropPhiLabel, sjSoftDropPhiHandle)
        sjSoftDropPhi  = sjSoftDropPhiHandle.product()
        event.getByLabel(sjSoftDropELabel, sjSoftDropEHandle)
        sjSoftDropE    = sjSoftDropEHandle.product()
        event.getByLabel(sjSoftDropCSVLabel, sjSoftDropCSVHandle)
        sjSoftDropCSV  = sjSoftDropCSVHandle.product()

        event.getByLabel (ak8JetAreaLabel, ak8JetAreaHandle)
        ak8JetAreas = ak8JetAreaHandle.product()
        event.getByLabel( ak8JetJECLabel, ak8JetJECHandle )
        ak8JetJECs = ak8JetJECHandle.product()

        event.getByLabel( ak8JetNeuHadEnergyFracLabel, ak8JetNeuHadEnergyFracHandle )
        event.getByLabel( ak8JetNeuEmEnergyFracLabel, ak8JetNeuEmEnergyFracHandle )
        event.getByLabel( ak8JetChHadEnergyFracLabel, ak8JetChHadEnergyFracHandle )
        event.getByLabel( ak8JetChEmEnergyFracLabel, ak8JetChEmEnergyFracHandle )
        event.getByLabel( ak8JetNeuMultiLabel, ak8JetNeuMultiHandle )
        event.getByLabel( ak8JetChMultiLabel, ak8JetChMultiHandle )

        ak8JetNeuHadEnergyFracs = ak8JetNeuHadEnergyFracHandle.product()
        ak8JetNeuEmEnergyFracs = ak8JetNeuEmEnergyFracHandle.product()
        ak8JetChHadEnergyFracs = ak8JetChHadEnergyFracHandle.product()
        ak8JetChEmEnergyFracs = ak8JetChEmEnergyFracHandle.product()
        ak8JetNeuMultis = ak8JetNeuMultiHandle.product()
        ak8JetChMultis = ak8JetChMultiHandle.product()

        event.getByLabel(ak8MatchedGenJetPtLabel, ak8MatchedGenJetPtHandle)
        ak8MatchedGenJetPts = ak8MatchedGenJetPtHandle.product()
        event.getByLabel(ak8MatchedGenJetEtaLabel, ak8MatchedGenJetEtaHandle)
        ak8MatchedGenJetEtas = ak8MatchedGenJetEtaHandle.product()
        event.getByLabel(ak8MatchedGenJetPhiLabel, ak8MatchedGenJetPhiHandle)
        ak8MatchedGenJetPhis = ak8MatchedGenJetPhiHandle.product()
        event.getByLabel(ak8MatchedGenJetEnergyLabel, ak8MatchedGenJetEnergyHandle)
        ak8MatchedGenJetEnergys = ak8MatchedGenJetEnergyHandle.product()
        event.getByLabel(ak8JERSFnomLabel, ak8JERSFnomHandle)
        ak8JERSFnoms = ak8JERSFnomHandle.product()
        event.getByLabel(ak8JERSFupLabel, ak8JERSFupHandle)
        ak8JERSFups = ak8JERSFupHandle.product()
        event.getByLabel(ak8JERSFdownLabel, ak8JERSFdownHandle)
        ak8JERSFdowns = ak8JERSFdownHandle.product()
        event.getByLabel(ak8JERLabel, ak8JERHandle)
        ak8JERs = ak8JERHandle.product()

        event.getByLabel(ak8JetKeysLabel, ak8JetKeysHandle)
        ak8JetKeys = ak8JetKeysHandle.product()

        # -------------------------------------------------------------------------------------
        # loop over AK8 jets
        # -------------------------------------------------------------------------------------

        # loop over jets
        for ijet in xrange( len(ak8JetPt) ) :
            
            # jet IDs must be calculated prior to JECs, so remove them 
            AK8jetP4Pre = ROOT.TLorentzVector()
            AK8jetP4Pre.SetPtEtaPhiE( ak8JetPt[ijet], ak8JetEta[ijet], ak8JetPhi[ijet], ak8JetE[ijet] )

            # get the raw jet energy
            AK8jetP4Raw = AK8jetP4Pre * ak8JetJECs[ijet]

            # calculate the neutral/charged hadron/em energy fractions
            AK8nhf = ak8JetNeuHadEnergyFracs[ijet]
            AK8nef = ak8JetNeuEmEnergyFracs[ijet]
            AK8chf = ak8JetChHadEnergyFracs[ijet]
            AK8cef = ak8JetChEmEnergyFracs[ijet]
            AK8nconstituents = ak8JetNeuMultis[ijet] + ak8JetChMultis[ijet]
            AK8nchmult = ak8JetChMultis[ijet] 

            # require loose jet ID
            if (AK8nhf >= 0.99 or AK8nef >= 0.99 or AK8chf <= 0.00 or AK8cef >= 0.99 or AK8nconstituents <= 1 or AK8nchmult <= 0) :
                continue

            # clean leptons from jets 
            zeroedEnergy = False

            AK8cleanedLepton = 0
            for ilep in xrange( len(lepCand) ) :
                if lepCand[ilep].DeltaR(AK8jetP4Raw) < 0.4 :
                    # check jet daughters close to the lepton
                    #pfcands = int(ak8JetNumDaughters[ijet])
                    #for ipf in range(0, pfcands) :
                    ak8daughters = ak8JetKeys[ijet]
                    for daughterKey in ak8daughters :
                        # if jet daughter matches lepton, remove lepton p4 from (raw) jet p4
                        #if ak8JetKeys[ijet][ipf] in lepCandKey[ilep]:
                        if daughterKey in lepCandKey[ilep]:
                            if options.debug:
                                print 'Event {0:d}, removing lepton with pt/eta/phi = {1:6.2f},{2:6.2f},{3:6.2f} from AK8 jet with pt/eta/phi = {4:6.2f},{5:6.2f},{6:6.2f}'.format(ntotal,lepCand[ilep].Perp(),lepCand[ilep].Eta(),lepCand[ilep].Phi(),AK8jetP4Raw.Perp(), AK8jetP4Raw.Eta(), AK8jetP4Raw.Phi())
                            if lepCand[ilep].E() > AK8jetP4Raw.E() :
                                zeroedEnergy = True
                            AK8cleanedLepton = lepCand[ilep]
                            AK8jetP4Raw -= lepCand[ilep]
                            break

            if zeroedEnergy: 
                continue

            # apply back jet energy corrections
            ak8JetCorrector.setJetEta( AK8jetP4Raw.Eta() )
            ak8JetCorrector.setJetPt( AK8jetP4Raw.Perp() )
            ak8JetCorrector.setJetE( AK8jetP4Raw.E() )
            ak8JetCorrector.setJetA( ak8JetAreas[ijet] )
            ak8JetCorrector.setRho( rho )
            ak8JetCorrector.setNPV( NPV )

            AK8newJEC = ak8JetCorrector.getCorrection()
            AK8jetP4 = AK8jetP4Raw*AK8newJEC

            if AK8cleanedLepton != 0: #Correct MET if JEC changed due to cleaning
                metv += AK8jetP4Pre
                metv -= AK8cleanedLepton
                metv -= AK8jetP4
            else:                     #Correct MET due to adjusting JEC
                metv += AK8jetP4Pre
                metv -= AK8jetP4

            UncertJetAK8.setJetEta(AK8jetP4Raw.Eta())
            UncertJetAK8.setJetPhi(AK8jetP4Raw.Phi())
            UncertJetAK8.setJetPt(AK8jetP4.Perp())
            unc = UncertJetAK8.getUncertainty(True)

            # JER
            AK8jetP4_JERup = ROOT.TLorentzVector
            AK8jetP4_JERdown = ROOT.TLorentzVector

            metv += AK8jetP4 #Now correct MET for JER
            
            # Scale jet pt if there is a matched gen jet
            if options.isMC:
                if ak8MatchedGenJetPts[ijet] > 0:
                    AK8genJetP4 = ROOT.TLorentzVector()
                    AK8genJetP4.SetPtEtaPhiE(ak8MatchedGenJetPts[ijet],ak8MatchedGenJetEtas[ijet],ak8MatchedGenJetPhis[ijet],ak8MatchedGenJetEnergys[ijet])
                    if AK8jetP4.DeltaR(AK8genJetP4) < 0.2 and abs(AK8jetP4.Perp() - AK8genJetP4.Perp()) < (3 * ak8JERs[ijet] * AK8jetP4.Perp()) : # Do matching requirement
                        AK8jetP4 -= AK8genJetP4
                        AK8jetP4 *= ak8JERSFnoms[ijet]
                        AK8jetP4 += AK8genJetP4
                        AK8jetP4_JERup = AK8jetP4 - AK8genJetP4
                        AK8jetP4_JERup *= ak8JERSFups[ijet]
                        AK8jetP4_JERup += AK8genJetP4
                        AK8jetP4_JERdown = AK8jetP4 - AK8genJetP4
                        AK8jetP4_JERdown *= ak8JERSFdowns[ijet]
                        AK8jetP4_JERdown += AK8genJetP4

                    else:
                        AK8jetP4 *= (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)))
                        AK8jetP4_JERup   = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFups[ijet]*ak8JERSFups[ijet]-1)))
                        AK8jetP4_JERdown = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFdowns[ijet]*ak8JERSFdowns[ijet]-1)))
                else:
                    AK8jetP4 *= (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)))
                    AK8jetP4_JERup   = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFups[ijet]*ak8JERSFups[ijet]-1)))
                    AK8jetP4_JERdown = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,ak8JERs[ijet])*math.sqrt(max(0.0,ak8JERSFdowns[ijet]*ak8JERSFdowns[ijet]-1)))

            metv -= AK8jetP4

            if (AK8jetP4.Perp() < TOP_PT_CUT or abs(AK8jetP4.Eta()) > MAX_JET_ETA):
                continue

            ak8jetPt.push_back(AK8jetP4.Perp())
            ak8jetEta.push_back(AK8jetP4.Eta())
            ak8jetPhi.push_back(AK8jetP4.Phi())
            ak8jetMass.push_back(AK8jetP4.M())
            if ak8JetPrunMass[ijet] >= 0.0 :
                ak8jetMassPruned.push_back(ak8JetPrunMass[ijet])
            else :
                ak8jetMassPruned.push_back(-10.0)
            if ak8JetFiltMass[ijet] >= 0.0 :
                ak8jetMassFiltered.push_back(ak8JetFiltMass[ijet])
            else :
                ak8jetMassFiltered.push_back(-10.0)
            if ak8JetTrimMass[ijet] >= 0.0 :
                ak8jetMassTrimmed.push_back(ak8JetTrimMass[ijet])
            else :
                ak8jetMassTrimmed.push_back(-10.0)
            ak8jetTau1.push_back(ak8JetTau1[ijet])
            ak8jetTau2.push_back(ak8JetTau2[ijet])
            ak8jetTau3.push_back(ak8JetTau3[ijet])
            ak8jetCSV.push_back(ak8JetCSV[ijet])
            if ak8JetSoftDropMass[ijet] >= 0.0 :
                ak8jetSDmass.push_back(ak8JetSoftDropMass[ijet])
            else :
                ak8jetSDmass.push_back(-10.0)
            if int(ak8JetSDSj0[ijet]) >= 0 and int(ak8JetSDSj0[ijet]) < len(sjSoftDropPt) :
                sj0P4 = ROOT.TLorentzVector()
                sj0P4.SetPtEtaPhiE(sjSoftDropPt[int(ak8JetSDSj0[ijet])],sjSoftDropEta[int(ak8JetSDSj0[ijet])],sjSoftDropPhi[int(ak8JetSDSj0[ijet])],sjSoftDropE[int(ak8JetSDSj0[ijet])])
                ak8jetSDsubjet0pt.push_back(sj0P4.Perp())
                ak8jetSDsubjet0eta.push_back(sj0P4.Eta())
                ak8jetSDsubjet0phi.push_back(sj0P4.Phi())
                ak8jetSDsubjet0mass.push_back(sj0P4.M())
                ak8jetSDsubjet0CSV.push_back(sjSoftDropCSV[int(ak8JetSDSj0[ijet])])
            else :
                ak8jetSDsubjet0pt.push_back(-10.0)
                ak8jetSDsubjet0eta.push_back(-10.0)
                ak8jetSDsubjet0phi.push_back(-10.0)
                ak8jetSDsubjet0mass.push_back(-10.0)
                ak8jetSDsubjet0CSV.push_back(-10.0)
            if int(ak8JetSDSj1[ijet]) >= 0 and int(ak8JetSDSj1[ijet]) < len(sjSoftDropPt) :
                sj1P4 = ROOT.TLorentzVector()
                sj1P4.SetPtEtaPhiE(sjSoftDropPt[int(ak8JetSDSj1[ijet])],sjSoftDropEta[int(ak8JetSDSj1[ijet])],sjSoftDropPhi[int(ak8JetSDSj1[ijet])],sjSoftDropE[int(ak8JetSDSj1[ijet])])
                ak8jetSDsubjet1pt.push_back(sj1P4.Perp())
                ak8jetSDsubjet1eta.push_back(sj1P4.Eta())
                ak8jetSDsubjet1phi.push_back(sj1P4.Phi())
                ak8jetSDsubjet1mass.push_back(sj1P4.M())
                ak8jetSDsubjet1CSV.push_back(sjSoftDropCSV[int(ak8JetSDSj1[ijet])])
            else :
                ak8jetSDsubjet1pt.push_back(-10.0)
                ak8jetSDsubjet1eta.push_back(-10.0)
                ak8jetSDsubjet1phi.push_back(-10.0)
                ak8jetSDsubjet1mass.push_back(-10.0)
                ak8jetSDsubjet1CSV.push_back(-10.0)
            ak8jetJECunc.push_back(unc)
            if options.isMC:
                ak8jetPtJERup.push_back(AK8jetP4_JERup.Perp())
                ak8jetPtJERdown.push_back(AK8jetP4_JERdown.Perp())
                ak8jetEtaJERup.push_back(AK8jetP4_JERup.Eta())
                ak8jetEtaJERdown.push_back(AK8jetP4_JERdown.Eta())
                ak8jetPhiJERup.push_back(AK8jetP4_JERup.Phi())
                ak8jetPhiJERdown.push_back(AK8jetP4_JERdown.Phi())
                ak8jetMassJERup.push_back(AK8jetP4_JERup.M())
                ak8jetMassJERdown.push_back(AK8jetP4_JERdown.M())

    if len(ak8jetPt) == 0 :
        passReco = False
        if not options.fullTruth :
            continue
    else :
        nPassAK8jet += 1
                

    metPt.push_back(metv.Perp())
    metPhi.push_back(metv.Phi())
    ht.push_back(HT)

    if passReco :
        nEventsPass += 1
        
    eventWeight_nom.push_back(weight_nom)
    if options.isMC:
        eventWeight_puUp.push_back(weight_puUp)
        eventWeight_puDown.push_back(weight_puDown)
        if options.semilep is not None:
            eventWeight_Q2Up.push_back(weight_Q2up)
            eventWeight_Q2Down.push_back(weight_Q2down)
            eventWeight_PDF.push_back(weight_PDF)
            eventWeight_alphaUp.push_back(weight_alphaUp)
            eventWeight_alphaDown.push_back(weight_alphaDown)

    if passReco :
        recoTree.Fill()

    if options.fullTruth and not passReco:
        trueTree.Fill()
    
# -------------------------------------------------------------------------------------
# END OF LOOPING OVER EVENTS!!!
# -------------------------------------------------------------------------------------

h_cutflow.Fill(1,ntotal)
h_cutflow.Fill(2,nPassNPV)
h_cutflow.Fill(3,nPassMetFilter)
h_cutflow.Fill(4,nPassRho)
h_cutflow.Fill(5,nPassLep)
h_cutflow.Fill(6,nPassAK4jet)
h_cutflow.Fill(7,nPassAK8jet)
h_cutflow.Fill(8,nEventsPass)

print 'Total Events:    ' + str(ntotal)
if not options.isMC :
    print 'Events in given IOV: ' + str(nInIOV)
print '-------------------------------'
print 'Parton level breakdown:'
print 'Muon channel: ' + str(nMuChParton)
print 'Electron channel: ' + str(nElChParton)
print 'Tau (only) channel: ' + str(nTauOnlyParton)
print 'Dilepton channel: ' + str(nDilepParton)
print 'Hadronic channel: ' + str(nHadParton)
print 'Hadronic top in wrong hemisphere: ' + str(nWrongHemisphere)
print '-------------------------------'
print 'Pass nPV:        ' + str(nPassNPV)
if options.semilep == 1:
    print 'Pass LHE:        ' + str(nPassLHE)
print 'Pass MET filter: ' + str(nPassMetFilter)
print 'Pass rho:        ' + str(nPassRho)
print '-------------------------------'
print 'Pass parton:     ' + str(nPassParton)
print '-------------------------------'
print 'Pass particle:   ' + str(nPassParticle)
print '-------------------------------'
print 'Pass lepton:     ' + str(nPassLep)
print 'Pass AK4 jet:    ' + str(nPassAK4jet)
print 'Pass AK8 jet:    ' + str(nPassAK8jet)
print 'Pass reco:       ' + str(nEventsPass)
print '-------------------------------'
#if nEleRaw != 0 :
#    print 'Fraction of all medium electrons not passing twiki cuts: ' + str(nMediumNotManualEl / nEleRaw)
#    print 'Fraction of all electrons passing twiki cuts not medium: ' + str(nManualNotMediumEl / nEleRaw)
#    print 'Fraction of all tight electrons not passing twiki cuts:  ' + str(nTightNotManualEl / nEleRaw)
#    print 'Fraction of all electrons passing twiki cuts not tight:  ' + str(nManualNotTightEl / nEleRaw)
#    print '-------------------------------'
#    print 'Fraction of all medium no iso electrons not passing twiki cuts: ' + str(nMediumNoIsoNotManualEl / nEleRaw)
#    print 'Fraction of all electrons passing twiki cuts not medium no iso: ' + str(nManualNotMediumNoIsoEl / nEleRaw)
#    print 'Fraction of all tight no iso electrons not passing twiki cuts:  ' + str(nTightNoIsoNotManualEl / nEleRaw)
#    print 'Fraction of all electrons passing twiki cuts not tight no iso:  ' + str(nManualNotTightNoIsoEl / nEleRaw)
#else :
#    print 'No electrons!'
#print '-------------------------------'
if nMuRaw != 0:
    print 'Fraction of all tight muons not passing twiki cuts:      ' + str(nTightNotManualMu / nMuRaw)
    print 'Fraction of all muons passing twiki cuts not tight:      ' + str(nManualNotTightMu / nMuRaw)
else:
    print 'No muons!'
print '-------------------------------'
    
f.cd()
f.Write()
f.Close()

