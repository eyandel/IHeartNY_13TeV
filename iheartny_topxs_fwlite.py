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

parser.add_option('--debug', metavar='M', action='store_true',
                  default=False,
                  dest='debug',
                  help='Print out debug statements')

parser.add_option('--fullTruth', metavar='M', action='store_true',
                  default=False,
                  dest='fullTruth',
                  help='Save truth info for all events')

parser.add_option('--puFile', metavar='F', type='string', action='store',
                  default=None,
                  dest='puFile',
                  help='Pileup reweighting file')


(options, args) = parser.parse_args()
argv = []

import ROOT
ROOT.gROOT.Macro("rootlogon.C")

from array import *

import sys
from DataFormats.FWLite import Events, Handle


# -------------------------------------------------------------------------------------
# jet energy corrections
# -------------------------------------------------------------------------------------

ROOT.gSystem.Load('libCondFormatsJetMETObjects')

jetname = "chs"
if options.usePuppi:
    jetname = "puppi"

if options.isMC : 
    L3JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L3Absolute_AK4PF'+jetname+'.txt')
    L3JetParAK4 = ROOT.JetCorrectorParameters(L3JecStrAK4)
    L2JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L2Relative_AK4PF'+jetname+'.txt')
    L2JetParAK4 = ROOT.JetCorrectorParameters(L2JecStrAK4)
    L1JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L1FastJet_AK4PF'+jetname+'.txt')
    L1JetParAK4 = ROOT.JetCorrectorParameters(L1JecStrAK4)
    UncJecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_MC_Uncertainty_AK4PF'+jetname+'.txt')
    UncertJetAK4 = ROOT.JetCorrectionUncertainty(UncJecStrAK4)
    L3JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L3Absolute_AK8PF'+jetname+'.txt')
    L3JetParAK8 = ROOT.JetCorrectorParameters(L3JecStrAK8)
    L2JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L2Relative_AK8PF'+jetname+'.txt')
    L2JetParAK8 = ROOT.JetCorrectorParameters(L2JecStrAK8)
    L1JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_MC_L1FastJet_AK8PF'+jetname+'.txt')
    L1JetParAK8 = ROOT.JetCorrectorParameters(L1JecStrAK8)
    UncJecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_MC_Uncertainty_AK8PF'+jetname+'.txt')
    UncertJetAK8 = ROOT.JetCorrectionUncertainty(UncJecStrAK8)
    
else :
    L3JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L3Absolute_AK4PF'+jetname+'.txt')
    L3JetParAK4 = ROOT.JetCorrectorParameters(L3JecStrAK4)
    L2JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L2Relative_AK4PF'+jetname+'.txt')
    L2JetParAK4 = ROOT.JetCorrectorParameters(L2JecStrAK4)
    L1JecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L1FastJet_AK4PF'+jetname+'.txt')
    L1JetParAK4 = ROOT.JetCorrectorParameters(L1JecStrAK4)
    ResJecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L2L3Residual_AK4PF'+jetname+'.txt')
    ResJetParAK4 = ROOT.JetCorrectorParameters(ResJecStrAK4)
    UncJecStrAK4 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_Uncertainty_AK4PF'+jetname+'.txt')
    UncertJetAK4 = ROOT.JetCorrectionUncertainty(UncJecStrAK4)
    L3JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L3Absolute_AK8PF'+jetname+'.txt')
    L3JetParAK8 = ROOT.JetCorrectorParameters(L3JecStrAK8)
    L2JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L2Relative_AK8PF'+jetname+'.txt')
    L2JetParAK8 = ROOT.JetCorrectorParameters(L2JecStrAK8)
    L1JecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L1FastJet_AK8PF'+jetname+'.txt')
    L1JetParAK8 = ROOT.JetCorrectorParameters(L1JecStrAK8)
    ResJecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_L2L3Residual_AK8PF'+jetname+'.txt')
    ResJetParAK8 = ROOT.JetCorrectorParameters(ResJecStrAK8)
    UncJecStrAK8 = ROOT.std.string('JECs/Spring16_25nsV3_DATA_Uncertainty_AK8PF'+jetname+'.txt')
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

myTree = ROOT.TTree("myTree", "myTree")

# parton level
if options.isMC and options.semilep is not None:
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

if options.isMC :
    muTrigPass             = ROOT.vector('int')()
    elTrigPass             = ROOT.vector('int')()
    if options.semilep == 1:
        truthChannel       = ROOT.vector('int')()

if options.isMC and options.semilep is not None:
    myTree.Branch('genTopPt'               , genTopPt               )
    myTree.Branch('genTopEta'              , genTopEta              )
    myTree.Branch('genTopPhi'              , genTopPhi              )
    myTree.Branch('genMuPt'                , genMuPt                )
    myTree.Branch('genMuEta'               , genMuEta               )
    myTree.Branch('genMuPhi'               , genMuPhi               )
    myTree.Branch('genElPt'                , genElPt                )
    myTree.Branch('genElEta'               , genElEta               )
    myTree.Branch('genElPhi'               , genElPhi               )
    myTree.Branch('genTTbarMass'           , genTTbarMass           )
    
    #myTree.Branch('genAK4jetPt'            , genAK4jetPt            )
    #myTree.Branch('genAK4jetEta'           , genAK4jetEta           )
    #myTree.Branch('genAK4jetPhi'           , genAK4jetPhi           )
    #myTree.Branch('genAK4jetMass'          , genAK4jetMass          )
    myTree.Branch('genAK8jetPt'            , genAK8jetPt            )
    myTree.Branch('genAK8jetEta'           , genAK8jetEta           )
    myTree.Branch('genAK8jetPhi'           , genAK8jetPhi           )
    myTree.Branch('genAK8jetMass'          , genAK8jetMass          )
    myTree.Branch('partMuPt'               , partMuPt               )
    myTree.Branch('partMuEta'              , partMuEta              )
    myTree.Branch('partMuPhi'              , partMuPhi              )
    myTree.Branch('partElPt'               , partElPt               )
    myTree.Branch('partElEta'              , partElEta              )
    myTree.Branch('partElPhi'              , partElPhi              )

myTree.Branch('metPt'                  , metPt                  )
myTree.Branch('metPhi'                 , metPhi                 )
myTree.Branch('ht'                     , ht                     )
myTree.Branch('muPt'                   , muPt                   )
myTree.Branch('muEta'                  , muEta                  )
myTree.Branch('muPhi'                  , muPhi                  )
myTree.Branch('muMiniIso'              , muMiniIso              )
myTree.Branch('muPtRelPt15'            , muPtRelPt15            )
myTree.Branch('muPtRelPt20'            , muPtRelPt20            )
myTree.Branch('muPtRelPt25'            , muPtRelPt25            )
myTree.Branch('muPtRelPt30'            , muPtRelPt30            )
myTree.Branch('muPtRelPt35'            , muPtRelPt35            )
myTree.Branch('muPtRelPt40'            , muPtRelPt40            )
myTree.Branch('muPtRelPt45'            , muPtRelPt45            )
myTree.Branch('mudRPt15'               , mudRPt15               )
myTree.Branch('mudRPt20'               , mudRPt20               )
myTree.Branch('mudRPt25'               , mudRPt25               )
myTree.Branch('mudRPt30'               , mudRPt30               )
myTree.Branch('mudRPt35'               , mudRPt35               )
myTree.Branch('mudRPt40'               , mudRPt40               )
myTree.Branch('mudRPt45'               , mudRPt45               )
myTree.Branch('muTight'                , muTight                )
myTree.Branch('elPt'                   , elPt                   )
myTree.Branch('elEta'                  , elEta                  )
myTree.Branch('elPhi'                  , elPhi                  )
myTree.Branch('elMiniIso'              , elMiniIso              )
myTree.Branch('elPtRelPt15'            , elPtRelPt15            )
myTree.Branch('elPtRelPt20'            , elPtRelPt20            )
myTree.Branch('elPtRelPt25'            , elPtRelPt25            )
myTree.Branch('elPtRelPt30'            , elPtRelPt30            )
myTree.Branch('elPtRelPt35'            , elPtRelPt35            )
myTree.Branch('elPtRelPt40'            , elPtRelPt40            )
myTree.Branch('elPtRelPt45'            , elPtRelPt45            )
myTree.Branch('eldRPt15'               , eldRPt15               )
myTree.Branch('eldRPt20'               , eldRPt20               )
myTree.Branch('eldRPt25'               , eldRPt25               )
myTree.Branch('eldRPt30'               , eldRPt30               )
myTree.Branch('eldRPt35'               , eldRPt35               )
myTree.Branch('eldRPt40'               , eldRPt40               )
myTree.Branch('eldRPt45'               , eldRPt45               )
myTree.Branch('elTight'                , elTight                )
myTree.Branch('ak4jetPt'               , ak4jetPt               )
myTree.Branch('ak4jetEta'              , ak4jetEta              )
myTree.Branch('ak4jetPhi'              , ak4jetPhi              )
myTree.Branch('ak4jetMass'             , ak4jetMass             )
myTree.Branch('ak4jetCSV'              , ak4jetCSV              )
myTree.Branch('ak4jetVtxMass'          , ak4jetVtxMass          )
myTree.Branch('ak4jetJECunc'           , ak4jetJECunc           )
if options.isMC:
    myTree.Branch('ak4jetHadronFlavour'    , ak4jetHadronFlavour    )
    myTree.Branch('ak4jetPtJERup'          , ak4jetPtJERup          )
    myTree.Branch('ak4jetPtJERdown'        , ak4jetPtJERdown        )
    myTree.Branch('ak4jetEtaJERup'         , ak4jetEtaJERup         )
    myTree.Branch('ak4jetEtaJERdown'       , ak4jetEtaJERdown       )
    myTree.Branch('ak4jetPhiJERup'         , ak4jetPhiJERup         )
    myTree.Branch('ak4jetPhiJERdown'       , ak4jetPhiJERdown       )
    myTree.Branch('ak4jetMassJERup'        , ak4jetMassJERup        )
    myTree.Branch('ak4jetMassJERdown'      , ak4jetMassJERdown      )
myTree.Branch('ak8jetPt'               , ak8jetPt               )
myTree.Branch('ak8jetEta'              , ak8jetEta              )
myTree.Branch('ak8jetPhi'              , ak8jetPhi              )
myTree.Branch('ak8jetMass'             , ak8jetMass             )
myTree.Branch('ak8jetMassPruned'       , ak8jetMassPruned       )
myTree.Branch('ak8jetMassFiltered'     , ak8jetMassFiltered     )
myTree.Branch('ak8jetMassTrimmed'      , ak8jetMassTrimmed      )
myTree.Branch('ak8jetTau1'             , ak8jetTau1             )
myTree.Branch('ak8jetTau2'             , ak8jetTau2             )
myTree.Branch('ak8jetTau3'             , ak8jetTau3             )
myTree.Branch('ak8jetCSV'              , ak8jetCSV              )
myTree.Branch('ak8jetSDmass'           , ak8jetSDmass           )
myTree.Branch('ak8jetSDsubjet0pt'      , ak8jetSDsubjet0pt      )
myTree.Branch('ak8jetSDsubjet0eta'     , ak8jetSDsubjet0eta     )
myTree.Branch('ak8jetSDsubjet0phi'     , ak8jetSDsubjet0phi     )
myTree.Branch('ak8jetSDsubjet0mass'    , ak8jetSDsubjet0mass    )
myTree.Branch('ak8jetSDsubjet0CSV'     , ak8jetSDsubjet0CSV     )
myTree.Branch('ak8jetSDsubjet1pt'      , ak8jetSDsubjet1pt      ) 
myTree.Branch('ak8jetSDsubjet1eta'     , ak8jetSDsubjet1eta     )
myTree.Branch('ak8jetSDsubjet1phi'     , ak8jetSDsubjet1phi     )
myTree.Branch('ak8jetSDsubjet1mass'    , ak8jetSDsubjet1mass    )
myTree.Branch('ak8jetSDsubjet1CSV'     , ak8jetSDsubjet1CSV     )
myTree.Branch('ak8jetJECunc'           , ak8jetJECunc           )
if options.isMC:
    myTree.Branch('ak8jetPtJERup'          , ak8jetPtJERup          )
    myTree.Branch('ak8jetPtJERdown'        , ak8jetPtJERdown        )
    myTree.Branch('ak8jetEtaJERup'         , ak8jetEtaJERup         )
    myTree.Branch('ak8jetEtaJERdown'       , ak8jetEtaJERdown       )
    myTree.Branch('ak8jetPhiJERup'         , ak8jetPhiJERup         )
    myTree.Branch('ak8jetPhiJERdown'       , ak8jetPhiJERdown       )
    myTree.Branch('ak8jetMassJERup'        , ak8jetMassJERup        )
    myTree.Branch('ak8jetMassJERdown'      , ak8jetMassJERdown      )

myTree.Branch('eventWeight_nom'           , eventWeight_nom           )

if options.isMC:
    myTree.Branch('eventWeight_puUp'          , eventWeight_puUp          )
    myTree.Branch('eventWeight_puDown'        , eventWeight_puDown        )
    myTree.Branch('muTrigPass'             , muTrigPass             )
    myTree.Branch('elTrigPass'             , elTrigPass             )
    if options.semilep == 1:
        myTree.Branch('truthChannel'            , truthChannel           )

# -------------------------------------------------------------------------------------
# define all variables to be read from input files
# -------------------------------------------------------------------------------------

events = Events (files)

jetname = "CHS"
if options.usePuppi:
    jetname = "Puppi"

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
HBHEfilterHandle    = Handle("bool")
HBHEfilterLabel     = ("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult")

# genParticles
if options.isMC and options.semilep is not None:
    genParticlesPtHandle     = Handle("std::vector<float>")
    genParticlesPtLabel      = ("genPart", "genPartPt")
    genParticlesEtaHandle    = Handle("std::vector<float>")
    genParticlesEtaLabel     = ("genPart", "genPartEta")
    genParticlesPhiHandle    = Handle("std::vector<float>")
    genParticlesPhiLabel     = ("genPart", "genPartPhi")
    genParticlesMassHandle   = Handle("std::vector<float>")
    genParticlesMassLabel    = ("genPart", "genPartMass")
    genParticlesPdgIdHandle  = Handle("std::vector<float>")
    genParticlesPdgIdLabel   = ("genPart", "genPartID")
    genParticlesStatusHandle = Handle("std::vector<float>")
    genParticlesStatusLabel  = ("genPart", "genPartStatus")
    genParticlesMom0StatusHandle = Handle("std::vector<float>")
    genParticlesMom0StatusLabel  = ("genPart", "genPartMom0Status")
    genParticlesMom0IDHandle     = Handle("std::vector<float>")
    genParticlesMom0IDLabel      = ("genPart", "genPartMom0ID")
    genParticlesDau0IDHandle     = Handle("std::vector<float>")
    genParticlesDau0IDLabel      = ("genPart", "genPartDau0ID")
    genParticlesDau1IDHandle     = Handle("std::vector<float>")
    genParticlesDau1IDLabel      = ("genPart", "genPartDau1ID")
    
    ## only have AK8 gen jets currently (in addition to gen jets matched to AK4/AK8 reco jets)
    ak8GenJetPtHandle   = Handle("std::vector<float>")
    ak8GenJetPtLabel    = ("genJetsAK8", "genJetsAK8Pt")
    ak8GenJetEtaHandle  = Handle("std::vector<float>")
    ak8GenJetEtaLabel   = ("genJetsAK8", "genJetsAK8Eta")
    ak8GenJetPhiHandle  = Handle("std::vector<float>")
    ak8GenJetPhiLabel   = ("genJetsAK8", "genJetsAK8Phi")
    ak8GenJetMassHandle = Handle("std::vector<float>")
    ak8GenJetMassLabel  = ("genJetsAK8", "genJetsAK8Mass")
    
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
electronTightHandle       = Handle("std::vector<float>")
electronTightLabel        = ("electrons" , "elvidTight" )
electronMediumHandle      = Handle("std::vector<float>")
electronMediumLabel       = ("electrons" , "elvidMedium" )
electronLooseHandle       = Handle("std::vector<float>")
electronLooseLabel        = ("electrons" , "elvidLoose" )
electronfull5x5sieeHandle = Handle("std::vector<float>")
electronfull5x5sieeLabel  = ( "electrons" , "elfull5x5siee")
electrondEtaInHandle      = Handle("std::vector<float>")
electrondEtaInLabel       = ( "electrons" , "eldEtaIn" )
electrondPhiInHandle      = Handle("std::vector<float>")
electrondPhiInLabel       = ( "electrons" , "eldPhiIn" )
electronHoEHandle         = Handle("std::vector<float>")
electronHoELabel          = ( "electrons" , "elHoE" )
electronooEmooPHandle     = Handle("std::vector<float>")
electronooEmooPLabel      = ( "electrons" , "elooEmooP")
electronD0Handle          = Handle("std::vector<float>")
electronD0Label           = ( "electrons" , "elD0" )
electronDzHandle          = Handle("std::vector<float>")
electronDzLabel           = ( "electrons" , "elDz" )
electronMiniIsoHandle     = Handle("std::vector<float>")
electronMiniIsoLabel      = ( "electrons" , "elMiniIso" )

elKeyHandle = Handle("std::vector<std::vector<int> >")
elKeyLabel = ( "electronKeys" )

# AK4 jet collection
ak4JetPtHandle   = Handle( "std::vector<float>" )
ak4JetPtLabel    = ("jetsAK4"+jetname, "jetAK4"+jetname+"Pt")
ak4JetEtaHandle  = Handle( "std::vector<float>" )
ak4JetEtaLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"Eta")
ak4JetPhiHandle  = Handle( "std::vector<float>" )
ak4JetPhiLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"Phi")
ak4JetMassHandle = Handle( "std::vector<float>" )
ak4JetMassLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"Mass")
ak4JetCSVHandle  = Handle( "std::vector<float>" )
ak4JetCSVLabel   = ("jetsAK4"+jetname, "jetAK4"+jetname+"CSVv2")
ak4JetVtxMassHandle = Handle( "std::vector<float>" )
ak4JetVtxMassLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"SV0mass")    
ak4JetAreaHandle = Handle( "std::vector<float>" )
ak4JetAreaLabel  = ("jetsAK4"+jetname, "jetAK4"+jetname+"jetArea")

ak4JetNeuHadEnergyHandle = Handle("std::vector<float>")
ak4JetNeuHadEnergyLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"neutralHadronEnergy")
ak4JetNeuEmEnergyHandle = Handle("std::vector<float>")
ak4JetNeuEmEnergyLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"neutralEmEnergy")
ak4JetChHadEnergyHandle = Handle("std::vector<float>")
ak4JetChHadEnergyLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"chargedHadronEnergy")
ak4JetChEmEnergyHandle = Handle("std::vector<float>")
ak4JetChEmEnergyLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"chargedEmEnergy")
ak4JetNumDaughterHandle = Handle("std::vector<float>")
ak4JetNumDaughterLabel = ("jetsAK4"+jetname , "jetAK4"+jetname+"numberOfDaughters")
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
ak8JetYHandle    = Handle( "std::vector<float>" )
ak8JetYLabel     = ("jetsAK8"+jetname, "jetAK8"+jetname+"Y" )
ak8JetMassHandle = Handle( "std::vector<float>" )
ak8JetMassLabel  = ("jetsAK8"+jetname, "jetAK8"+jetname+"Mass")
ak8JetTrimMassHandle = Handle("std::vector<float>")
ak8JetTrimMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"trimmedMass" )
ak8JetPrunMassHandle = Handle("std::vector<float>")
ak8JetPrunMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"prunedMass" )
ak8JetFiltMassHandle = Handle("std::vector<float>")
ak8JetFiltMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"filteredMass" )
ak8JetTau1Handle = Handle("std::vector<float>")
ak8JetTau1Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau1" )
ak8JetTau2Handle = Handle("std::vector<float>")
ak8JetTau2Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau2" )
ak8JetTau3Handle = Handle("std::vector<float>")
ak8JetTau3Label = ("jetsAK8"+jetname, "jetAK8"+jetname+"tau3" )
ak8JetCSVHandle = Handle("std::vector<float>")               
ak8JetCSVLabel = ( "jetsAK8"+jetname , "jetAK8"+jetname+"CSVv2" )

ak8JetSoftDropMassHandle = Handle("std::vector<float>")
ak8JetSoftDropMassLabel = ("jetsAK8"+jetname, "jetAK8"+jetname+"softDropMass" )
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
sjSoftDropMassHandle            = Handle( "std::vector<float>")
sjSoftDropMassLabel             = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"Mass")
sjSoftDropYHandle               = Handle( "std::vector<float>")
sjSoftDropYLabel                = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"Y")
sjSoftDropCSVHandle           = Handle( "std::vector<float>")
sjSoftDropCSVLabel            = ("subjetsAK8"+jetname, "subjetAK8"+jetname+"CSVv2")

# JETID variables
ak8JetNeuHadEnergyHandle = Handle("std::vector<float>")
ak8JetNeuHadEnergyLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"neutralHadronEnergy")
ak8JetNeuEmEnergyHandle = Handle("std::vector<float>")
ak8JetNeuEmEnergyLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"neutralEmEnergy")
ak8JetChHadEnergyHandle = Handle("std::vector<float>")
ak8JetChHadEnergyLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"chargedHadronEnergy")
ak8JetChEmEnergyHandle = Handle("std::vector<float>")
ak8JetChEmEnergyLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"chargedEmEnergy")
ak8JetNumDaughterHandle = Handle("std::vector<float>")
ak8JetNumDaughterLabel = ("jetsAK8"+jetname , "jetAK8"+jetname+"numberOfDaughters")
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

# -------------------------------------------------------------------------------------
# Histograms to save
# -------------------------------------------------------------------------------------

h_NtrueIntPU     = ROOT.TH1D("h_NtrueIntPU"    , "", 50,0,50 )
h_NPV_noweight   = ROOT.TH1F("h_NPV_noweight"  , "", 50,0,50 )
h_NPV            = ROOT.TH1F("h_NPV"           , "", 50,0,50 )
h_muPrescale     = ROOT.TH1F("h_muPrescale"    , "", 50,0,50 )
h_elPrescale     = ROOT.TH1F("h_elPrescale"    , "", 50,0,50 )
h_JERSF_AK4      = ROOT.TH1F("h_JERSF_AK4"       , "", 122,0.99,1.60)
h_JERSF_AK8      = ROOT.TH1F("h_JERSF_AK8"       , "", 122,0.99,1.60)

# -------------------------------------------------------------------------------------
# Get pileup weights
# -------------------------------------------------------------------------------------

if options.isMC and options.puFile is not None:
    fPUweight      = ROOT.TFile(options.puFile)
    hPUweight_nom  = fPUweight.Get("PUweight_true")
    hPUweight_up   = fPUweight.Get("PUweight_up")
    hPUweight_down = fPUweight.Get("PUweight_down")

# -------------------------------------------------------------------------------------
# reset various counters
# -------------------------------------------------------------------------------------

ntotal = 0       # total number of events

# Event quality
nPassNPV = 0
nPassMetFilter = 0
nPassRho = 0
nPassSemiLep = 0

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

# -------------------------------------------------------------------------------------
# start looping over events
# -------------------------------------------------------------------------------------

print "Start looping over events!"

for event in events :

    if options.isMC and options.semilep is not None:
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
        muTrigPass.clear()
        elTrigPass.clear()
        if options.semilep == 1:
            truthChannel.clear()
    
    weight_nom = 1.0 #event weight
    weight_puUp = 1.0
    weight_puDown = 1.0

    mttbarGen = -1.0
    
    if ntotal % 1000 == 0 :
      print  '--------- Processing Event ' + str(ntotal)
    ntotal += 1

    if options.debug :
        if ntotal > 1000:
            continue
        print 'Event ' + str(ntotal)

    filledEvent = False

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
        continue
    nPassNPV += 1

    # -------------------------------
    # Require met filter
    # -------------------------------
    
    metFilt = False
    hbheFilt = False

    gotName = event.getByLabel( metFilterNameLabel, metFilterNameHandle )
    gotBits = event.getByLabel( metFilterBitsLabel, metFilterBitsHandle )
    gotHBHE = event.getByLabel( HBHEfilterLabel, HBHEfilterHandle )

    if gotName == False or gotBits == False  :
        continue

    filterNameStrings = metFilterNameHandle.product()
    filterBits = metFilterBitsHandle.product()

    hbheFilt = HBHEfilterHandle.product()[0]

    for itrig in xrange(0, len(filterNameStrings) ) :
        if any(s in filterNameStrings[itrig] for s in ("globalTightHalo2016Filter","goodVertices","eeBadScFilter","EcalDeadCellTriggerPrimitiveFilter")) and filterBits[itrig] == 1:
            if options.debug :
                print 'MET filter: ' + filterNameStrings[itrig]
            metFilt = True

        if "HBHE" in filterNameStrings[itrig] and filterBits[itrig] == 1: #Basically OR of our HBHE filter and the one stored in METUserData
            if options.debug :
                print 'MET filter: ' + filterNameStrings[itrig]
            hbheFilt = True

    if hbheFilt == False or metFilt == False :
        continue
    nPassMetFilter += 1
                
    # -------------------------------
    # Require a trigger if MC
    # -------------------------------

    passMuTrig = True
    passElTrig = True

    # 80X does not have trigger stored! Turn this off for now
    #
    #if options.isMC :
    #    passMuTrig = False
    #    passElTrig = False
    #    prescale = 1.0
    #        
    #    event.getByLabel( trigNameLabel, trigNameHandle )
    #    event.getByLabel( trigBitsLabel, trigBitsHandle )
    #    event.getByLabel( trigPrescalesLabel, trigPrescalesHandle )
    #
    #    triggerNames = trigNameHandle.product()
    #    triggerBits = trigBitsHandle.product()
    #    triggerPrescales = trigPrescalesHandle.product()
    #
    #    trigToRun = None
    #    for itrig in xrange(0, len(triggerBits) ) :
    #        if triggerBits[itrig] != 1 :
    #            continue
    #        trigName = triggerNames[itrig]
    #        if "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50" in trigName :
    #            passElTrig = True
    #            prescale = prescale * triggerPrescales[itrig]
    #            h_elPrescale.Fill(triggerPrescales[itrig])
    #        if "HLT_Mu45_eta2p1" in trigName :
    #            passMuTrig = True
    #            prescale = prescale * triggerPrescales[itrig]
    #            h_muPrescale.Fill(triggerPrescales[itrig])
    #        
    #    weight_nom = weight_nom * prescale #Currently prescale has both mu and el prescales if both triggers fired
    #    weight_puUp = weight_puUp * prescale 
    #    weight_puDown = weight_puDown * prescale 

    # -------------------------------------------------------------------------------------
    # read event rho value
    # -------------------------------------------------------------------------------------

    event.getByLabel( rhoLabel, rhoHandle )
    if len(rhoHandle.product()) == 0 :
        print "Event has no rho values."
        continue
    rho = rhoHandle.product()[0]
    nPassRho += 1
    
    # -------------------------------------------------------------------------------------
    # Get parton-level info
    # -------------------------------------------------------------------------------------

    if options.isMC and options.semilep is not None :

        topQuarks = []
        genMuons = []
        genElectrons = []
        hadTop = None
        lepTop = None
        isSemiLeptonicGen = True
        isMuon = False
        isElectron = False
        
        event.getByLabel( genParticlesPtLabel, genParticlesPtHandle )
        event.getByLabel( genParticlesEtaLabel, genParticlesEtaHandle )
        event.getByLabel( genParticlesPhiLabel, genParticlesPhiHandle )
        event.getByLabel( genParticlesMassLabel, genParticlesMassHandle )
        event.getByLabel( genParticlesPdgIdLabel, genParticlesPdgIdHandle )
        event.getByLabel( genParticlesStatusLabel, genParticlesStatusHandle )
        event.getByLabel( genParticlesMom0IDLabel, genParticlesMom0IDHandle)
        event.getByLabel( genParticlesDau0IDLabel, genParticlesDau0IDHandle)
        event.getByLabel( genParticlesDau1IDLabel, genParticlesDau1IDHandle)
        
        if genParticlesPtHandle.isValid() == False :
            continue
        
        genParticlesPt  = genParticlesPtHandle.product()
        genParticlesEta = genParticlesEtaHandle.product()
        genParticlesPhi = genParticlesPhiHandle.product()
        genParticlesMass   = genParticlesMassHandle.product()
        genParticlesPdgId  = genParticlesPdgIdHandle.product()
        genParticlesStatus = genParticlesStatusHandle.product()
        genParticlesMom0ID = genParticlesMom0IDHandle.product()
        genParticlesDau0ID = genParticlesDau0IDHandle.product()
        genParticlesDau1ID = genParticlesDau1IDHandle.product()
        
        p4Top = ROOT.TLorentzVector()
        p4Antitop = ROOT.TLorentzVector()
        topDecay = 0        # 0 = hadronic, 1 = leptonic
        antitopDecay = 0    # 0 = hadronic, 1 = leptonic
        tauToEl = False
        tauToMu = False
        
        # loop over gen particles
        for igen in xrange( len(genParticlesPt) ) :
            
            # Find tops -- |pdgID| = 6, status 22
            if genParticlesPdgId[igen] == 6 and genParticlesStatus[igen] == 22 :
                p4Top.SetPtEtaPhiM( genParticlesPt[igen], genParticlesEta[igen], genParticlesPhi[igen], genParticlesMass[igen] )                    
            if genParticlesPdgId[igen] == -6 and genParticlesStatus[igen] == 22:
                p4Antitop.SetPtEtaPhiM( genParticlesPt[igen], genParticlesEta[igen], genParticlesPhi[igen], genParticlesMass[igen] )
            
            # If there is an antilepton (e+, mu+, tau+) in the W+ decay then the top is leptonic
            if ( genParticlesPdgId[igen] == -11 or genParticlesPdgId[igen] == -13 or genParticlesPdgId[igen] == -15) :
                if genParticlesMom0ID[igen] == 24 :
                    topDecay = 1
            # If there is an lepton (e-, mu-, tau-) in the W- decay then the antitop is leptonic
            if ( genParticlesPdgId[igen] == 11 or genParticlesPdgId[igen] == 13 or genParticlesPdgId[igen] == 15) :
                if genParticlesMom0ID[igen] == -24 :
                    antitopDecay = 1
            
            # Special checks for ttbar -> tau+jets -> e/mu+jets decay chain
            if (abs(genParticlesPdgId[igen]) == 15 and abs(genParticlesMom0ID[igen]) == 24 and (abs(genParticlesDau0ID[igen]) == 11 or abs(genParticlesDau1ID[igen]) == 11)) :
                tauToEl = True
                if options.debug :
                    print 'Tau decay to electron'
            if (abs(genParticlesPdgId[igen]) == 15 and abs(genParticlesMom0ID[igen]) == 24 and (abs(genParticlesDau0ID[igen]) == 13 or abs(genParticlesDau1ID[igen]) == 13)) :
                tauToMu = True
                if options.debug :
                    print 'Tau decay to muon'

            # Get muon
            if (abs(genParticlesPdgId[igen]) == 13 and (abs(genParticlesMom0ID[igen]) == 24 or (abs(genParticlesMom0ID[igen]) == 15 and tauToMu))) :
                isMuon = True
                p4Muon = ROOT.TLorentzVector()
                p4Muon.SetPtEtaPhiM( genParticlesPt[igen], genParticlesEta[igen], genParticlesPhi[igen], genParticlesMass[igen] )
                genMuons.append(p4Muon)
                if options.debug :
                    if abs(genParticlesMom0ID[igen]) == 15 and tauToMu :
                        print 'Got muon from tau'
                
            if (abs(genParticlesPdgId[igen]) == 11 and (abs(genParticlesMom0ID[igen]) == 24 or (abs(genParticlesMom0ID[igen]) == 15 and tauToEl))) :
                isElectron = True
                p4Electron = ROOT.TLorentzVector()
                p4Electron.SetPtEtaPhiM( genParticlesPt[igen], genParticlesEta[igen], genParticlesPhi[igen], genParticlesMass[igen] )
                genElectrons.append(p4Electron)
                if options.debug :
                    if abs(genParticlesMom0ID[igen]) == 15 and tauToEl :
                        print 'Got electron from tau'
                
        topQuarks.append( GenTopQuark( 6, p4Top, topDecay) )
        topQuarks.append( GenTopQuark( -6, p4Antitop, antitopDecay) )
        
        if (topDecay + antitopDecay == 1) and (isMuon == True) and (isElectron == False) :
            isSemiLeptonicGen = True
        elif (topDecay + antitopDecay == 1) and (isMuon == False) and (isElectron == True) :
            isSemiLeptonicGen = True
        else :
            isSemiLeptonicGen = False

        if options.semilep > 0 and isSemiLeptonicGen == False:
            continue
        if options.semilep < 0 and isSemiLeptonicGen == True:
            continue

        # Store channel for semileptonic events
        if options.semilep == 1:
            if isMuon:
                truthChannel.push_back(0)
            elif isElectron:
                truthChannel.push_back(1)
            else:
                truthChannel.push_back(2)
        
        if topDecay == 0 :
            hadTop = topQuarks[0]
            lepTop = topQuarks[1]
        else :
            hadTop = topQuarks[1]
            lepTop = topQuarks[0]

            
        # cut on generated m(ttbar) if stitching sample
        ttbarGen = hadTop.p4 + lepTop.p4
        mttbarGen = ttbarGen.M()
        
        if options.mttGenMax is not None :
            if mttbarGen > options.mttGenMax :
                continue

        # ------------------------------------------------------------
        # Store parton-level information if event good at parton level
        # ------------------------------------------------------------
        
        if hadTop.p4.Perp() > TOP_PT_CUT :
            nPassParton += 1
            filledEvent = True
            genTopPt.push_back(hadTop.p4.Perp())
            genTopEta.push_back(hadTop.p4.Eta())
            genTopPhi.push_back(hadTop.p4.Phi())
            genTTbarMass.push_back(ttbarGen.M())
            
            if len(genMuons) != 0:
                genMuPt.push_back(genMuons[0].Perp())
                genMuEta.push_back(genMuons[0].Eta())
                genMuPhi.push_back(genMuons[0].Phi())
                
            if len(genElectrons) != 0:
                genElPt.push_back(genElectrons[0].Perp())
                genElEta.push_back(genElectrons[0].Eta())
                genElPhi.push_back(genElectrons[0].Phi())   

    nPassSemiLep += 1
            
    # -------------------------------------------------------------------------------------
    # read gen jets
    # -------------------------------------------------------------------------------------
        
    if options.isMC and options.semilep is not None :

        # Get AK4 gen jet information
        #genBJets = []
        #event.getByLabel( ak4MatchedGenJetPtLabel, ak4MatchedGenJetPtHandle )
        #if ak4MatchedGenJetPtHandle.isValid() :
        #    event.getByLabel( ak4MatchedGenJetEtaLabel, ak4MatchedGenJetEtaHandle )
        #    event.getByLabel( ak4MatchedGenJetPhiLabel, ak4MatchedGenJetPhiHandle )
        #    event.getByLabel( ak4MatchedGenJetEnergyLabel, ak4MatchedGenJetEnergyHandle )
        #    
        #    ak4GenJetPt   = ak4MatchedGenJetPtHandle.product()
        #    ak4GenJetEta  = ak4MatchedGenJetEtaHandle.product()
        #    ak4GenJetPhi  = ak4MatchedGenJetPhiHandle.product()
        #    ak4GenJetE    = ak4MatchedGenJetEnergyHandle.product()
        #
        #    # loop over AK4 gen jets
        #    if len(ak4GenJetPt) != 0:
        #        for iak4 in xrange( len(ak4GenJetPt) ) :
        #            if ak4GenJetPt[iak4] > MIN_JET_PT and abs(ak4GenJetEta[iak4]) < MAX_JET_ETA :
        #                p4 = ROOT.TLorentzVector()
        #                p4.SetPtEtaPhiE( ak4GenJetPt[iak4], ak4GenJetEta[iak4], ak4GenJetPhi[iak4], ak4GenJetE[iak4] )
        #                genBJets.append(p4)

        # Get AK8 gen jet information
        genTops = []
        event.getByLabel( ak8GenJetPtLabel, ak8GenJetPtHandle )
        if ak8GenJetPtHandle.isValid() :
            event.getByLabel( ak8GenJetEtaLabel, ak8GenJetEtaHandle )
            event.getByLabel( ak8GenJetPhiLabel, ak8GenJetPhiHandle )
            event.getByLabel( ak8GenJetMassLabel, ak8GenJetMassHandle )
            
            ak8GenJetPt   = ak8GenJetPtHandle.product()
            ak8GenJetEta  = ak8GenJetEtaHandle.product()
            ak8GenJetPhi  = ak8GenJetPhiHandle.product()
            ak8GenJetMass = ak8GenJetMassHandle.product()
                
            # loop over AK8 gen jets
            if len(ak8GenJetPt) != 0:
                for iak8 in xrange( len(ak8GenJetPt) ) :
                    if ak8GenJetPt[iak8] > TOP_PT_CUT and abs(ak8GenJetEta[iak8]) < MAX_JET_ETA and ak8GenJetMass[iak8] > 140. and ak8GenJetMass[iak8] < 250.:
                        p4 = ROOT.TLorentzVector()
                        p4.SetPtEtaPhiM( ak8GenJetPt[iak8], ak8GenJetEta[iak8], ak8GenJetPhi[iak8], ak8GenJetMass[iak8] )
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
            filledEvent = True
            #if options.fullTruth:
            #    genAK8jetPt.push_back(genTops[0].Perp())
            #    genAK8jetEta.push_back(genTops[0].Eta())
            #    genAK8jetPhi.push_back(genTops[0].Phi())
            #    genAK8jetMass.push_back(genTops[0].M())
            #    if len(partMu) != 0:
            #        partMuPt.push_back(partMu[0].Perp())
            #    if len(partEl) != 0:
            #        partElPt.push_back(partEl[0].Perp())

            ##else :
            #for ibjet in genBJets :
            #    genAK4jetPt.push_back(ibjet.Perp())
            #    genAK4jetEta.push_back(ibjet.Eta())
            #    genAK4jetPhi.push_back(ibjet.Phi())
            #    genAK4jetMass.push_back(ibjet.M())

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
        event.getByLabel (electronD0Label, electronD0Handle)
        electronD0s = electronD0Handle.product()
        event.getByLabel (electronDzLabel, electronDzHandle)
        electronDzs = electronDzHandle.product()
        event.getByLabel (electronMiniIsoLabel, electronMiniIsoHandle)
        electronMiniIsos = electronMiniIsoHandle.product()

        event.getByLabel (elKeyLabel, elKeyHandle)
        elKeys = elKeyHandle.product()

        for ielectronPt in range(0,len(electronPts)) :
            electronPt = electronPts[ielectronPt]
            electronEta = electronEtas[ielectronPt]
            electronPhi = electronPhis[ielectronPt]
            electronMass = 0.0
            if (electronPt < MIN_EL_PT or abs(electronEta) > MAX_EL_ETA ) :
                continue
            nEleRaw += 1.0
            manualEisMedium = False
            if abs( electronEta ) <= 1.479 :
                if abs(electronDEtaIns[ielectronPt]) < 0.0103 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.0336 :
                        if electronFull5x5siees[ielectronPt] < 0.0101 :
                            if electronHoEs[ielectronPt] <  0.0876 :
                                if abs(electronD0s[ielectronPt]) < 0.0118 :
                                    if abs(electronDzs[ielectronPt]) <  0.373 :
                                        if electronooEmooPs[ielectronPt] <  0.0174 :
                                            manualEisMedium = True
            else :
                if abs(electronDEtaIns[ielectronPt]) < 0.00733 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.114 :
                        if electronFull5x5siees[ielectronPt] < 0.0283 :
                            if electronHoEs[ielectronPt] <  0.0678 :
                                if abs(electronD0s[ielectronPt]) < 0.0739 :
                                    if abs(electronDzs[ielectronPt]) <  0.602 :
                                        if electronooEmooPs[ielectronPt] <  0.0898 :
                                            manualEisMedium = True
                                            
            manualEisTight = False
            if abs( electronEta ) <= 1.479 :
                if abs(electronDEtaIns[ielectronPt]) < 0.00926 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.0336 :
                        if electronFull5x5siees[ielectronPt] < 0.0101 :
                            if electronHoEs[ielectronPt] <  0.0597 :
                                if abs(electronD0s[ielectronPt]) < 0.0111 :
                                    if abs(electronDzs[ielectronPt]) <  0.0466 :
                                        if electronooEmooPs[ielectronPt] <  0.012 :
                                            manualEisTight = True
            else :
                if abs(electronDEtaIns[ielectronPt]) < 0.00724 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.0918 :
                        if electronFull5x5siees[ielectronPt] < 0.0279 :
                            if electronHoEs[ielectronPt] <  0.0615 :
                                if abs(electronD0s[ielectronPt]) < 0.0351 :
                                    if abs(electronDzs[ielectronPt]) <  0.417 :
                                        if electronooEmooPs[ielectronPt] <  0.00999 :
                                            manualEisTight = True
                                            
            eIsTight = isElectronTight[ielectronPt]
            eIsMedium = isElectronMedium[ielectronPt]
            
            if eIsMedium and not manualEisMedium :
                nMediumNotManualEl += 1.0
            if manualEisMedium and not eIsMedium :
                nManualNotMediumEl += 1.0
            if eIsTight and not manualEisTight :
                nTightNotManualEl += 1.0
            if manualEisTight and not eIsTight :
                nManualNotTightEl += 1.0

            if not manualEisMedium :
                continue

            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( electronPt, electronEta, electronPhi, electronMass )

            elPt.push_back(electronPt)
            elEta.push_back(electronEta)
            elPhi.push_back(electronPhi)
            elMiniIso.push_back(electronMiniIsos[ielectronPt])

            elCand.append(p4)
            lepCand.append(p4)
            elCandKey.append(elKeys[ielectronPt])
            lepCandKey.append(elKeys[ielectronPt])
                
            if manualEisTight :
                elTight.push_back(1)
            else :
                elTight.push_back(0)
                
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
    # check that we have at least one lepton candidate
    # -------------------------------------------------------------------------------------

    if not ((passMuTrig and len(muCand) == 1 and len(elCand) == 0) or (passElTrig and len(elCand) == 1 and len(muCand) == 0)):
        passReco = False
        if not options.fullTruth :
            continue
    else :
        nPassLep += 1

    #Store which trigger passes for later use
    if options.isMC:
        if passMuTrig :
            muTrigPass.push_back(1)
        else :
            muTrigPass.push_back(0)

        if passElTrig :
            elTrigPass.push_back(1)
        else :
            elTrigPass.push_back(0)

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
        event.getByLabel (ak4JetMassLabel, ak4JetMassHandle)
        ak4JetMasss = ak4JetMassHandle.product()
        event.getByLabel (ak4JetCSVLabel, ak4JetCSVHandle)
        ak4JetCSVs = ak4JetCSVHandle.product()
        event.getByLabel (ak4JetVtxMassLabel, ak4JetVtxMassHandle)
        ak4JetVtxMasses = ak4JetVtxMassHandle.product()
        event.getByLabel (ak4JetAreaLabel, ak4JetAreaHandle)
        ak4JetAreas = ak4JetAreaHandle.product()

        # variables needed for jet ID 
        event.getByLabel( ak4JetNeuHadEnergyLabel, ak4JetNeuHadEnergyHandle )
        event.getByLabel( ak4JetNeuEmEnergyLabel, ak4JetNeuEmEnergyHandle )
        event.getByLabel( ak4JetChHadEnergyLabel, ak4JetChHadEnergyHandle )
        event.getByLabel( ak4JetChEmEnergyLabel, ak4JetChEmEnergyHandle )
        event.getByLabel( ak4JetNumDaughterLabel, ak4JetNumDaughterHandle )
        event.getByLabel( ak4JetChMultiLabel, ak4JetChMultiHandle )

        ak4JetNeuHadEnergys = ak4JetNeuHadEnergyHandle.product()
        ak4JetNeuEmEnergys = ak4JetNeuEmEnergyHandle.product()
        ak4JetChHadEnergys = ak4JetChHadEnergyHandle.product()
        ak4JetChEmEnergys = ak4JetChEmEnergyHandle.product()
        ak4JetNumDaughters = ak4JetNumDaughterHandle.product()
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
            jetP4Pre.SetPtEtaPhiM( ak4JetPts[ijet], ak4JetEtas[ijet], ak4JetPhis[ijet], ak4JetMasss[ijet] )

            # get the raw jet energy
            jetP4Raw = jetP4Pre * ak4JetJECs[ijet]

            # calculate the neutral/charged hadron/em energy fractions
            nhf = ak4JetNeuHadEnergys[ijet] / jetP4Raw.E()
            nef = ak4JetNeuEmEnergys[ijet] / jetP4Raw.E()
            chf = ak4JetChHadEnergys[ijet] / jetP4Raw.E()
            cef = ak4JetChEmEnergys[ijet] / jetP4Raw.E()
            nconstituents = ak4JetNumDaughters[ijet]
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
                    pfcands = int(ak4JetNumDaughters[ijet])
                    for ipf in range(0, pfcands) :
                        # if jet daughter matches lepton, remove lepton p4 from (raw) jet p4
                        if ak4JetKeys[ijet][ipf] in lepCandKey[ilep]: 
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
            else :
                if abs(jetP4.Perp() - jetP4Pre.Perp()) > 0.001 :
                    print 'Error: AK4 jet pt changed after undoing / redoing JEC w/o cleaning! Old pt ' + str(jetP4Pre.Perp()) + ' new pt ' + str(jetP4.Perp())
                    
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
                h_JERSF_AK4.Fill(ak4JERSFnoms[ijet])
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
                        jetP4 *= (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)*ak4JERs[ijet]))
                        jetP4_JERup   = jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFups[ijet]*ak4JERSFups[ijet]-1)*ak4JERs[ijet]))
                        jetP4_JERdown = jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFdowns[ijet]*ak4JERSFdowns[ijet]-1)*ak4JERs[ijet]))
                else:
                    jetP4 *= (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)*ak4JERs[ijet]))
                    jetP4_JERup   = jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFups[ijet]*ak4JERSFups[ijet]-1)*ak4JERs[ijet]))
                    jetP4_JERdown = jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak4JERSFdowns[ijet]*ak4JERSFdowns[ijet]-1)*ak4JERs[ijet]))

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
        event.getByLabel (ak8JetYLabel, ak8JetYHandle)
        ak8JetY = ak8JetYHandle.product()
        event.getByLabel (ak8JetMassLabel, ak8JetMassHandle)
        ak8JetMass = ak8JetMassHandle.product()
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
        event.getByLabel(sjSoftDropYLabel, sjSoftDropYHandle)
        sjSoftDropY    = sjSoftDropYHandle.product()
        event.getByLabel(sjSoftDropMassLabel, sjSoftDropMassHandle)
        sjSoftDropMass = sjSoftDropMassHandle.product()
        event.getByLabel(sjSoftDropCSVLabel, sjSoftDropCSVHandle)
        sjSoftDropCSV  = sjSoftDropCSVHandle.product()

        event.getByLabel (ak8JetAreaLabel, ak8JetAreaHandle)
        ak8JetAreas = ak8JetAreaHandle.product()
        event.getByLabel( ak8JetJECLabel, ak8JetJECHandle )
        ak8JetJECs = ak8JetJECHandle.product()

        event.getByLabel( ak8JetNeuHadEnergyLabel, ak8JetNeuHadEnergyHandle )
        event.getByLabel( ak8JetNeuEmEnergyLabel, ak8JetNeuEmEnergyHandle )
        event.getByLabel( ak8JetChHadEnergyLabel, ak8JetChHadEnergyHandle )
        event.getByLabel( ak8JetChEmEnergyLabel, ak8JetChEmEnergyHandle )
        event.getByLabel( ak8JetNumDaughterLabel, ak8JetNumDaughterHandle )
        event.getByLabel( ak8JetChMultiLabel, ak8JetChMultiHandle )

        ak8JetNeuHadEnergys = ak8JetNeuHadEnergyHandle.product()
        ak8JetNeuEmEnergys = ak8JetNeuEmEnergyHandle.product()
        ak8JetChHadEnergys = ak8JetChHadEnergyHandle.product()
        ak8JetChEmEnergys = ak8JetChEmEnergyHandle.product()
        ak8JetNumDaughters = ak8JetNumDaughterHandle.product()
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
            AK8jetP4Pre.SetPtEtaPhiM( ak8JetPt[ijet], ak8JetEta[ijet], ak8JetPhi[ijet], ak8JetMass[ijet] )

            # get the raw jet energy
            AK8jetP4Raw = AK8jetP4Pre * ak8JetJECs[ijet]

            # calculate the neutral/charged hadron/em energy fractions
            AK8nhf = ak8JetNeuHadEnergys[ijet] / AK8jetP4Raw.E()
            AK8nef = ak8JetNeuEmEnergys[ijet] / AK8jetP4Raw.E()
            AK8chf = ak8JetChHadEnergys[ijet] / AK8jetP4Raw.E()
            AK8cef = ak8JetChEmEnergys[ijet] / AK8jetP4Raw.E()
            AK8nconstituents = ak8JetNumDaughters[ijet]
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
                    pfcands = int(ak8JetNumDaughters[ijet])
                    for ipf in range(0, pfcands) :
                        # if jet daughter matches lepton, remove lepton p4 from (raw) jet p4
                        if ak8JetKeys[ijet][ipf] in lepCandKey[ilep]:
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
            else :
                if abs(AK8jetP4.Perp() - AK8jetP4Pre.Perp()) > 0.001 :
                    print 'Error: AK8 jet pt changed after undoing / redoing JEC w/o cleaning! Old pt ' + str(jetP4Pre.Perp()) + ' new pt ' + str(jetP4.Perp())

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
                h_JERSF_AK8.Fill(ak8JERSFnoms[ijet])
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
                        AK8jetP4 *= (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)*ak8JERs[ijet]))
                        AK8jetP4_JERup   = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFups[ijet]*ak8JERSFups[ijet]-1)*ak8JERs[ijet]))
                        AK8jetP4_JERdown = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFdowns[ijet]*ak8JERSFdowns[ijet]-1)*ak8JERs[ijet]))
                else:
                    AK8jetP4 *= (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)*ak8JERs[ijet]))
                    AK8jetP4_JERup   = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFups[ijet]*ak8JERSFups[ijet]-1)*ak8JERs[ijet]))
                    AK8jetP4_JERdown = AK8jetP4 * (1.0 + smearfunc.Gaus(0.0,math.sqrt(ak8JERSFdowns[ijet]*ak8JERSFdowns[ijet]-1)*ak8JERs[ijet]))

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
                ak8jetSDsubjet0pt.push_back(sjSoftDropPt[int(ak8JetSDSj0[ijet])])
                ak8jetSDsubjet0eta.push_back(sjSoftDropEta[int(ak8JetSDSj0[ijet])])
                ak8jetSDsubjet0phi.push_back(sjSoftDropPhi[int(ak8JetSDSj0[ijet])])
                ak8jetSDsubjet0mass.push_back(sjSoftDropMass[int(ak8JetSDSj0[ijet])])
                ak8jetSDsubjet0CSV.push_back(sjSoftDropCSV[int(ak8JetSDSj0[ijet])])
            else :
                ak8jetSDsubjet0pt.push_back(-10.0)
                ak8jetSDsubjet0eta.push_back(-10.0)
                ak8jetSDsubjet0phi.push_back(-10.0)
                ak8jetSDsubjet0mass.push_back(-10.0)
                ak8jetSDsubjet0CSV.push_back(-10.0)
            if int(ak8JetSDSj1[ijet]) >= 0 and int(ak8JetSDSj1[ijet]) < len(sjSoftDropPt) :
                ak8jetSDsubjet1pt.push_back(sjSoftDropPt[int(ak8JetSDSj1[ijet])])
                ak8jetSDsubjet1eta.push_back(sjSoftDropEta[int(ak8JetSDSj1[ijet])])
                ak8jetSDsubjet1phi.push_back(sjSoftDropPhi[int(ak8JetSDSj1[ijet])])
                ak8jetSDsubjet1mass.push_back(sjSoftDropMass[int(ak8JetSDSj1[ijet])])
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
        filledEvent = True

    if options.fullTruth : # Only store information relevant for unfolding -- differential quantities and things needed for selection

        ak8jetTau1.clear() 
        ak8jetMassPruned.clear()
        ak8jetMassFiltered.clear()
        ak8jetMassTrimmed.clear()   
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
    
    if filledEvent :
        eventWeight_nom.push_back(weight_nom)
        if options.isMC:
            eventWeight_puUp.push_back(weight_puUp)
            eventWeight_puDown.push_back(weight_puDown)
        myTree.Fill()
    
# -------------------------------------------------------------------------------------
# END OF LOOPING OVER EVENTS!!!
# -------------------------------------------------------------------------------------

print 'Total Events:    ' + str(ntotal)
print '-------------------------------'
print 'Pass nPV:        ' + str(nPassNPV)
print 'Pass MET filter: ' + str(nPassMetFilter)
print 'Pass rho:        ' + str(nPassRho)
print 'Pass semilep:    ' + str(nPassSemiLep)
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
print 'Fraction of all medium electrons not passing twiki cuts: ' + str(nMediumNotManualEl / nEleRaw)
print 'Fraction of all electrons passing twiki cuts not medium: ' + str(nManualNotMediumEl / nEleRaw)
print 'Fraction of all tight electrons not passing twiki cuts:  ' + str(nTightNotManualEl / nEleRaw)
print 'Fraction of all electrons passing twiki cuts not tight:  ' + str(nManualNotTightEl / nEleRaw)
print 'Fraction of all tight muons not passing twiki cuts:      ' + str(nTightNotManualMu / nMuRaw)
print 'Fraction of all muons passing twiki cuts not tight:      ' + str(nManualNotTightMu / nMuRaw)
print '-------------------------------'

f.cd()
f.Write()
f.Close()

