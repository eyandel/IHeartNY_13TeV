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

TOP_PT_CUT = 400.0

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

parser.add_option('--debug', metavar='M', action='store_true',
                  default=False,
                  dest='debug',
                  help='Print out debug statements')

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
from DataFormats.FWLite import Events, Handle, Runs
    
# -------------------------------------------------------------------------------------
# jet energy corrections
# -------------------------------------------------------------------------------------

ROOT.gSystem.Load('libCondFormatsJetMETObjects')

L3JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt')
L3JetParAK4 = ROOT.JetCorrectorParameters(L3JecStrAK4)
L2JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt')
L2JetParAK4 = ROOT.JetCorrectorParameters(L2JecStrAK4)
L1JecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt')
L1JetParAK4 = ROOT.JetCorrectorParameters(L1JecStrAK4)
UncJecStrAK4 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt')
UncertJetAK4 = ROOT.JetCorrectionUncertainty(UncJecStrAK4)
L3JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt')
L3JetParAK8 = ROOT.JetCorrectorParameters(L3JecStrAK8)
L2JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt')
L2JetParAK8 = ROOT.JetCorrectorParameters(L2JecStrAK8)
L1JecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt')
L1JetParAK8 = ROOT.JetCorrectorParameters(L1JecStrAK8)
UncJecStrAK8 = ROOT.std.string('JECs/Summer16_23Sep2016V3_MC_Uncertainty_AK8PFchs.txt')
UncertJetAK8 = ROOT.JetCorrectionUncertainty(UncJecStrAK8)
    
#  load JetCorrectorParameter objects into vector (order matters!)
vParJecAK4 = ROOT.std.vector(ROOT.JetCorrectorParameters)()
vParJecAK4.push_back(L1JetParAK4)
vParJecAK4.push_back(L2JetParAK4)
vParJecAK4.push_back(L3JetParAK4)

ak4JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK4)

vParJecAK8 = ROOT.std.vector(ROOT.JetCorrectorParameters)()
vParJecAK8.push_back(L1JetParAK8)
vParJecAK8.push_back(L2JetParAK8)
vParJecAK8.push_back(L3JetParAK8)

ak8JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK8)

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

# reco level
elPt                   = ROOT.vector('float')()
elEta                  = ROOT.vector('float')()
elPhi                  = ROOT.vector('float')()
elFull5x5siee          = ROOT.vector('float')()
eldEtaIn               = ROOT.vector('float')()
eldPhiIn               = ROOT.vector('float')()
elHoE                  = ROOT.vector('float')()
elooEmooP              = ROOT.vector('float')()
elSCEta                = ROOT.vector('float')()
elMissHits             = ROOT.vector('float')()
elConVeto              = ROOT.vector('float')()
elDxy                  = ROOT.vector('float')()
elDz                   = ROOT.vector('float')()
elMiniIso              = ROOT.vector('float')()
ak4jetPt               = ROOT.vector('float')()
ak4jetEta              = ROOT.vector('float')()
ak4jetPhi              = ROOT.vector('float')()
ak4jetMass             = ROOT.vector('float')()
ak8jetPt               = ROOT.vector('float')()
ak8jetEta              = ROOT.vector('float')()
ak8jetPhi              = ROOT.vector('float')()
ak8jetMass             = ROOT.vector('float')()
metPt                  = ROOT.vector('float')()
metPhi                 = ROOT.vector('float')()
eventWeight_nom        = ROOT.vector('float')()

recoTree.Branch('elPt'                   , elPt                   )
recoTree.Branch('elEta'                  , elEta                  )
recoTree.Branch('elPhi'                  , elPhi                  )
recoTree.Branch('elFull5x5siee'          , elFull5x5siee          )
recoTree.Branch('eldEtaIn'               , eldEtaIn               )
recoTree.Branch('eldPhiIn'               , eldPhiIn               )
recoTree.Branch('elHoE'                  , elHoE                  )
recoTree.Branch('elooEmooP'              , elooEmooP              )
recoTree.Branch('elSCEta'                , elSCEta                )
recoTree.Branch('elMissHits'             , elMissHits             )
recoTree.Branch('elConVeto'              , elConVeto              )
recoTree.Branch('elDxy'                  , elDxy                  )
recoTree.Branch('elDz'                   , elDz                   )
recoTree.Branch('elMiniIso'              , elMiniIso              )
recoTree.Branch('ak4jetPt'               , ak4jetPt               )
recoTree.Branch('ak4jetEta'              , ak4jetEta              )
recoTree.Branch('ak4jetPhi'              , ak4jetPhi              )
recoTree.Branch('ak4jetMass'             , ak4jetMass             )
recoTree.Branch('ak8jetPt'               , ak8jetPt               )
recoTree.Branch('ak8jetEta'              , ak8jetEta              )
recoTree.Branch('ak8jetPhi'              , ak8jetPhi              )
recoTree.Branch('ak8jetMass'             , ak8jetMass             )
recoTree.Branch('metPt'                  , metPt                  )
recoTree.Branch('metPhi'                 , metPhi                 )
recoTree.Branch('eventWeight_nom'        , eventWeight_nom        )

# -------------------------------------------------------------------------------------
# define all variables to be read from input files
# -------------------------------------------------------------------------------------

events = Events (files)
runs = Runs ( files)

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
muMediumHandle    = Handle("std::vector<float>")
muMediumLabel     = ("muons" , "muIsMediumMuon" )

electronPtHandle          = Handle( "std::vector<float>")
electronPtLabel           = ("electrons", "elPt")
electronEtaHandle         = Handle( "std::vector<float>")
electronEtaLabel          = ("electrons", "elEta")
electronPhiHandle         = Handle( "std::vector<float>")
electronPhiLabel          = ("electrons", "elPhi")
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
electronDxyHandle         = Handle("std::vector<float>")
electronDxyLabel          = ( "electrons" , "elDxy")
electronDzHandle          = Handle("std::vector<float>")
electronDzLabel           = ( "electrons" , "elDz")

elKeyHandle = Handle("std::vector<std::vector<int> >")
elKeyLabel = ( "electronKeys" )

# AK4 jet collection
ak4JetPtHandle   = Handle( "std::vector<float>" )
ak4JetPtLabel    = ("jetsAK4CHS", "jetAK4CHSPt")
ak4JetEtaHandle  = Handle( "std::vector<float>" )
ak4JetEtaLabel   = ("jetsAK4CHS", "jetAK4CHSEta")
ak4JetPhiHandle  = Handle( "std::vector<float>" )
ak4JetPhiLabel   = ("jetsAK4CHS", "jetAK4CHSPhi")
ak4JetEHandle = Handle( "std::vector<float>" )
ak4JetELabel  = ("jetsAK4CHS", "jetAK4CHSE")
ak4JetCSVHandle  = Handle( "std::vector<float>" )
ak4JetCSVLabel   = ("jetsAK4CHS", "jetAK4CHSCSVv2")
ak4JetVtxMassHandle = Handle( "std::vector<float>" )
ak4JetVtxMassLabel  = ("jetsAK4CHS", "jetAK4CHSSV0mass")    
ak4JetAreaHandle = Handle( "std::vector<float>" )
ak4JetAreaLabel  = ("jetsAK4CHS", "jetAK4CHSjetArea")

ak4JetNeuHadEnergyFracHandle = Handle("std::vector<float>")
ak4JetNeuHadEnergyFracLabel = ("jetsAK4CHS" , "jetAK4CHSneutralHadronEnergyFrac")
ak4JetNeuEmEnergyFracHandle = Handle("std::vector<float>")
ak4JetNeuEmEnergyFracLabel = ("jetsAK4CHS" , "jetAK4CHSneutralEmEnergyFrac")
ak4JetChHadEnergyFracHandle = Handle("std::vector<float>")
ak4JetChHadEnergyFracLabel = ("jetsAK4CHS" , "jetAK4CHSchargedHadronEnergyFrac")
ak4JetChEmEnergyFracHandle = Handle("std::vector<float>")
ak4JetChEmEnergyFracLabel = ("jetsAK4CHS" , "jetAK4CHSchargedEmEnergyFrac")
ak4JetNeuMultiHandle = Handle("std::vector<float>")
ak4JetNeuMultiLabel = ("jetsAK4CHS" , "jetAK4CHSneutralMultiplicity")
ak4JetChMultiHandle = Handle("std::vector<float>")
ak4JetChMultiLabel = ("jetsAK4CHS" , "jetAK4CHSchargedMultiplicity")

ak4JetJECHandle = Handle("std::vector<float>")
ak4JetJECLabel = ("jetsAK4CHS" , "jetAK4CHSjecFactor0")

# JER variables
ak4MatchedGenJetPtHandle   = Handle("std::vector<float>")
ak4MatchedGenJetPtLabel    = ("jetsAK4CHS", "jetAK4CHSGenJetPt")
ak4MatchedGenJetEtaHandle  = Handle("std::vector<float>")
ak4MatchedGenJetEtaLabel   = ("jetsAK4CHS", "jetAK4CHSGenJetEta")
ak4MatchedGenJetPhiHandle  = Handle("std::vector<float>")
ak4MatchedGenJetPhiLabel   = ("jetsAK4CHS", "jetAK4CHSGenJetPhi")
ak4MatchedGenJetEnergyHandle = Handle("std::vector<float>")
ak4MatchedGenJetEnergyLabel  = ("jetsAK4CHS", "jetAK4CHSGenJetE")
ak4JERSFnomHandle  = Handle("std::vector<float>")
ak4JERSFnomLabel   = ("jetsAK4CHS", "jetAK4CHSJERSF")
ak4JERHandle       = Handle("std::vector<float>")
ak4JERLabel        = ("jetsAK4CHS", "jetAK4CHSPtResolution")

ak4JetKeysHandle = Handle("std::vector<std::vector<int> >")
ak4JetKeysLabel = ( "jetKeysAK4CHS" , "" )

# AK8 jet collection
ak8JetPtHandle   = Handle( "std::vector<float>" )
ak8JetPtLabel    = ("jetsAK8CHS", "jetAK8CHSPt")
ak8JetEtaHandle  = Handle( "std::vector<float>" )
ak8JetEtaLabel   = ("jetsAK8CHS", "jetAK8CHSEta")
ak8JetPhiHandle  = Handle( "std::vector<float>" )
ak8JetPhiLabel   = ("jetsAK8CHS", "jetAK8CHSPhi")
ak8JetEHandle = Handle( "std::vector<float>" )
ak8JetELabel  = ("jetsAK8CHS", "jetAK8CHSE")

# JETID variables
ak8JetNeuHadEnergyFracHandle = Handle("std::vector<float>")
ak8JetNeuHadEnergyFracLabel = ("jetsAK8CHS" , "jetAK8CHSneutralHadronEnergyFrac")
ak8JetNeuEmEnergyFracHandle = Handle("std::vector<float>")
ak8JetNeuEmEnergyFracLabel = ("jetsAK8CHS" , "jetAK8CHSneutralEmEnergyFrac")
ak8JetChHadEnergyFracHandle = Handle("std::vector<float>")
ak8JetChHadEnergyFracLabel = ("jetsAK8CHS" , "jetAK8CHSchargedHadronEnergyFrac")
ak8JetChEmEnergyFracHandle = Handle("std::vector<float>")
ak8JetChEmEnergyFracLabel = ("jetsAK8CHS" , "jetAK8CHSchargedEmEnergyFrac")
ak8JetNeuMultiHandle = Handle("std::vector<float>")
ak8JetNeuMultiLabel = ("jetsAK8CHS" , "jetAK8CHSneutralMultiplicity")
ak8JetChMultiHandle = Handle("std::vector<float>")
ak8JetChMultiLabel = ("jetsAK8CHS" , "jetAK8CHSchargedMultiplicity")

#JEC
ak8JetJECHandle = Handle("std::vector<float>")
ak8JetJECLabel = ("jetsAK8CHS" , "jetAK8CHSjecFactor0")
ak8JetAreaHandle = Handle( "std::vector<float>" )
ak8JetAreaLabel  = ("jetsAK8CHS", "jetAK8CHSjetArea")

#JER
ak8MatchedGenJetPtHandle   = Handle("std::vector<float>")
ak8MatchedGenJetPtLabel    = ("jetsAK8CHS", "jetAK8CHSGenJetPt")
ak8MatchedGenJetEtaHandle  = Handle("std::vector<float>")
ak8MatchedGenJetEtaLabel   = ("jetsAK8CHS", "jetAK8CHSGenJetEta")
ak8MatchedGenJetPhiHandle  = Handle("std::vector<float>")
ak8MatchedGenJetPhiLabel   = ("jetsAK8CHS", "jetAK8CHSGenJetPhi")
ak8MatchedGenJetEnergyHandle = Handle("std::vector<float>")
ak8MatchedGenJetEnergyLabel  = ("jetsAK8CHS", "jetAK8CHSGenJetE")
ak8JERSFnomHandle  = Handle("std::vector<float>")
ak8JERSFnomLabel   = ("jetsAK8CHS", "jetAK8CHSJERSF")
ak8JERHandle       = Handle("std::vector<float>")
ak8JERLabel        = ("jetsAK8CHS", "jetAK8CHSPtResolution")

ak8JetKeysHandle = Handle("std::vector<std::vector<int> >")
ak8JetKeysLabel = ( "jetKeysAK8CHS" , "" )

rhoHandle = Handle("double")
rhoLabel = ("fixedGridRhoFastjetAll", "")

puNtrueIntHandle = Handle("std::int")
puNtrueIntLabel = ( "eventUserData" , "puNtrueInt" )

runnumHandle = Handle("uint")
runnumLabel = ("eventInfo", "evtInfoRunNumber")

# -------------------------------------------------------------------------------------
# Get pileup weights
# -------------------------------------------------------------------------------------

fPUweight      = ROOT.TFile(options.puFile)
hPUweight_nom  = fPUweight.Get("PUweight_true")

# -------------------------------------------------------------------------------------
# reset various counters
# -------------------------------------------------------------------------------------

ntotal = 0       # total number of events

# Event quality
nPassNPV = 0
nPassLHE = 0
nPassMetFilter = 0
nPassRho = 0

nPassLep = 0
nPassAK4jet = 0
nPassAK8jet = 0
nEventsPass = 0

smearfunc = ROOT.TRandom3()

# -------------------------------                                                                                                               
# Get MET filter names, indices of mu and el triggers   
# -------------------------------
elTrigIndices = {}
metFiltIndices = {}

for run in runs:

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
        itrig += 1

# -------------------------------------------------------------------------------------
# start looping over events
# -------------------------------------------------------------------------------------

print "Start looping over events!"

for event in events :

    elPt.clear()
    elEta.clear()
    elPhi.clear()
    elFull5x5siee.clear()
    eldEtaIn.clear()
    eldPhiIn.clear()
    elHoE.clear()
    elooEmooP.clear()
    elSCEta.clear()
    elMissHits.clear()
    elConVeto.clear()
    elDxy.clear()
    elDz.clear()
    elMiniIso.clear()
    ak4jetPt.clear()
    ak4jetEta.clear()
    ak4jetPhi.clear()
    ak4jetMass.clear()
    ak8jetPt.clear()
    ak8jetEta.clear()
    ak8jetPhi.clear()
    ak8jetMass.clear()
    metPt.clear()
    metPhi.clear()
    eventWeight_nom.clear()

    weight_nom = 1.0 #event weight
    
    if ntotal % 1000 == 0 :
      print  '--------- Processing Event ' + str(ntotal)
    ntotal += 1

    if options.debug :
        if ntotal > 1000:
            continue
        print 'Event ' + str(ntotal)

    event.getByLabel(runnumLabel,runnumHandle)
    runnumber = runnumHandle.product()[0]

    # -------------------
    #
    # R E C O   L E V E L
    #
    # -------------------

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
    
    event.getByLabel(puNtrueIntLabel, puNtrueIntHandle)
    puNTrueInt = puNtrueIntHandle.product()[0] 

    if options.puFile is not None :
        weight_nom  *= hPUweight_nom.GetBinContent( hPUweight_nom.GetXaxis().FindBin( puNTrueInt ) )
        if options.debug :
            print 'Event weight after pileup reweigting is : ' + str(weight_nom)

    if NPV == 0 :
        continue
    nPassNPV += 1

    # -------------------------------
    # Require met filter
    # -------------------------------
    
    metFilt = True
    gotBits = event.getByLabel( metFilterBitsLabel, metFilterBitsHandle )

    if gotBits == False  :
        print 'No MET filter bits!'
        continue

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
        continue
    nPassMetFilter += 1
                
    # -------------------------------
    # Require a trigger
    # -------------------------------

    passElTrig = False
    prescale = 1.0
    
    event.getByLabel( trigBitsLabel, trigBitsHandle )
    event.getByLabel( trigPrescalesLabel, trigPrescalesHandle )
    
    triggerBits = trigBitsHandle.product()
    triggerPrescales = trigPrescalesHandle.product()
    
    if triggerBits[elTrigIndices[runnumber]] == 1 :
        passElTrig = True
        prescale = prescale * triggerPrescales[elTrigIndices[runnumber]]
            
    weight_nom = weight_nom * prescale
        
    if not passElTrig :
        continue
    
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
    # get electrons
    # -------------------------------------------------------------------------------------

    elCand = []
    elCandKey = []

    event.getByLabel (electronPtLabel, electronPtHandle)
    if electronPtHandle.isValid() :
        electronPts = electronPtHandle.product()
        event.getByLabel (electronEtaLabel, electronEtaHandle)
        electronEtas = electronEtaHandle.product()
        event.getByLabel (electronPhiLabel, electronPhiHandle)
        electronPhis = electronPhiHandle.product()
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
        event.getByLabel (electronDxyLabel, electronDxyHandle)
        electronDxys = electronDxyHandle.product()
        event.getByLabel (electronDzLabel, electronDzHandle)
        electronDzs = electronDzHandle.product()

        event.getByLabel (elKeyLabel, elKeyHandle)
        elKeys = elKeyHandle.product()

        for ielectronPt in range(0,len(electronPts)) :
            electronPt = electronPts[ielectronPt]
            electronEta = electronEtas[ielectronPt]
            electronPhi = electronPhis[ielectronPt]
            electronMass = 0.0
            if (electronPt < MIN_EL_PT or abs(electronEta) > MAX_EL_ETA ) :
                continue

            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( electronPt, electronEta, electronPhi, electronMass )

            elPt.push_back(electronPt)
            elEta.push_back(electronEta)
            elPhi.push_back(electronPhi)
            elFull5x5siee.push_back(electronFull5x5siees[ielectronPt])
            eldEtaIn.push_back(electronDEtaIns[ielectronPt])
            eldPhiIn.push_back(electronDPhiIns[ielectronPt])
            elHoE.push_back(electronHoEs[ielectronPt])
            elooEmooP.push_back(electronooEmooPs[ielectronPt])
            elSCEta.push_back(electronSCEtas[ielectronPt])
            elMissHits.push_back(electronMissingHits[ielectronPt])
            elConVeto.push_back(isElectronConVeto[ielectronPt])
            elDxy.push_back(electronDxys[ielectronPt])
            elDz.push_back(electronDzs[ielectronPt])
            elMiniIso.push_back(electronMiniIsos[ielectronPt])

            manualEisMedium = False
            if abs( electronSCEtas[ielectronPt] ) <= 1.479 :
                if abs(electronDEtaIns[ielectronPt]) < 0.00311 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.103 :
                        if electronFull5x5siees[ielectronPt] < 0.00998 :
                            if electronHoEs[ielectronPt] <  0.253 :
                                if electronooEmooPs[ielectronPt] <  0.134 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        if not isElectronConVeto[ielectronPt] :
                                            manualEisMedium = True
            else :
                if abs(electronDEtaIns[ielectronPt]) < 0.00609 :
                    if abs(electronDPhiIns[ielectronPt] ) < 0.045 :
                        if electronFull5x5siees[ielectronPt] < 0.0298 :
                            if electronHoEs[ielectronPt] <  0.0878 :
                                if electronooEmooPs[ielectronPt] <  0.13 :
                                    if electronMissingHits[ielectronPt] <= 1:
                                        if not isElectronConVeto[ielectronPt] :
                                            manualEisMedium = True

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

            if manualEisMedium and passD0 and passDz :
                elCand.append(p4)
                elCandKey.append(elKeys[ielectronPt])
                
    # --------------------------
    # get muons
    # --------------------------

    muCand = []

    event.getByLabel (muonPtLabel, muonPtHandle)
    if muonPtHandle.isValid() : 
        muonPts = muonPtHandle.product()
        event.getByLabel (muonEtaLabel, muonEtaHandle)
        muonEtas = muonEtaHandle.product()
        event.getByLabel (muonPhiLabel, muonPhiHandle)
        muonPhis = muonPhiHandle.product()
        event.getByLabel (muMediumLabel, muMediumHandle)
        isMediumMuon = muMediumHandle.product()

        for imuonPt in range(0,len(muonPts)) :
            muonPt = muonPts[imuonPt]
            muonEta = muonEtas[imuonPt]
            muonPhi = muonPhis[imuonPt]
            muonMass = 0.105658
            if muonPt < MIN_MU_PT or abs(muonEta) > MAX_MU_ETA :
                continue

            muIsMedium = isMediumMuon[imuonPt]
            if not muIsMedium :
                continue
            
            p4 = ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( muonPt, muonEta, muonPhi, muonMass )
            
            muCand.append(p4)
                
    # -------------------------------------------------------------------------------------
    # check that we have exactly one lepton candidate
    # -------------------------------------------------------------------------------------

    if not (elPt.size() >= 1 and len(muCand) == 0): # Veto on Medium ID muons, and at least one PF electron passing pt/eta cuts
        continue
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
        event.getByLabel(ak4JERLabel, ak4JERHandle)
        ak4JERs = ak4JERHandle.product()

        # jet constituents 
        event.getByLabel( ak4JetKeysLabel, ak4JetKeysHandle )
        ak4JetKeys = ak4JetKeysHandle.product()
        
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
            for ilep in xrange( len(elCand) ) :
                if elCand[ilep].DeltaR(jetP4Raw) < 0.4 :
                    ak4daughters = ak4JetKeys[ijet]
                    for daughterKey in ak4daughters :
                        if daughterKey in elCandKey[ilep]:
                            if options.debug:
                                print 'Event {0:d}, removing lepton with pt/eta/phi = {1:6.2f},{2:6.2f},{3:6.2f} from AK4 jet with pt/eta/phi = {4:6.2f},{5:6.2f},{6:6.2f}'.format(ntotal,elCand[ilep].Perp(),elCand[ilep].Eta(),elCand[ilep].Phi(),jetP4Raw.Perp(), jetP4Raw.Eta(), jetP4Raw.Phi())
                            if elCand[ilep].E() > jetP4Raw.E() :
                                zeroedEnergy = True
                            cleanedLepton = elCand[ilep]
                            jetP4Raw -= elCand[ilep]
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

            metv += jetP4 #Now correct MET for JER
            
            # Scale jet pt if there is a matched gen jet
            smearfacAK4 = smearfunc.Gaus(0.0,ak4JERs[ijet])
            if ak4MatchedGenJetPts[ijet] > 0:
                genJetP4 = ROOT.TLorentzVector()
                genJetP4.SetPtEtaPhiE(ak4MatchedGenJetPts[ijet],ak4MatchedGenJetEtas[ijet],ak4MatchedGenJetPhis[ijet],ak4MatchedGenJetEnergys[ijet])
                if jetP4.DeltaR(genJetP4) < 0.2 and abs(jetP4.Perp() - genJetP4.Perp()) < (3 * ak4JERs[ijet] * jetP4.Perp()) : # Do matching requirement
                    jetP4 -= genJetP4
                    jetP4 *= ak4JERSFnoms[ijet]
                    jetP4 += genJetP4
                else:
                    jetP4 *= (1.0 + smearfacAK4*math.sqrt(max(0.0,ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)))
            else:
                jetP4 *= (1.0 + smearfacAK4*math.sqrt(max(0.0,ak4JERSFnoms[ijet]*ak4JERSFnoms[ijet]-1)))
            
            metv -= jetP4

            if jetP4.Perp() < MIN_JET_PT or abs(jetP4.Eta()) > MAX_JET_ETA:
                continue

            ak4jetPt.push_back(jetP4.Perp())
            ak4jetEta.push_back(jetP4.Eta())
            ak4jetPhi.push_back(jetP4.Phi())
            ak4jetMass.push_back(jetP4.M())

    # -------------------------------------------------------
    # If not storing full truth information, require AK4 jet
    # -------------------------------------------------------

    if len(ak4jetPt) == 0 :
        continue

    nPassAK4jet += 1

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
            for ilep in xrange( len(elCand) ) :
                if elCand[ilep].DeltaR(AK8jetP4Raw) < 0.4 :
                    ak8daughters = ak8JetKeys[ijet]
                    for daughterKey in ak8daughters :
                        if daughterKey in elCandKey[ilep]:
                            if options.debug:
                                print 'Event {0:d}, removing lepton with pt/eta/phi = {1:6.2f},{2:6.2f},{3:6.2f} from AK8 jet with pt/eta/phi = {4:6.2f},{5:6.2f},{6:6.2f}'.format(ntotal,elCand[ilep].Perp(),elCand[ilep].Eta(),elCand[ilep].Phi(),AK8jetP4Raw.Perp(), AK8jetP4Raw.Eta(), AK8jetP4Raw.Phi())
                            if elCand[ilep].E() > AK8jetP4Raw.E() :
                                zeroedEnergy = True
                            AK8cleanedLepton = elCand[ilep]
                            AK8jetP4Raw -= elCand[ilep]
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

            metv += AK8jetP4 #Now correct MET for JER
            
            # Scale jet pt if there is a matched gen jet
            smearfacAK8 = smearfunc.Gaus(0.0,ak8JERs[ijet])
            if ak8MatchedGenJetPts[ijet] > 0:
                AK8genJetP4 = ROOT.TLorentzVector()
                AK8genJetP4.SetPtEtaPhiE(ak8MatchedGenJetPts[ijet],ak8MatchedGenJetEtas[ijet],ak8MatchedGenJetPhis[ijet],ak8MatchedGenJetEnergys[ijet])
                if AK8jetP4.DeltaR(AK8genJetP4) < 0.2 and abs(AK8jetP4.Perp() - AK8genJetP4.Perp()) < (3 * ak8JERs[ijet] * AK8jetP4.Perp()) : # Do matching requirement
                    AK8jetP4 -= AK8genJetP4
                    AK8jetP4 *= ak8JERSFnoms[ijet]
                    AK8jetP4 += AK8genJetP4
                else:
                    AK8jetP4 *= (1.0 + smearfacAK8*math.sqrt(max(0.0,ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)))
            else:
                AK8jetP4 *= (1.0 + smearfacAK8*math.sqrt(max(0.0,ak8JERSFnoms[ijet]*ak8JERSFnoms[ijet]-1)))
                
            metv -= AK8jetP4

            if (AK8jetP4.Perp() < TOP_PT_CUT or abs(AK8jetP4.Eta()) > MAX_JET_ETA):
                continue

            ak8jetPt.push_back(AK8jetP4.Perp())
            ak8jetEta.push_back(AK8jetP4.Eta())
            ak8jetPhi.push_back(AK8jetP4.Phi())
            ak8jetMass.push_back(AK8jetP4.M())

    if len(ak8jetPt) == 0 :
        continue
    nPassAK8jet += 1

    #if metv.Perp() < 35.0 :
    #    continue
    nEventsPass += 1
        
    eventWeight_nom.push_back(weight_nom)
    metPt.push_back(metv.Perp())
    metPhi.push_back(metv.Phi())

    recoTree.Fill()
    
# -------------------------------------------------------------------------------------
# END OF LOOPING OVER EVENTS!!!
# -------------------------------------------------------------------------------------

print 'Total Events:    ' + str(ntotal)
print '-------------------------------'
print 'Pass nPV:        ' + str(nPassNPV)
print 'Pass MET filter: ' + str(nPassMetFilter)
print 'Pass rho:        ' + str(nPassRho)
print '-------------------------------'
print 'Pass lepton:     ' + str(nPassLep)
print 'Pass AK4 jet:    ' + str(nPassAK4jet)
print 'Pass AK8 jet:    ' + str(nPassAK8jet)
print 'Pass reco:       ' + str(nEventsPass)
    
f.cd()
f.Write()
f.Close()

