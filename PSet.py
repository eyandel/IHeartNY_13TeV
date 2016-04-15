import FWCore.ParameterSet.Config as cms

process = cms.Process('NoSplit')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring([
     '/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_1.root',
'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_10.root',
'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_100.root',
            ]))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
