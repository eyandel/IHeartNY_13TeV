import FWCore.ParameterSet.Config as cms

process = cms.Process('NoSplit')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring([
'/store/group/lpctlbsm/B2GAnaFW_80X_V1p0/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160609_173944/0000/B2GEDMNtuple_1.root',
'/store/group/lpctlbsm/B2GAnaFW_80X_V1p0/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160609_173944/0000/B2GEDMNtuple_10.root',
'/store/group/lpctlbsm/B2GAnaFW_80X_V1p0/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160609_173944/0000/B2GEDMNtuple_100.root',
            ]))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
