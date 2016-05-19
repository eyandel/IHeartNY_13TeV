import FWCore.ParameterSet.Config as cms

process = cms.Process('NoSplit')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring([
     #'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_1.root',
#'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_10.root',
#'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0000/B2GEDMNtuple_100.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/SingleElectron/B2GAnaFW_76X_V1p1_Run2015D-16Dec2015-v1/160401_164525/0000/B2GEDMNtuple_1.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/SingleMuon/B2GAnaFW_76X_V1p1_Run2015C_25ns-16Dec2015-v1/160401_164545/0000/B2GEDMNtuple_10.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/SingleElectron/B2GAnaFW_76X_V1p1_Run2015D-16Dec2015-v1/160401_164525/0000/B2GEDMNtuple_1.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/SingleMuon/B2GAnaFW_76X_V1p1_Run2015D-16Dec2015-v1/160401_164605/0000/B2GEDMNtuple_1.root',
#'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0001/B2GEDMNtuple_1997.root', 
#'/store/group/phys_b2g/B2GAnaFW_76X_V1p2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/B2GAnaFW_RunIIFall15MiniAODv2_25ns_v76x_v1p2/160408_145006/0001/B2GEDMNtuple_1998.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_76X_V1p1_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160401_100101/0000/B2GEDMNtuple_1.root',
#'/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_76X_V1p1_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160401_100956/0000/B2GEDMNtuple_1.root',
            ]))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
