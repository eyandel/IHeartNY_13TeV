from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'WJets_HT-100To200'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs', 'pileup_reweight_mu.root']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_MC.sh'

config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/skinnari-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1_B2GAnaFW_80X_V2p4-bfea8033ab2d179bbb8e0faf6e2dc0cf/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.totalUnits = 1000
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/skinnari/13TeV_full2016'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
