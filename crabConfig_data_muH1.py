from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Data_muH1'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_dataH_mu.sh'

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/vorobiev-Run2016H-PromptReco-v1_B2GAnaFW_80X_v2p4-376a23645e94877b22a7f32873431514/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 1000
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/dittmer/13TeV_full2016'
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.whitelist = ["T2_HU_Budapest"]
