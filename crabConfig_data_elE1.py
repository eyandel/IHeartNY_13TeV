from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Data_elE1'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_dataEF_el.sh'

config.section_("Data")
config.Data.inputDataset = '/SingleElectron/vorobiev-SingleElectron_Run2016E-23Sep2016-v1_B2GAnaFW_80X_v2p4-961c7d882d8721e72fac616aaa90ecc1/USER'
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
