from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Data_2015C_el'
config.General.workArea = 'test'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.inputFiles = ['FrameworkJobReport.xml', 'execute_for_crab.py', 'iheartny_topxs_fwlite.py', 'JECs']
config.JobType.outputFiles = ['test_iheartNY.root']
config.JobType.scriptExe = 'execute_iheartNY_data.sh'

config.section_("Data")
config.Data.inputDataset = '/SingleElectron/jkarancs-B2GAnaFW_76X_V1p1_Run2015C_25ns-16Dec2015-v1-69b00753dd36562e8813bc06510c861e/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 1000
config.Data.ignoreLocality = True
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/dittmer/13TeV'
# This string is used to construct the output dataset name

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.whitelist = ["T2_HU_Budapest"]
