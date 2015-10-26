from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ZMuMu_M50_120_STARTUP' 
config.General.workArea = 'NTuples_ZMuMu_M50_120_STARTUP'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonPogCrabNtuplesStartupMC_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIISpring15DR74-Startup_EXOReReco_74X_Spring15_mcRun2_startup_v0-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 200000
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
#config.Data.runRange = '203777-208686'

config.Data.publication = True
config.Data.publishDataName = 'ZMuMu_M50_120_STARTUP' 
config.Data.outLFNDirBase = '/store/user/federica/NTuplesMuonPOG/MC' 

config.Site.storageSite = 'T2_CH_CERN'

