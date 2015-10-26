from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'WJets' 
config.General.workArea = 'MetNTuples_WJets'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonPogCrabNtuplesMC_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'


config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 1000000
config.Data.publication = True
config.Data.publishDataName = 'MuonPOG_WJets'
config.Data.outLFNDirBase = '/store/user/federica/NTuplesMuonPOG/MC' 

config.Site.storageSite = 'T2_CH_CERN'

