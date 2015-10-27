from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'tWantitop' 
config.General.workArea = 'MetNTuples_tWantitop'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonPogCrabNtuplesMC_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 100000
config.Data.publication = True
config.Data.publishDataName = 'MuonPOG_tWantitop'
config.Data.outLFNDirBase = '/store/user/federica/NTuplesMuonPOG/MC' 
config.Data.ignoreLocality = True

config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T2_CN_Beijing']
