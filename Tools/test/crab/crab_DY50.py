from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DY50_reReco' 
config.General.workArea = 'NTuples_DY50_reReco'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../muonPogCrabNtuplesMC_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset =  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' 
config.Data.unitsPerJob = 100000

config.Data.outLFNDirBase    = '/store/user/battilan/NTuplesMuonPOG/reRECO/MC/DY' 
config.Data.ignoreLocality   = True


config.Site.storageSite = 'T3_IT_Bologna'

