from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DYJets_Spring16'
config.General.workArea = 'NTuples_Spring16'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../muonPogNtuples_cfg.py'   
config.JobType.pyCfgParams = ['hltPathFilter=IsoMu20']
config.JobType.priority = 1

config.Data.inputDataset =  '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/AODSIM'
config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = 10


config.Data.outLFNDirBase    = '/store/group/phys_muon/Ntuples/2016'


config.Site.storageSite = 'T2_CH_CERN'
