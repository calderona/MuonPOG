from CRABClient.UserUtilities import config
config = config()

config.General.requestName  = 'SingleMu2015D_25ns_Dataset_reReco'
config.General.workArea     = 'NTuples_SingleMu2015D_25ns_reReco'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName   = '../muonPogCrabNtuplesData2015D_reReco_cfg.py'   
config.JobType.priority   = 1

config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/AOD'
config.Data.inputDBS    = 'global'
config.Data.splitting   = 'LumiBased' 
config.Data.unitsPerJob = 100

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'

config.Data.outputDatasetTag = 'SingleMu2015D_25ns_reReco' 
config.Data.allowNonValidInputDataset = True

config.Data.outLFNDirBase    = '/store/user/battilan/NTuplesMuonPOG/reRECO/Data/' 
config.Data.ignoreLocality   = True


config.Site.storageSite = 'T3_IT_Bologna'

