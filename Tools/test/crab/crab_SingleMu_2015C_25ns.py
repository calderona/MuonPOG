from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMu2015C_25ns-AllLumiS' 
config.General.workArea = 'NTuples_SingleMu2015C_25ns'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonPogCrabNtuplesData2015C_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset =  '/SingleMuon/Run2015C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased' 
config.Data.unitsPerJob = 100

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt'

#config.Data.runRange = '203777-208686'

config.Data.publication = True
config.Data.publishDataName = 'SingleMu2015C_25ns-AllLumiS' 
config.Data.outLFNDirBase = '/store/user/federica/NTuplesMuonPOG/Data' 

config.Site.storageSite = 'T2_CH_CERN'

