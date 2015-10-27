from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'SingleMu2015D_25ns-LumiS_256629-256843' 
config.General.requestName = 'SingleMu2015D_25ns_Dataset-v3_AllLumiS'
config.General.workArea = 'NTuples_SingleMu2015D_25ns_v3'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'muonPogCrabNtuplesData2015D_cfg.py'   
config.JobType.priority = 1

config.Data.inputDataset = '/SingleMuon/Run2015D-PromptReco-v3/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased' 
config.Data.unitsPerJob = 100

#old json 19 Sept
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#new json 12 Oct
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt'

#config.Data.runRange = '203777-208686'

config.Data.publication = True
config.Data.publishDataName = 'SingleMu2015D_25ns_Dataset-v3_AllLumiS' 
config.Data.outLFNDirBase = '/store/user/federica/NTuplesMuonPOG/Data' 
config.Data.ignoreLocality = True


config.Site.whitelist = ["T2_IT_Bari"]
config.Site.storageSite = 'T2_CH_CERN'

