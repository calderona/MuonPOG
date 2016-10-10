from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

runRange = '273017-273450'

config.General.requestName = 'SingleMuon_Prompt_v2_' + runRange 
config.General.workArea = 'NTuples_Prompt16'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['minMuPt=10.','minNMu=2','globalTag=80X_dataRun2_Prompt_v8','runOnMC=False','ntupleName=./muonPOGNtuple.root']

config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100

#Information on the most up to date JSON is available here: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.runRange = runRange

#config.Data.outLFNDirBase = '/store/user/battilan/MuonPOG/Ntuples/'
config.Data.outLFNDirBase = '/store/group/phys_muon/llinwei/Ntuples/'
config.Data.publication = False
config.Data.outputDatasetTag = 'SingleMuon_Prompt_v2_' + runRange

config.Site.storageSite = 'T2_CH_CERN' 

