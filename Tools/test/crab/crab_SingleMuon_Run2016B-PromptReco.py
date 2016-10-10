from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'SingleMuonRun2016B_PromptReco_Run273299'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_Prompt_v9',
                              'ntupleName=muonPOGNtuple_SingleMuonRun2016B_PromptReco.root',
                              'nEvents=-1',
                              'runOnMC=False',
                              'hltPathFilter=all',
                              'minMuPt=0.0',
                              'minNMu=0'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/AOD'

config.Data.runRange = '273299'

config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 150  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/group/phys_muon/calderon/NTuplesMuonPOG/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

