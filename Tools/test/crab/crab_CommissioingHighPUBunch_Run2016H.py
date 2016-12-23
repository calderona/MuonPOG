from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'CommissioningHighPU_Run2016H'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_SPECIALHIGHPUFILL_v0',
                              'ntupleName=muonPOGNtuple_commissioningHighPU_Run2016H.root',
                              'nEvents=-1',
                              'runOnMC=False',
                              'hltPathFilter=all'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
#config.Data.inputDataset = '/ZeroBiasBunchTrains5/Run2016H-09Nov2016-v1/AOD'

config.Data.inputDataset = '/Commissioning/Run2016H-09Nov2016-v1/AOD'

#config.Data.runRange = '273299'
config.Data.allowNonValidInputDataset = True

config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 150  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/group/phys_muon/calderon/NTuplesMuonPOG/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

