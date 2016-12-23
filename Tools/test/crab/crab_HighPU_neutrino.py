from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HighPU_Neutrino'

config.section_('JobType')
config.JobType.pluginName  = 'ANALYSIS'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6',
                              'ntupleName=muonPOGNtuple_highPU_Neutrino.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all'
               ]

config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
#config.Data.inputDataset = '/ZeroBiasBunchTrains5/Run2016H-09Nov2016-v1/AOD'

config.Data.inputDataset = '/SingleNeutrino/RunIIHighPUTrainsDR-HighPUTrains_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'

#config.Data.runRange = '273299'
#config.Data.allowNonValidInputDataset = True

config.Data.splitting    = 'LumiBased'
config.Data.unitsPerJob  = 150  # Since files based, 10 files per job
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outLFNDirBase  = '/store/group/phys_muon/calderon/NTuplesMuonPOG/'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

