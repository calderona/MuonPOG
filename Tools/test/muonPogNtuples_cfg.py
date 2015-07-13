import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/244/00000/1A55C9D8-5327-E511-BA1A-02163E0133D1.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/244/00000/5C0D9766-6727-E511-A61F-02163E014729.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/244/00000/BE92D76E-6E27-E511-96D2-02163E0133A7.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/251/00000/3004473B-8B27-E511-AD12-02163E0120B3.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/252/00000/283B8530-8D27-E511-8310-02163E014289.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/252/00000/70B9C1B4-9227-E511-ACC2-02163E01413E.root',
            '/store/data/Run2015B/DoubleMuon/AOD/PromptReco-v1/000/251/252/00000/7859893B-9C27-E511-BE3A-02163E01414A.root'
        ),
        secondaryFileNames = cms.untracked.vstring()

)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "74X_dataRun2_Prompt_v0"

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple
appendMuonPogNtuple(process,False,"HLT","ntuple_DoubleMuon_251244_251252.root")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
