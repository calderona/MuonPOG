import FWCore.ParameterSet.Config as cms

process = cms.Process("PRINT")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
        "file:///afs/cern.ch/work/e/escalant/public/EvDisplaysPromptReco/EvDisplaysPromptRecoMu.root"
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "74X_dataRun2_Prompt_v0"

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

#from SLHCUpgradeSimulations.Configuration.postLS1Customs import *
#process = customise_HLT( process )

process.muonEventDumper = cms.EDAnalyzer("MuonEventDumper",
                             TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                             TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                             #TrigResultsTag = cms.untracked.InputTag("none"),
                             #TrigSummaryTag = cms.untracked.InputTag("none"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             GenTag = cms.untracked.InputTag("none"),
                             PileUpInfoTag = cms.untracked.InputTag("none")
                             )

process.AOutput = cms.EndPath(process.muonEventDumper)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
