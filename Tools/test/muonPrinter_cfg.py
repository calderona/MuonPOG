import FWCore.ParameterSet.Config as cms

process = cms.Process("PRINT")

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_7_6_0/RelValZMM_13/GEN-SIM-RECO/PU25ns_76X_mcRun2_asymptotic_v11-v1/00000/2E482387-FA7F-E511-9F8E-0025905A60A6.root"
             ),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = "76X_mcRun2_asymptotic_v11"

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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
