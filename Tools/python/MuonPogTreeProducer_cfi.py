import FWCore.ParameterSet.Config as cms

MuonPogTree = cms.EDAnalyzer("MuonPogTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                             TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             PFMetTag         = cms.untracked.InputTag("pfMet"), 
                             PFChMetTag       = cms.untracked.InputTag("pfChMet"), 
                             CaloMetTag       = cms.untracked.InputTag("caloMet"), 

                             GenTag = cms.untracked.InputTag("prunedGenParticles"), # pruned
                             PileUpInfoTag = cms.untracked.InputTag("addPileupInfo"),
                             GenInfoTag = cms.untracked.InputTag("generator")
                             )



