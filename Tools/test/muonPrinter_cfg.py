import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("PRINT")

options = VarParsing.VarParsing()

options.register('globalTag',
                 '76X_dataRun2_v5', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('inputFile',
                 'file://./260373_249111510_MINIAOD_76X.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Input file (local or remote)")

options.register('runOnMiniAOD',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on AOD or miniAOD")

options.register('runOnMC',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on DATA or MC")

options.register('plotTriggerInfo',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Plot Trigger information (for AOD)")

options.register('adHocTrackCollection',
                 'none', #default value                                                                                                                              
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Print a specific track collection (e.g. sta not updated at vtx, or general tracks)")

options.parseArguments()


process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(options.inputFile),
        secondaryFileNames = cms.untracked.vstring()
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options.globalTag

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")

process.muonEventDumper = cms.EDAnalyzer("MuonEventDumper",
                             TrigResultsTag = cms.untracked.InputTag("none"),
                             TrigSummaryTag = cms.untracked.InputTag("none"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             GenTag = cms.untracked.InputTag("none"),
                             PileUpInfoTag = cms.untracked.InputTag("none"),

                             AdHocTrackTag = cms.untracked.InputTag(options.adHocTrackCollection)
                             )

if options.plotTriggerInfo and not options.runOnMiniAOD :
    process.muonEventDumper.TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT")
    process.muonEventDumper.TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT")

if options.runOnMiniAOD :
    process.muonEventDumper.MuonTag          = cms.untracked.InputTag("slimmedMuons")
    process.muonEventDumper.PrimaryVertexTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
    
if options.runOnMC :
    process.muonEventDumper.GenTag = cms.untracked.InputTag("genParticles")
    process.muonEventDumper.PileUpInfoTag = cms.untracked.InputTag("addPileupInfo")


process.AOutput = cms.EndPath(process.muonEventDumper)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
