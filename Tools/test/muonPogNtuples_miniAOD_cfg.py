import FWCore.ParameterSet.Config as cms

import subprocess

runOn761 = False
pathCut   = "all"
filterCut = "all"
#pathCut   = "HLT_IsoMu20_v"
#filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()

)

if runOn761 :
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v15')
    sourcefilesfolder = "/store/relval/CMSSW_7_6_1/DoubleMuon/MINIAOD/76X_dataRun2_v15_rerecoGT_RelVal_doubMu2015D-v1/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]    
else :
    process.GlobalTag.globaltag = cms.string('76X_dataRun2_v10')
    sourcefilesfolder = "/store/relval/CMSSW_7_6_0/DoubleMuon/MINIAOD/76X_dataRun2_v10_RelVal_doubMu2015D-v1/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple

if runOn761 :
    ntupleName = "ntuples_761.root"
else :
    ntupleName = "ntuples_760.root"
    
appendMuonPogNtuple(process,False,"HLT",ntupleName,pathCut,filterCut)

process.MuonPogTree.MuonTag = cms.untracked.InputTag("slimmedMuons")
process.MuonPogTree.PrimaryVertexTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
process.MuonPogTree.TrigResultsTag = cms.untracked.InputTag("none")
process.MuonPogTree.TrigSummaryTag = cms.untracked.InputTag("none")
process.MuonPogTree.PFMetTag = cms.untracked.InputTag("none")
process.MuonPogTree.PFChMetTag = cms.untracked.InputTag("none")
process.MuonPogTree.CaloMetTag = cms.untracked.InputTag("none")


