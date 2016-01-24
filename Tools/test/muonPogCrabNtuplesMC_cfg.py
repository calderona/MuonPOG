import FWCore.ParameterSet.Config as cms

import subprocess

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

pathCut   = "HLT_IsoMu20_v"
filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
ntupleName = "MuonPOGntuples.root"
runOnMC = True     

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.source = cms.Source("PoolSource",                           
                            fileNames = cms.untracked.vstring('file:PlaceHolder.root')
                            )

process.GlobalTag.globaltag = cms.string('76X_mcRun2_asymptotic_v12')

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple

appendMuonPogNtuple(process,runOnMC,"HLT",ntupleName,pathCut,filterCut)
