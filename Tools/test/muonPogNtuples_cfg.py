import FWCore.ParameterSet.Config as cms

import subprocess

runOnMC = False
#pathCut   = "all"
#filterCut = "all"
pathCut   = "HLT_IsoMu20_v"
filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))



process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()

)

if runOnMC :
    process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9A')
    sourcefilesfolder = " /store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/60000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]    
else :
    process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v4')
    sourcefilesfolder = "/store/data/Run2015D/SingleMuon/AOD/PromptReco-v3/000/258/158/00000"
    files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", sourcefilesfolder ])
    process.source.fileNames = [ sourcefilesfolder+"/"+f for f in files.split() ]

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from MuonPOG.Tools.MuonPogNtuples_cff import appendMuonPogNtuple

if runOnMC :
    ntupleName = "ntuples_DY_NLO.root"
else :
    ntupleName = "ntuples_SingleMu.root"
    
appendMuonPogNtuple(process,runOnMC,"HLT",ntupleName,pathCut,filterCut)
