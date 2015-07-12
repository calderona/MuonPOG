import FWCore.ParameterSet.Config as cms


def appendMuonPogNtuple(process, runOnMC, processTag="HLT", ntupleFileName="MuonPogTree.root") :

    process.load("MuonPOG.Tools.MuonPogTreeProducer_cfi")

    if processTag != "HLT" :
        print "[MuonPogNtuples]: Customising process tag for TriggerResults / Summary to :", processTag
        process.MuonPogTree.TrigResultsTag = "TriggerResults::"+processTag
        process.MuonPogTree.TrigSummaryTag = "hltTriggerSummaryAOD::"+processTag

    if runOnMC :
        process.load("MuonPOG.Tools.PrunedGenParticles_cfi")
        process.muonPogNtuple = cms.Sequence(process.prunedGenParticles + process.MuonPogTree)
    else :
        process.muonPogNtuple = cms.Sequence(process.MuonPogTree)
        PileUpInfoTag = cms.untracked.InputTag("none")
        GenTag = cms.untracked.InputTag("none")
        
    process.TFileService = cms.Service('TFileService',
        fileName = cms.string(ntupleFileName)
    )

    if hasattr(process,"AOutput") :
        print "[MuonPogNtuples]: EndPath AOutput found, appending ntuples"
        process.AOutput.replace(process.hltOutputA, process.hltOutputA + process.muonPogNtuple)
    else :
        print "[MuonPogNtuples]: EndPath AOuptput not found, creating it for ntuple sequence"
        process.AOutput = cms.EndPath(process.muonPogNtuple)
