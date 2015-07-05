import FWCore.ParameterSet.Config as cms


def appendMuonPogNtuple(process, runOnMC, processTag="HLT", ntupleFileName="MuonPogTree.root") :

    process.load("MuonPOG.Tools.MuonPogTreeProducer_cfi")

    if processTag != "HLT" :
        print "[MuonPogNtuples]: Customising process tag for TriggerResults / Summary to :", processTag
        process.MuonHltTree.TrigResultsTag = "TriggerResults::"+processTag
        process.MuonHltTree.TrigSummaryTag = "hltTriggerSummaryAOD::"+processTag

    if runOnMC :
        process.load("MuonPOG.Tools.PrunedGenParticles_cfi")
        process.muonHltNtuple = cms.Sequence(process.prunedGenParticles + process.MuonPogTree)
    else :
        process.muonHltNtuple = cms.Sequence(process.MuonPogTree)

    process.TFileService = cms.Service('TFileService',
        fileName = cms.string(ntupleFileName)
    )

    if hasattr(process,"AOutput") :
        print "[MuonPogNtuples]: EndPath AOutput found, appending ntuples"
        process.AOutput.replace(process.hltOutputA, process.hltOutputA + process.muonPogNtuple)
    else :
        print "[MuonPogNtuples]: EndPath AOuptput not found, creating it for ntuple sequence"
        process.AOutput = cms.EndPath(process.muonPogNtuple)
