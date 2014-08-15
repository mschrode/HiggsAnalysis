import FWCore.ParameterSet.Config as cms

HbbTreeMaker = cms.EDAnalyzer(
    "HbbTreeMaker",

    VertexCollection = cms.InputTag("goodVertices"),                       
    GenParticleCollection = cms.InputTag("genParticles"),
    PFJetCollection = cms.InputTag("ak5PFJetsCHS"),
    GenJetCollection = cms.InputTag("ak5GenJets"),
)
