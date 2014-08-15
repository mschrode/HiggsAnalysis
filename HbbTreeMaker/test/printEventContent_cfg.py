import FWCore.ParameterSet.Config as cms

process = cms.Process("EventContentPrinter")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                '/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/084DD595-47CC-E311-8D97-848F69FD291F.root'
                            )
                    )

process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(process.dump)
