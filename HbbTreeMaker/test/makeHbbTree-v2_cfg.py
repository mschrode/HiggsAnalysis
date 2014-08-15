import FWCore.ParameterSet.Config as cms

## --- The process -----------------------------------------------------
from PhysicsTools.PatAlgos.patTemplate_cfg import *


## --- Log output ------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(
    placeholder = cms.untracked.bool(True)
)
process.MessageLogger.cout = cms.untracked.PSet(
    INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(1))
)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
        
        
## --- Files to process ------------------------------------------------
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
#        "/store/mc/Spring14dr/QCD_Pt-470to600_TuneZ2star_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/22306849-03C8-E311-991F-003048679150.root",
#        "/store/mc/Spring14dr/QCD_Pt-470to600_TuneZ2star_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/2291A129-29C8-E311-BB67-003048D15DB6.root",
#        "/store/mc/Spring14dr/QCD_Pt-470to600_TuneZ2star_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/26B139EC-09C8-E311-B9FD-003048FFCBFC.root",

#        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/084DD595-47CC-E311-8D97-848F69FD291F.root",
#        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/20AD10CC-03CC-E311-8051-00266CFAE7C4.root",

        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/084DD595-47CC-E311-8D97-848F69FD291F.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/0CB20439-93CC-E311-ADE0-7845C4FC346A.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/20AD10CC-03CC-E311-8051-00266CFAE7C4.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/241263C7-60CC-E311-9760-00266CF9B684.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/4C4F3F73-5ACC-E311-9E39-848F69FD29CA.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/5CF1C4A4-B4CC-E311-B79E-00A0D1EEE660.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/800657F5-07CC-E311-8F41-7845C4FC3635.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/8E687BFD-5BCC-E311-94CF-F04DA275C2FE.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/9091EEDF-03CC-E311-B08A-00266CFAE8C4.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/942086EB-5ECC-E311-8EA3-7845C4FC39D1.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/96D5BF7A-57CC-E311-81AE-00A0D1EEE68C.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/C4CC9A6A-03CC-E311-B8AE-7845C4FC39AD.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E6928625-65CC-E311-9555-7845C4FC37B5.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/E8AFA640-74CC-E311-ACB7-7845C4FC346A.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/ECAE813E-43CC-E311-B5AC-00266CFAE228.root",
        "/store/mc/Spring14dr/SUSYBBHToBB_M-500_13TeV-pythia6-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/F49E2E14-50CC-E311-A751-008CFA001DB8.root",
    )
)


## --- Conditions ------------------------------------------------------
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("POSTLS170_V7::All") # https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?redirectedfrom=CMS.SWGuideFrontierConditions#CSA14_exercise


## --- Output file -----------------------------------------------------
process.outpath.remove(process.out) # we do not want to write a PAT tuple
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("HbbTree.root")
)


## --- Primary vertices --------------------------------------------------------
process.goodPrimaryVertices = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2"),
    filter = cms.bool(True),    # rejects event if no vertex passing cut
)


## --- PAT configuration -------------------------------------------------------

# verbose flags for the PF2PAT modules
process.options.allowUnscheduled = cms.untracked.bool(True)

# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
jetAlgo="AK5"
usePF2PAT(
    process,
    runPF2PAT=True,
    jetAlgo=jetAlgo,
    runOnMC=True,
    postfix=postfix,
    jetCorrections=('AK5PFchs',['L1FastJet', 'L2Relative', 'L3Absolute'],'None'),
    pvCollection=cms.InputTag('goodPrimaryVertices'),
)

process.objects = cms.Sequence(    
    process.goodPrimaryVertices *
    getattr(process,"selectedPatJets"+postfix)
)

# this doesn't work yet
applyPostfix(process,"patJets",postfix).tagInfoSources = cms.VInputTag(
    "impactParameterTagInfosAODPFlow",
    "secondaryVertexTagInfosAODPFlow",
    "softMuonTagInfosAODPFlow"
)



## --- Configure the TreeMaker -----------------------------------------
from HiggsAnalysis.HbbTreeMaker.HbbTreeMaker_cfi import HbbTreeMaker
process.HbbTree = HbbTreeMaker.clone(
    VertexCollection = cms.InputTag("goodPrimaryVertices"),
    GenParticleCollection = cms.InputTag("genParticles"),
    PFJetCollection = cms.InputTag("selectedPatJets"+postfix),
    GenJetCollection = cms.InputTag("ak5GenJets"),
)


process.dump = cms.EDAnalyzer("EventContentAnalyzer")
process.hbbTreeMaker = cms.Path(
    process.objects *
    #    process.dump *
    process.HbbTree
)
