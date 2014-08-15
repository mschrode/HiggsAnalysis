#ifndef HIGGS_ANALYSIS_HBB_TREE_MAKER_H
#define HIGGS_ANALYSIS_HBB_TREE_MAKER_H

// -*- C++ -*-
//
// Package:    HiggsAnalysis/HbbTreeMaker
// Class:      HbbTreeMaker
// 
/**\class HbbTreeMaker HiggsAnalysis/HbbTreeMaker/sr/HbbTreeMaker.cc

 Description: Produce ntuple for the MSSM Hbb search.

*/
//
// Original Author:  Matthias Schroeder
//         Created:  Fri, 18 Jul 2014 09:00:37 GMT
//
//


// system include files
#include <memory>

// user include files
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TTree.h"


class HbbTreeMaker : public edm::EDAnalyzer {
public:
  explicit HbbTreeMaker(const edm::ParameterSet&);
  ~HbbTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  TTree* tree_;

  // Event information 
  Float_t xsec_;
  bool xsecIsStored_;
  UInt_t runNum_;      
  UInt_t lumiBlockNum_;
  UInt_t evtNum_;
  edm::InputTag vertexCollectionTag_;
  UShort_t nVtx_;

  // Gen-particle information
  edm::InputTag genParticleCollectionTag_;

  // Higgs boson
  Int_t   higgsPdg_;   
  Float_t higgsEta_;   
  Float_t higgsPhi_;   
  Float_t higgsPt_;    
  Float_t higgsEnergy_;
  Float_t higgsMass_;

  // Decay products of Higgs
  Int_t   maxNumPartonsFromHiggs_;
  Int_t   numPartonsFromHiggs_;  
  Int_t   partonsFromHiggsPdg_[5];   
  Float_t partonsFromHiggsPt_[5];    
  Float_t partonsFromHiggsEta_[5];   
  Float_t partonsFromHiggsPhi_[5];   
  Float_t partonsFromHiggsEnergy_[5];

  // Gen jets
  edm::InputTag genJetCollectionTag_;
  Int_t   maxNumGenJets_;
  Int_t   numGenJets_;  
  Float_t genJetsPt_[50];    
  Float_t genJetsEta_[50];   
  Float_t genJetsPhi_[50];   
  Float_t genJetsEnergy_[50];

  // Reco jets
  edm::InputTag patJetCollectionTag_;
  Int_t   maxNumRecoJets_;
  Int_t   numRecoJets_;  
  // kinematic info
  Float_t recoJetsPt_[50];    
  Float_t recoJetsEta_[50];   
  Float_t recoJetsPhi_[50];   
  Float_t recoJetsEnergy_[50];
  // b-tagging info
  Float_t combSVBJetTags_[50];
  // parton info (number matched partons)
  Int_t numUDSGInRecoJet_[50];
  Int_t numCInRecoJet_[50];
  Int_t numBInRecoJet_[50];



  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void setBranchVariablesToDefault();
  void fillHiggsAndPartonsFromHiggsInfo(const edm::Event& iEvent);
  void fillGenJetsInfo(const edm::Event& iEvent);
  void fillRecoJetsInfo(const edm::Event& iEvent);
  void fillPartonsInJet(const reco::Candidate& jet, const int recoJetIdx, const edm::Handle<reco::GenParticleCollection>& genParticles, const double deltaRMax);
  float deltaPhi(const reco::Candidate& cand1, const reco::Candidate& cand2) const;
  float deltaR(const reco::Candidate& cand1, const reco::Candidate& cand2) const;
};

#endif
