// -*- C++ -*-
//
// Package:    HiggsAnalysis/HbbTreeMaker
// Class:      HbbTreeMaker
// 
/**\class HbbTreeMaker HbbTreeMaker.cc HiggsAnalysis/HbbTreeMaker/src/HbbTreeMaker.cc

 Description: Produce ntuple for the MSSM Hbb search.

*/
//
// Original Author:  Matthias Schroeder
//         Created:  Fri, 18 Jul 2014 09:00:37 GMT
//
//


// system include files
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TVector2.h"

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "HiggsAnalysis/HbbTreeMaker/interface/HbbTreeMaker.h"


//
// constructors and destructor
//
HbbTreeMaker::HbbTreeMaker(const edm::ParameterSet& iConfig) {
  vertexCollectionTag_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
  genParticleCollectionTag_ = iConfig.getParameter<edm::InputTag>("GenParticleCollection");
  patJetCollectionTag_ = iConfig.getParameter<edm::InputTag>("PFJetCollection");
  genJetCollectionTag_ = iConfig.getParameter<edm::InputTag>("GenJetCollection");

  maxNumPartonsFromHiggs_ = 5;
  maxNumGenJets_ = 50;
  maxNumRecoJets_ = 50;

  // xsec is only retrieved once and the same value written to all events
  xsec_ = 0.;
  xsecIsStored_ = false;
}


HbbTreeMaker::~HbbTreeMaker() {
  delete tree_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
HbbTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // gen particles (associated to gen jets + higgs?)
  // gen jets
  // reco jets

  // Set default values for all branch variables
  setBranchVariablesToDefault();
  
  // Cross section (only read once from the event)
  if( !xsecIsStored_ ) {
    edm::Handle<GenRunInfoProduct> genRunInfoProduct;
    iEvent.getRun().getByLabel("generator",genRunInfoProduct);
    if( genRunInfoProduct.isValid() ) {
      xsec_ = genRunInfoProduct->crossSection();
    }
    xsecIsStored_ = true;
  }

  // Event information
  edm::EventAuxiliary aux = iEvent.eventAuxiliary();
  runNum_       = aux.run();
  lumiBlockNum_ = aux.luminosityBlock();
  evtNum_       = aux.event();
  
  // Number of vertices
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexCollectionTag_,vertices);
  if( vertices.isValid() ) {
    nVtx_ = vertices->size();
  }

  // Parton-level information of Higgs and its decay products
  fillHiggsAndPartonsFromHiggsInfo(iEvent);

  // Jet information
  fillGenJetsInfo(iEvent);
  fillRecoJetsInfo(iEvent);

  // Store event in tree
  tree_->Fill();
}



// store reco-jet info
void
HbbTreeMaker::fillRecoJetsInfo(const edm::Event& iEvent) {
  
  edm::Handle< edm::View<pat::Jet> > recoJets; // need to use PAT jets to get b-tag info easily
  iEvent.getByLabel(patJetCollectionTag_,recoJets);
  if( recoJets.isValid() ) {

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(genParticleCollectionTag_,genParticles);
    const bool hasGenParticles = genParticles.isValid();

    numRecoJets_ = 0;
    for(size_t iJet = 0; iJet < recoJets->size() && numRecoJets_ < maxNumRecoJets_; ++iJet) {
      // use separate counter in case of cuts
      
      // kinematic properties
      recoJetsPt_[numRecoJets_] = recoJets->at(iJet).pt();
      recoJetsEta_[numRecoJets_] = recoJets->at(iJet).eta();
      recoJetsPhi_[numRecoJets_] = recoJets->at(iJet).phi();
      recoJetsEnergy_[numRecoJets_] = recoJets->at(iJet).energy();

      // b-tagging information
      const std::vector< std::pair< std::string, float > > 
	discriminatorLabelValuePairs = recoJets->at(iJet).getPairDiscri();
      for(size_t iPair = 0; iPair < discriminatorLabelValuePairs.size(); ++iPair) {
	const std::string label = discriminatorLabelValuePairs.at(iPair).first;
	const float value = discriminatorLabelValuePairs.at(iPair).second;
	// if( iJet == 0 ) {
	//   std::cout << ">> " << label << " : " << value << std::endl;
	// }
	if( label == "combinedSecondaryVertexBJetTags" ) {
	  combSVBJetTags_[numRecoJets_] = value;
	}
      }	// end of loop over discriminatorLabelValuePairs

      // parton matching
      if( hasGenParticles ) {
	fillPartonsInJet(recoJets->at(iJet),numRecoJets_,genParticles,0.3);
      }

	// increment jet counter
      ++numRecoJets_;
    } // end of loop over reco jets
  }
}



// store gen-jet info
void
HbbTreeMaker::fillGenJetsInfo(const edm::Event& iEvent) {

  edm::Handle< edm::View<reco::Candidate> > genJets;
  iEvent.getByLabel(genJetCollectionTag_,genJets);
  if( genJets.isValid() ) {
    numGenJets_ = 0;
    for(size_t j = 0; j < genJets->size(); ++j) {
      genJetsPt_[numGenJets_] = genJets->at(j).pt();
      genJetsEta_[numGenJets_] = genJets->at(j).eta();
      genJetsPhi_[numGenJets_] = genJets->at(j).phi();
      genJetsEnergy_[numGenJets_] = genJets->at(j).energy();
      ++numGenJets_;
      if( numGenJets_ == maxNumGenJets_ ) break;
    }
  }
}



// store parton-level info of Higgs and its decay products
void
HbbTreeMaker::fillHiggsAndPartonsFromHiggsInfo(const edm::Event& iEvent) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleCollectionTag_,genParticles);
  if( genParticles.isValid() ) {
    for(size_t iP = 0; iP < genParticles->size(); ++iP) {
      const reco::GenParticleRef genParticle(genParticles,iP);
      int genPartPdg = genParticle->pdgId();
      if( genPartPdg==25 || genPartPdg==35 || genPartPdg==36 ) { // one of the Higgs bosons 
	// h (pdg=25) , H (pdg=35) or A (pdg=36)
	higgsPdg_    = genPartPdg;
	higgsEta_    = genParticle->eta();
	higgsPhi_    = genParticle->phi();
	higgsPt_     = genParticle->pt();
	higgsEnergy_ = genParticle->energy();
	higgsMass_   = genParticle->mass();

	// accessing daughters --->
	// this is a trick to get partons/quarks from Higgs decay
	bool partonsFromHiggsNotYetFound = true;
	bool isFirstIteration = true;
	const reco::Candidate* intermediateDaughter = NULL;
	numPartonsFromHiggs_ = 0; // number of partons is initially set to 0
	while( partonsFromHiggsNotYetFound ) {
                        
	  int numberOfDaughters = 0;
                        
	  if( isFirstIteration ) numberOfDaughters = genParticle->numberOfDaughters();
	  else                   numberOfDaughters = intermediateDaughter->numberOfDaughters();
                        
	  if(        numberOfDaughters == 0 ) { // at end of decay chain
	    partonsFromHiggsNotYetFound = false; // exit condition

	  } else if( numberOfDaughters == 1 ) { // some intermediate decay step
	    const reco::Candidate* daughter = NULL;
	    if( isFirstIteration ) daughter = genParticle->daughter(0);
	    else                   daughter = intermediateDaughter->daughter(0);
	    intermediateDaughter = daughter;
	    isFirstIteration = false;

	  } else {		// these are the final partons from the Higgs decay
	    for(int iDaughter=0; iDaughter<numberOfDaughters; ++iDaughter) {
	      const reco::Candidate* daughter = NULL;
	      if( isFirstIteration ) daughter = genParticle->daughter(iDaughter);
	      else                   daughter = intermediateDaughter->daughter(iDaughter);

	      const int pdgId = daughter->pdgId();
	      if( (pdgId != higgsPdg_ ) && ( numPartonsFromHiggs_ < maxNumPartonsFromHiggs_ ) ) {
		partonsFromHiggsPdg_[numPartonsFromHiggs_]    = daughter->pdgId();
		partonsFromHiggsEta_[numPartonsFromHiggs_]    = daughter->eta();
		partonsFromHiggsPhi_[numPartonsFromHiggs_]    = daughter->phi();
		partonsFromHiggsPt_[numPartonsFromHiggs_]     = daughter->pt();		   
		partonsFromHiggsEnergy_[numPartonsFromHiggs_] = daughter->energy();
		numPartonsFromHiggs_++;
	      }
	    }
	    partonsFromHiggsNotYetFound = false;
	  }
                        
	} // end of while( partonsFromHiggsNotYetFound )

	break;		// found the Higgs, leave the loop
	 
      } // end if is one of the Higgs bosons
    } // end of loop over gen paticles
	
  } // end valid GenParticleCollection

}


// store the number of UDSG, C, B partons inside deltaRMax
// of jet. Only status-2 partons, i.e. directly before the
// hadronisation step (which have string or cluster daughters)
// are considered.
void
HbbTreeMaker::fillPartonsInJet(const reco::Candidate& jet, const int recoJetIdx, const edm::Handle<reco::GenParticleCollection>& genParticles, const double deltaRMax) {
  for(size_t iP = 0; iP < genParticles->size(); ++iP) {
    const reco::GenParticleRef particle(genParticles,iP);
    const int partPdg = std::abs(particle->pdgId());
 
    //FIXME: Should we require here isthep == 2 or isthep == 3?
    //But the latter will probably work only for Pythia not for e.g. Herwig ...

    // select only u,d,s,c,b, and g
    if( (partPdg >= 1 && partPdg <= 5) || partPdg == 21 ) {
      const double dR = deltaR(*particle,jet);
      if( dR > deltaRMax ) continue;

      bool daughter_is_string = false;
      bool daughter_is_cluster = false;
      for(size_t iD = 0; iD < particle->numberOfDaughters(); ++iD) {
        const reco::Candidate* const daughter = particle->daughter(iD);
        const int pdgDaughter = std::abs(daughter->pdgId());
        if(      pdgDaughter == 92 ) daughter_is_string = true;
        else if( pdgDaughter == 91 ) daughter_is_cluster = true;
      }
      if( daughter_is_cluster || daughter_is_string ) {
	// this it the end of the parton showering because

	if(      partPdg == 5 ) numBInRecoJet_[recoJetIdx]++;
	else if( partPdg == 4 ) numCInRecoJet_[recoJetIdx]++;
	else                    numUDSGInRecoJet_[recoJetIdx]++;
      }
    } // end is u,d,s,c,b, or g
  } // end of loop over gen particles
}



void 
HbbTreeMaker::setBranchVariablesToDefault() {
  runNum_ = 0;      
  lumiBlockNum_ = 0;
  evtNum_ = 0;      
  nVtx_ = 0;

  higgsPdg_ = 0;   
  higgsEta_ = -99999.;   
  higgsPhi_ = -99999.;   
  higgsPt_ = -99999.;  
  higgsEnergy_ = -99999.;
  higgsMass_ = -99999.;

  numPartonsFromHiggs_ = 0;  
  for(int i = 0; i < maxNumPartonsFromHiggs_; ++i) {
    partonsFromHiggsPdg_[i] = 0;   
    partonsFromHiggsPt_[i] = -99999.;    
    partonsFromHiggsEta_[i] = -99999.;   
    partonsFromHiggsPhi_[i] = -99999.;   
    partonsFromHiggsEnergy_[i] = -99999.;
  }

  numGenJets_ = 0;
  for(int i = 0; i < maxNumRecoJets_; ++i) {
    genJetsPt_[i] = -99999.;    
    genJetsEta_[i] = -99999.;   
    genJetsPhi_[i] = -99999.;   
    genJetsEnergy_[i] = -99999.;
  }

  numRecoJets_ = 0;
  for(int i = 0; i < maxNumRecoJets_; ++i) {
    recoJetsPt_[i] = -99999.;    
    recoJetsEta_[i] = -99999.;   
    recoJetsPhi_[i] = -99999.;   
    recoJetsEnergy_[i] = -99999.;
    combSVBJetTags_[i] = -99999.;
    numUDSGInRecoJet_[i] = 0;
    numCInRecoJet_[i] = 0;
    numBInRecoJet_[i] = 0;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
HbbTreeMaker::beginJob() {
  edm::Service<TFileService> fs;
  if( !fs ) {
    throw edm::Exception(edm::errors::Configuration,
			 "TFile Service is not registered in cfg file");
  }
  tree_ = fs->make<TTree>("MSSMHbbTree","Ntuple for MSSM Hbb Analysis");
  // tree_->SetAutoSave(10000000000);
  // tree_->SetAutoFlush(1000000);

  tree_->Branch("XSec",&xsec_,"XSec/F");
  tree_->Branch("RunNum",&runNum_,"RunNum/i");
  tree_->Branch("LumiBlockNum",&lumiBlockNum_,"LumiBlockNum/i");
  tree_->Branch("EvtNum",&evtNum_,"EvtNum/i");
  tree_->Branch("NVtx",&nVtx_,"NVtx/s");

  // Generator-level information on Higgs boson kinematics and its decay products
  tree_->Branch( "HiggsPdg",    &higgsPdg_,    "HiggsPdg/I"    );
  tree_->Branch( "HiggsEta",    &higgsEta_,    "HiggsEta/F"    );
  tree_->Branch( "HiggsPhi",    &higgsPhi_,    "HiggsPhi/F"    );
  tree_->Branch( "HiggsPt",     &higgsPt_,     "HiggsPt/F"     );
  tree_->Branch( "HiggsEnergy", &higgsEnergy_, "HiggsEnergy/F" );
  tree_->Branch( "HiggsMass",   &higgsMass_,   "HiggsMass/F"   );
  tree_->Branch( "NumPartonsFromHiggs",   &numPartonsFromHiggs_,   "NumPartonsFromHiggs/I"                     );
  tree_->Branch( "PartonsFromHiggsPdg",    partonsFromHiggsPdg_,    "PartonsFromHiggsPdg[NumPartonsFromHiggs]/I"    );
  tree_->Branch( "PartonsFromHiggsPt",     partonsFromHiggsPt_,     "PartonsFromHiggsPt[NumPartonsFromHiggs]/F"     );
  tree_->Branch( "PartonsFromHiggsEta",    partonsFromHiggsEta_,    "PartonsFromHiggsEta[NumPartonsFromHiggs]/F"    );
  tree_->Branch( "PartonsFromHiggsPhi",    partonsFromHiggsPhi_,    "PartonsFromHiggsPhi[NumPartonsFromHiggs]/F"    );
  tree_->Branch( "PartonsFromHiggsEnergy", partonsFromHiggsEnergy_, "PartonsFromHiggsEnergy[NumPartonsFromHiggs]/F" );

  // Gen-jet info
  tree_->Branch( "NumGenJets",   &numGenJets_,    "NumGenJets/I"                );
  tree_->Branch( "GenJetsPt",     genJetsPt_,     "GenJetsPt[NumGenJets]/F"     );
  tree_->Branch( "GenJetsEta",    genJetsEta_,    "GenJetsEta[NumGenJets]/F"    );
  tree_->Branch( "GenJetsPhi",    genJetsPhi_,    "GenJetsPhi[NumGenJets]/F"    );
  tree_->Branch( "GenJetsEnergy", genJetsEnergy_, "GenJetsEnergy[NumGenJets]/F" );

  // Reco-jet info
  tree_->Branch( "NumRecoJets",   &numRecoJets_,    "NumRecoJets/I"                 );
  tree_->Branch( "RecoJetsPt",     recoJetsPt_,     "RecoJetsPt[NumRecoJets]/F"     );
  tree_->Branch( "RecoJetsEta",    recoJetsEta_,    "RecoJetsEta[NumRecoJets]/F"    );
  tree_->Branch( "RecoJetsPhi",    recoJetsPhi_,    "RecoJetsPhi[NumRecoJets]/F"    );
  tree_->Branch( "RecoJetsEnergy", recoJetsEnergy_, "RecoJetsEnergy[NumRecoJets]/F" );
  tree_->Branch( "CombSVBJetTags", combSVBJetTags_, "CombSVBJetTags[NumRecoJets]/F" );
  tree_->Branch( "NumUDSGInRecoJet", numUDSGInRecoJet_, "NumUDSGInRecoJet[NumRecoJets]/I" );
  tree_->Branch( "NumCInRecoJet",    numCInRecoJet_,    "NumCInRecoJet[NumRecoJets]/I"    );
  tree_->Branch( "NumBInRecoJet",    numBInRecoJet_,    "NumBInRecoJet[NumRecoJets]/I"    );
}


float HbbTreeMaker::deltaPhi(const reco::Candidate& cand1, const reco::Candidate& cand2) const {
  return TVector2::Phi_mpi_pi(cand1.phi()-cand2.phi());
}


float HbbTreeMaker::deltaR(const reco::Candidate& cand1, const reco::Candidate& cand2) const {
  const float dPhi = deltaPhi(cand1,cand2);
  const float dEta = cand1.eta()-cand2.eta();

  return sqrt( dPhi*dPhi + dEta*dEta );
}


// ------------ method called once each job just after ending the event loop  ------------
void 
HbbTreeMaker::endJob() {
  tree_->Write();
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HbbTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HbbTreeMaker);
