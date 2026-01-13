#ifndef egammaSeedingEfficiencyStudy_H
#define egammaSeedingEfficiencyStudy_H

/**\class egammaSeedingEfficiencyStudy egammaSeedingEfficiencyStudy.cc egammaGPUdevelopment/egammaSeedEfficiency/plugins/egammaSeedingEfficiencyStudy.cc
 Description: Macro that looks into the pixel tracker hit reconstruction efficiency of electrons. Goal is to
	  		  understand why we loose efficiency when doubles are removed from the electron seeding.
 Implementation:
    - Find the TrackingParticles that are electrons
	- Find the simHits & RecHits associated with these tracking particles
	- Look into the number & pattern of reconstructed pixel hits for each TP
*/

// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"

// CMSSW framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// CMSSW DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/TrackBase.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

// CMSSW sim DataFormats
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/UniqueSimTrackId.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other includes 
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"

#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

class egammaSeedingEfficiencyStudy : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
  public:
		explicit egammaSeedingEfficiencyStudy(const edm::ParameterSet&);
		~egammaSeedingEfficiencyStudy();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		virtual void initialize();

		reco::GenParticle get_lastcopy_prefsrSN(reco::GenParticle part);
		reco::GenParticle get_lastcopySN(reco::GenParticle part);
		reco::GenParticleCollection get_genpartsSN(reco::GenParticleCollection, int);
		std::pair<reco::GenParticle*, int> match_to_genSN(double, double, reco::GenParticleCollection, double);
		reco::GenParticle* matchGenSN(reco::Electron, const reco::GenParticleCollection&);

		const edm::EDGetTokenT<reco::ElectronCollection> electronToken;
		const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
		const edm::EDGetTokenT<reco::ElectronSeedCollection> electronSeedToken;
		// const edm::EDGetTokenT<std::vector<SimTrack>> simtracksToken;
		// const edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken;
		// const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
		// const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
		// std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHit_;
		// const edm::EDGetTokenT<ClusterTPAssociation>  clusterTPAssocToken_;
		// const edm::EDGetTokenT<SiPixelRecHitCollection>  pixelRecHitToken_;

		TTree* mctree_;
		TTree* prtree_;
		TTree* prseedtree_;

		double 	mcElePt_;
		double 	mcElePhi_;
		double 	mcEleEta_;
		double 	mcEleE_;
		ushort  nFound_;
		bool		isFound_;

		double 	prElePt_;
		double 	prElePhi_;
		double 	prEleEta_;
		double 	prEleE_;
		ushort	nMatched_;
		bool		isMatched_;
		double 	prEleSeedPt_;
		double 	prEleSeedPhi_;
		double 	prEleSeedEta_;
		double 	prEleSeedE_;
		unsigned int prEleSeedUniqueID_;

		double minDR_;

		int run_, lumi_, event_;
		// bool eventPassed_ = false;

		bool verbose_;
		double DeltaR_;

    // Define either Zee or ttbar for the process
    bool isZee_ = true;
    bool isttbar_ = false;
    uint8_t eMotherPDG_;
};

//Constructor
egammaSeedingEfficiencyStudy::egammaSeedingEfficiencyStudy(const edm::ParameterSet& iConfig): 
	electronToken     (consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electron"))),
	genParticlesToken (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
	electronSeedToken (consumes<reco::ElectronSeedCollection>(iConfig.getParameter<edm::InputTag>("electronSeed"))),
	// trackingParticlesToken (consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
	// topoToken_(esConsumes()),
	// geomToken_(esConsumes()),
	// simHit_(),
	// clusterTPAssocToken_(consumes<ClusterTPAssociation>(iConfig.getParameter<edm::InputTag>("cluster2TPSrc"))),
	// pixelRecHitToken_(consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("pixelRecHits"))),
	verbose_(iConfig.getParameter<bool>("verbose")),
	DeltaR_(iConfig.getParameter<double>("deltaR")),
	isZee_(iConfig.getParameter<bool>("isZee")),
	isttbar_(iConfig.getParameter<bool>("isttbar"))
{
	initialize();
	usesResource("TFileService");

	// std::vector<edm::InputTag> tags = iConfig.getParameter<std::vector<edm::InputTag>>("simHitSrc");
	// simHit_.reserve(tags.size());
	// for (auto const &tag : tags) 
	// {
	// 	simHit_.emplace_back(consumes<edm::PSimHitContainer>(tag));
	// }

  assert((isZee_ != isttbar_) && "Can't set both processes Z -> ee and ttbar or none of them, chose one.");
  // If we are working with Zee, then the mother PDG id is 23 (Z), for ttbar it is 24 (W^+-)
  eMotherPDG_ = isZee_ ? 23 : 24;
}

//Destructor
egammaSeedingEfficiencyStudy::~egammaSeedingEfficiencyStudy() {}

void egammaSeedingEfficiencyStudy::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void egammaSeedingEfficiencyStudy::initialize() {

	mcElePt_    = 0.;
	mcElePhi_   = 0.;
	mcEleEta_   = 0.;
	mcEleE_     = 0.;
	nFound_     = 0;
	isFound_    = false;

	prElePt_    = 0.;
	prElePhi_   = 0.;
	prEleEta_   = 0.;
	prEleE_     = 0.;
	nMatched_   = 0;
	isMatched_  = false;
	prEleSeedPt_    = 0.;
	prEleSeedPhi_   = 0.;
	prEleSeedEta_   = 0.;
	prEleSeedE_     = 0.;
	prEleSeedUniqueID_ = 0;

	minDR_ 			= std::numeric_limits<double>::max();

  // eventPassed_ = false;
	run_    = 0;
	lumi_   = 0;
	event_  = 0;
}


void egammaSeedingEfficiencyStudy::beginJob() 
{

	// Access the TFileService
	edm::Service<TFileService> fs; 

	// Create the TTree
	mctree_ = fs->make<TTree>("mctree"  , "mctree");

	mctree_->Branch("event",    &event_,  "event/I");
	mctree_->Branch("lumi",     &lumi_,   "lumi/I");
	mctree_->Branch("run",      &run_,    "run/I");
	// mctree_->Branch("passed",	&eventPassed_, 	"passed/O");

	mctree_->Branch("mcElePt",  &mcElePt_, 	"mcElePt/D");
	mctree_->Branch("mcElePhi", &mcElePhi_, "mcElePhi/D");
	mctree_->Branch("mcEleEta", &mcEleEta_, "mcEleEta/D");
	mctree_->Branch("mcEleE",   &mcEleE_, 	"mcEleE/D");
	mctree_->Branch("nFound",   &nFound_, 	"nFound/s");
	mctree_->Branch("isFound"   &isFound_, 	"isFound/O");

	mctree_->Branch("minDR",    &minDR_,    "minDR/D");
	mctree_->Branch("prElePt",  &prElePt_,  "prElePt/D");
	mctree_->Branch("prElePhi", &prElePhi_, "prElePhi/D");
	mctree_->Branch("prEleEta", &prEleEta_, "prEleEta/D");
	mctree_->Branch("prEleE",   &prEleE_, 	"prEleE/D");
	mctree_->Branch("prEleSeedPt", 	&prEleSeedPt_, 	"prEleSeedPt/D");
	mctree_->Branch("prEleSeedPhi",	&prEleSeedPhi_, "prEleSeedPhi/D");
	mctree_->Branch("prEleSeedEta",	&prEleSeedEta_, "prEleSeedEta/D");
	mctree_->Branch("prEleSeedE",   &prEleSeedE_,   "prEleSeedE/D");
	mctree_->Branch("prEleSeedUniqueID",  &prEleSeedUniqueID_,  "prEleSeedUniqueID/i");

	// Create the TTree
	prtree_ = fs->make<TTree>("prtree"  , "prtree");

	prtree_->Branch("event",    &event_,  "event/I");
	prtree_->Branch("lumi",     &lumi_,   "lumi/I");
	prtree_->Branch("run",      &run_,    "run/I");
	// prtree_->Branch("passed",	&eventPassed_, 	"passed/O");

	prtree_->Branch("mcElePt",  &mcElePt_, 	"mcElePt/D");
	prtree_->Branch("mcElePhi", &mcElePhi_, "mcElePhi/D");
	prtree_->Branch("mcEleEta", &mcEleEta_, "mcEleEta/D");
	prtree_->Branch("mcEleE",   &mcEleE_, 	"mcEleE/D");

	prtree_->Branch("minDR",    &minDR_,    "minDR/D");
	prtree_->Branch("prElePt",  &prElePt_,  "prElePt/D");
	prtree_->Branch("prElePhi", &prElePhi_, "prElePhi/D");
	prtree_->Branch("prEleEta", &prEleEta_, "prEleEta/D");
	prtree_->Branch("prEleE",   &prEleE_,   "prEleE/D");
	prtree_->Branch("nMatched", &nMatched_, "nMatched/s");
	prtree_->Branch("isMatched",&isMatched_,"isMatched/O");
	prtree_->Branch("prEleSeedPt", 	&prEleSeedPt_, 	"prEleSeedPt/D");
	prtree_->Branch("prEleSeedPhi",	&prEleSeedPhi_, "prEleSeedPhi/D");
	prtree_->Branch("prEleSeedEta",	&prEleSeedEta_, "prEleSeedEta/D");
	prtree_->Branch("prEleSeedE",   &prEleSeedE_,   "prEleSeedE/D");
	prtree_->Branch("prEleSeedUniqueID",  &prEleSeedUniqueID_,  "prEleSeedUniqueID/i");
}


void egammaSeedingEfficiencyStudy::endJob() {}
void egammaSeedingEfficiencyStudy::endRun(edm::Run const&, edm::EventSetup const&) {}

using P = std::pair<OmniClusterRef, TrackingParticleRef>;
bool compare(const P &i, const P &j) { return i.second.index() > j.second.index(); }

void egammaSeedingEfficiencyStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	using namespace edm;
	using namespace std;

	initialize();

	// const TrackerTopology* tTopo = &iSetup.getData(topoToken_);
	// const TrackerGeometry* tGeom = &iSetup.getData(geomToken_);

	edm::Handle<reco::ElectronCollection> electronH;
	iEvent.getByToken(electronToken, electronH);
	edm::Handle<reco::GenParticleCollection> genParticlesH;
	iEvent.getByToken(genParticlesToken, genParticlesH);
	edm::Handle<reco::ElectronSeedCollection> electronSeedH;
	iEvent.getByToken(electronSeedToken, electronSeedH);

	//-------------- Event Info -----------------------------------
	run_    = iEvent.id().run();
	event_  = iEvent.id().event();
	lumi_   = iEvent.id().luminosityBlock();

	//-------------- Gen particle info -----------------------------------
	// https://cmssdt.cern.ch/dxr/CMSSW/source/DataFormats/HepMCCandidate/interface/GenParticle.h

	std::vector<reco::GenParticle> genElectrons;

	for (auto genItr = genParticlesH->begin(); genItr != genParticlesH->end(); ++genItr) {
		const auto& genPart = *genItr;
		if (abs(genItr->pdgId())==11) {
			if (std::abs(genPart.mother()->pdgId()) == eMotherPDG_ && genPart.isPromptFinalState())	{
				if (genPart.status()==1)
					genElectrons.push_back(genPart);
				else if (genPart.status()==2)
					genElectrons.push_back(get_lastcopy_prefsrSN(genPart));
				else if (genPart.status()==3)
					genElectrons.push_back(get_lastcopySN(genPart));
			}
			else {
				if (genPart.numberOfMothers() == 0)
					genElectrons.push_back(genPart);
			}
		}
	}

	if (verbose_) {
		std::cout<<" --- Collecion of gen electrons ---- "<<std::endl;
		for(size_t i = 0; i != genElectrons.size(); ++i) {
			std::cout<<" Gen electron "<< i <<" 4-momentum :("<< genElectrons.at(i).pt() <<","<<genElectrons.at(i).eta()<<","<<genElectrons.at(i).phi() <<","<< genElectrons.at(i).energy()<<")"<<std::endl;	
		}
		std::cout<<" ------------------------------------"<<std::endl;
	}

	// edm::Handle<TrackingParticleCollection> TrackingParticleH;
	// iEvent.getByToken(trackingParticlesToken, TrackingParticleH);
	// const TrackingParticleCollection& trackingParticles = *TrackingParticleH;
	
	auto sq = [](float x) { return x * x; };
	auto dr = [sq](float x1,float x2,float y1,float y2) {return std::sqrt(sq(x1 - x2) + sq(y1 - y2));};
	
	//-------------- hltGsfElectrons -----------------------------------~
	// DR matching with gen electron identified previously
	// https://cmssdt.cern.ch/dxr/CMSSW/source/DataFormats/EgammaCandidates/interface/GsfElectron.h

	if (electronH.isValid()) {

		// For each MC electron check whether it has a reco counterpart. If so, consider it found.
		// If no reco counterpart is found, the electron is considered missing
		for (const auto& genEle : genElectrons) {

			mcElePt_ 	= genEle.pt();
			mcEleEta_ = genEle.eta();
			mcElePhi_ = genEle.phi();
			mcEleE_   = genEle.energy();
			nFound_		= 0;
			isFound_  = false;

			minDR_ = std::numeric_limits<double>::max();

			if (electronH->size() == 0) {
				prElePt_  = -1.;
				prEleEta_ = -9.;
				prElePhi_ = -4.;
				prEleE_   = -1.;
				prEleSeedUniqueID_  = -1;
				minDR_    = -1.;
				mctree_->Fill();
				continue;
			}

			for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) {

				const float deltaR = dr(genEle.eta(), eleItr->eta(), genEle.phi(), eleItr->phi());
			
				isFound_ |= deltaR < DeltaR_;
				// There can be multiple reco electrons close to the MC electron. Only store the information of
				// the closest reco electrons, but count the reco electrons that would be considered found.
				if (deltaR < DeltaR_ and deltaR < minDR_) {
					nFound_++;
					prElePt_  = eleItr->pt();
					prEleEta_ = eleItr->eta();
					prElePhi_ = eleItr->phi();
					prEleE_   = eleItr->energy();
					const auto& seed = eleItr->gsfTrack()->seedRef();
					prEleSeedUniqueID_  = seed->uniqueID();
					minDR_    = deltaR;
				}
			}

			// No matching reco electron found, find the closest PR electron anyway
			if (not isFound_) {

				minDR_ = std::numeric_limits<double>::max();

				for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) {

					const float deltaR = dr(genEle.eta(), eleItr->eta(), genEle.phi(), eleItr->phi());
				
					if (deltaR < minDR_) {
						prElePt_  = eleItr->pt();
						prEleEta_ = eleItr->eta();
						prElePhi_ = eleItr->phi();
						prEleE_   = eleItr->energy();
						const auto& seed = eleItr->gsfTrack()->seedRef();
						prEleSeedUniqueID_ = seed->uniqueID();
						minDR_    = deltaR;
					}
				}
			}
			mctree_->Fill();
		}

		// For each reco electron check whether it has an MC counterpart. If so, consider it matched.
		// If no MC counterpart is found, the reco electron is considered a fake.
		for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) {
			const auto& seed = eleItr->gsfTrack()->seedRef();

			prElePt_    = eleItr->pt();
			prEleEta_   = eleItr->eta();
			prElePhi_   = eleItr->phi();
			prEleE_     = eleItr->energy();
			prEleSeedUniqueID_ = seed->uniqueID();
			nMatched_   = 0;
			isMatched_  = false;

			minDR_ = std::numeric_limits<double>::max();

			if (genElectrons.size() == 0) {
				mcElePt_  = -1.;
				mcEleEta_ = -9.;
				mcElePhi_ = -4.;
				mcEleE_   = -1.;
				minDR_    = -1.;
				prtree_->Fill();
				continue;
			}

			for(const auto& genEle : genElectrons) {
			
				const float deltaR = dr(genEle.eta(), eleItr->eta(), genEle.phi(), eleItr->phi());
			
				isMatched_ |= deltaR < DeltaR_;
				// There can be multiple MC electrons close to the reco electron. Only store the information of
				// the closest MC electrons, but count the MC electrons that would be considered matched.
				if (deltaR < DeltaR_ and deltaR < minDR_) {
					nMatched_++;
					mcElePt_ 	= genEle.pt();
					mcEleEta_ = genEle.eta();
					mcElePhi_ = genEle.phi();
					mcEleE_   = genEle.energy();
					minDR_ 		= deltaR;
				}
			}

			// No matching MC electron found, find the closest MC electron anyway
			if (not isMatched_) {

				minDR_ = std::numeric_limits<double>::max();

				for(const auto& genEle : genElectrons) {
				
					const float deltaR = dr(genEle.eta(), eleItr->eta(), genEle.phi(), eleItr->phi());
				
					if (deltaR < DeltaR_ and deltaR < minDR_) {
						mcElePt_ 	= genEle.pt();
						mcEleEta_ = genEle.eta();
						mcElePhi_ = genEle.phi();
						mcEleE_   = genEle.energy();
						minDR_ 		= deltaR;
					}
				}
			}
			prtree_->Fill();
		}

	}
}

// Taken from here : https://github.com/waredjeb/ElectronPixelMatching/blob/add_gen_particles/SeedFromTrack/SeedFromTrackAnalyzer/plugins/ElectronMatchSeed.cc#L206

reco::GenParticle egammaSeedingEfficiencyStudy::get_lastcopy_prefsrSN(reco::GenParticle part) {
	auto daughters = part.daughterRefVector();
	if (daughters.size() == 1 && daughters.at(0)->pdgId() == part.pdgId()) {
		return get_lastcopy_prefsrSN(*(daughters.at(0)));
	}
	return part;
}

reco::GenParticle egammaSeedingEfficiencyStudy::get_lastcopySN(reco::GenParticle part) {
	auto daughters = part.daughterRefVector();
	for (size_t p = 0; p != daughters.size(); ++p) {
		if (daughters.at(p)->pdgId() == part.pdgId()) {
			return get_lastcopySN(*(daughters.at(p)));
		}
	}
	return part;
}

reco::GenParticleCollection egammaSeedingEfficiencyStudy::get_genpartsSN(reco::GenParticleCollection genparts, int status = 2) {
	std::vector<reco::GenParticle> selected;
	for (size_t i = 0; i != genparts.size(); ++i) {
		auto part = genparts.at(i);
		if (abs(part.pdgId()) == 11 and fabs(part.mother()->pdgId()) == eMotherPDG_) {
			if (part.isHardProcess()) {
				if (status == 1) {
					selected.push_back(part);
				} else if (status == 2) {
					selected.push_back(get_lastcopy_prefsrSN(part));
				} else if (status == 3) {
					selected.push_back(get_lastcopySN(part));
				} else {
					throw std::runtime_error("error status not implemented");
				}
			}
		} else {
			if (part.numberOfMothers() == 0) {
				selected.push_back(part);
			}
		}

	}

	return selected;
}


std::pair<reco::GenParticle*, int> egammaSeedingEfficiencyStudy::match_to_genSN(double eta_reco, double phi_reco, reco::GenParticleCollection genParticles, double max_dr = 0.1)
{
	reco::GenParticle* best_match = nullptr;
	double best_dr2 = max_dr * max_dr;
	auto selected_parts = get_genpartsSN(genParticles);
	//std::cout << "Number selected particles: " << selected_parts.size() << std::endl;
	for (size_t i = 0; i != selected_parts.size(); ++i) {
		auto s_part = selected_parts.at(i);

		auto sq = [](float x) { return x * x; };
		auto dr2 = sq(eta_reco - s_part.eta()) + sq(phi_reco - s_part.phi());
		//std::cout<<" dr2 : "<< dr2<<std::endl;
		if (dr2 < best_dr2) {
			best_match = &(selected_parts.at(i));
			best_dr2 = dr2;
		}
	}
	// return best_match
	return std::make_pair(best_match, selected_parts.size());
}

reco::GenParticle* egammaSeedingEfficiencyStudy::matchGenSN(reco::Electron electron, const reco::GenParticleCollection& genParticle ) {
	auto eta = electron.eta();
	auto phi = electron.phi();

	auto matchedEle = match_to_genSN(eta,phi, genParticle).first;

	return matchedEle;
}


void egammaSeedingEfficiencyStudy::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {}
void egammaSeedingEfficiencyStudy::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void egammaSeedingEfficiencyStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(egammaSeedingEfficiencyStudy);

#endif 
