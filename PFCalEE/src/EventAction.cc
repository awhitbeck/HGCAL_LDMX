#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"
#include "DetectorConstruction.hh"

#include "HGCSSInfo.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "TVector3.h"
#include "TMath.h"
//
EventAction::EventAction() {
	runAct = (RunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
	eventMessenger = new EventActionMessenger(this);
	printModulo = 100;
	outF_ = TFile::Open("PFcal.root", "RECREATE");
	outF_->cd();

	double xysize =
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->GetCalorSizeXY();

	//save some info
	HGCSSInfo *info = new HGCSSInfo();
	info->calorSizeXY(xysize);
	info->cellSize(CELL_SIZE_X);
	info->model(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getModel());
	info->version(
			((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->getVersion());

	std::cout << " -- check Info: version = " << info->version() << " model = "
			<< info->model() << std::endl;
	outF_->WriteObjectAny(info, "HGCSSInfo", "Info");
	initLayer = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->initLayer();
	steelThick = ((DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction())->steelThickness();
	//honeycomb
	geomConv_ = new HGCSSGeometryConversion(info->model(), CELL_SIZE_X);
	geomConv_->initialiseHoneyComb(xysize, CELL_SIZE_X);
	//square map for BHCAL
	geomConv_->initialiseSquareMap(xysize, 10.);

	tree_ = new TTree("HGCSSTree", "HGC Standalone simulation tree");
	tree_->Branch("HGCSSEvent", "HGCSSEvent", &event_);
	tree_->Branch("HGCSSSamplingSectionVec",
			"std::vector<HGCSSSamplingSection>", &ssvec_);
	tree_->Branch("HGCSSSimHitVec", "std::vector<HGCSSSimHit>", &hitvec_);
	//Branch containing incident particle and secondaries from the initial W target
	tree_->Branch("HGCSSTargetVec", "std::vector<HGCSSGenParticle>",
			&targetvec_);
	//Branch containing (''long lasting'') hadronic tracks
	tree_->Branch("HGCSSHadronVec", "std::vector<HGCSSGenParticle>", &hadronvec_);
	// }
}

//
EventAction::~EventAction() {
	outF_->cd();
	tree_->Write();
	outF_->Close();
	//fout_.close();
	delete eventMessenger;
}

//
void EventAction::BeginOfEventAction(const G4Event* evt) {
	evtNb_ = evt->GetEventID();
	if (evtNb_ % printModulo == 0) {
		G4cout << "\n---> Begin of event: " << evtNb_ << G4endl;
		CLHEP::HepRandom::showEngineStatus();
	}
	//fout_ << "Event " << evtNb_ << std::endl;

}

//
void EventAction::Detect(G4double eng, G4double edep, G4double stepl,
		G4double globalTime, G4int pdgId, G4VPhysicalVolume *volume,
		const G4ThreeVector & position, G4int trackID, G4int parentID,
		const HGCSSGenParticle & genPart, G4bool targetParticle) {
	G4bool inc_ = genPart.isIncoming();
	for (size_t i = 0; i < detector_->size(); i++)
		(*detector_)[i].add(eng, edep, stepl, globalTime, pdgId, volume,
				position, trackID, parentID, i,inc_);

	if (inc_){
		if (targetParticle){
			targetvec_.push_back(genPart);
		}
		else{
			hadronvec_.push_back(genPart);
		}
	}
}

//
void EventAction::EndOfEventAction(const G4Event* g4evt) {
	//return;
	bool debug(evtNb_ % printModulo == 0);
	hitvec_.clear();

	event_.eventNumber(evtNb_);
	event_.vtx_x(g4evt->GetPrimaryVertex(0)->GetX0());
	event_.vtx_y(g4evt->GetPrimaryVertex(0)->GetY0());
	event_.vtx_z(g4evt->GetPrimaryVertex(0)->GetZ0());
	event_.steelThick(steelThick);

	ssvec_.clear();
	ssvec_.reserve(detector_->size());
	//Changing initLayer because initial layers contain tracking sections.
	for (size_t i = initLayer; i < detector_->size(); i++) {
		HGCSSSamplingSection lSec;
		lSec.volNb(i);

		lSec.volX0trans((*detector_)[i].getAbsorberX0());
		lSec.voldEdx((*detector_)[i].getAbsorberdEdx());
		lSec.volLambdatrans((*detector_)[i].getAbsorberLambda());

		//Measurements now follow

		lSec.absorberDep((*detector_)[i].getAbsorbedEnergy());
		lSec.sensDep((*detector_)[i].getMeasuredEnergy(false));
		lSec.totalDep((*detector_)[i].getTotalEnergy());

		lSec.gDepFrac((*detector_)[i].getPhotonFraction());
		lSec.eDepFrac((*detector_)[i].getElectronFraction());
		lSec.muDepFrac((*detector_)[i].getMuonFraction());
		lSec.neutronDepFrac((*detector_)[i].getNeutronFraction());
		lSec.hadDepFrac((*detector_)[i].getHadronicFraction());

		lSec.avgTime((*detector_)[i].getAverageTime());
		lSec.nSiHits((*detector_)[i].getTotalSensHits());

		//Gen level quantities now follow

		lSec.neutronKinFlux((*detector_)[i].getKinNeutron());
		lSec.hadKinFlux((*detector_)[i].getKinHadron());
		lSec.muKinFlux((*detector_)[i].getKinMuon());

		lSec.neutronCount((*detector_)[i].getNeutronCount());
		lSec.hadCount((*detector_)[i].getHadronCount());
		lSec.muCount((*detector_)[i].getMuonCount());

		if (evtNb_ == 1)
			std::cout << "if (layer==" << i << ") return " << lSec.voldEdx()
					<< ";" << std::endl;

		TVector3 eleCentroid_, hadCentroid_, neutCentroid_, muCentroid_;
		double eleShowerSizeE_ = 0, muShowerSizeE_ = 0, neutShowerSizeE_ = 0, hadShowerSizeE_ = 0, eleEWgt_ = 0,
				hadEWgt_ = 0, neutEWgt_ = 0, muEWgt_ = 0;
		bool is_scint = (*detector_)[i].hasScintillator;
		for (unsigned idx(0); idx < (*detector_)[i].n_sens_elements; ++idx) {
			std::map<unsigned, HGCSSSimHit> lHitMap;
			std::pair<std::map<unsigned, HGCSSSimHit>::iterator, bool> isInserted;

			for (unsigned iSiHit(0);
					iSiHit < (*detector_)[i].getSiHitVec(idx).size();
					++iSiHit) {
				G4SiHit lSiHit = (*detector_)[i].getSiHitVec(idx)[iSiHit];
				HGCSSSimHit lHit(lSiHit, idx,
						is_scint ?
								geomConv_->squareMap() :
								geomConv_->hexagonMap());

				isInserted = lHitMap.insert(
						std::pair<unsigned, HGCSSSimHit>(lHit.cellid(), lHit));
				if (!isInserted.second)
					isInserted.first->second.Add(lSiHit);
				if (idx == 0) {
					TVector3 v1(lSiHit.hit_x, lSiHit.hit_y,0);
					Int_t pdgId_ = lSiHit.pdgId;
					Double_t parentEng = lSiHit.parentKE;
					if ((abs(pdgId_) == 11) || (abs(pdgId_) == 22)) {
						eleCentroid_ += v1 * parentEng;
						eleEWgt_ += parentEng;
					} else if (abs(pdgId_) == 13) {
						muCentroid_ += v1 * parentEng;
						muEWgt_ += parentEng;

					} else if (pdgId_ == 2112) {
						neutCentroid_ += v1 * parentEng;
						neutEWgt_ += parentEng;
					} else if ((abs(pdgId_) != 111) && (abs(pdgId_) != 310)
							&& (pdgId_ != -2212)) {
						hadCentroid_ += v1 * parentEng;
						hadEWgt_ += parentEng;
					}
				}
			}

			std::map<unsigned, HGCSSSimHit>::iterator lIter = lHitMap.begin();
			hitvec_.reserve(hitvec_.size() + lHitMap.size());
			for (; lIter != lHitMap.end(); ++lIter) {
				(lIter->second).calculateTime();
				hitvec_.push_back(lIter->second);
			}

		} //loop on sensitive layers

		if ((*detector_)[i].n_sens_elements > 0) {

			unsigned absSize_ = (*detector_)[i].getSiHitVec(0).size();

			for (unsigned jAbsHit(0); jAbsHit < absSize_; ++jAbsHit) {
				G4SiHit lAbsHit = (*detector_)[i].getSiHitVec(0)[jAbsHit];
				TVector3 v1(lAbsHit.hit_x, lAbsHit.hit_y, 0);
				Int_t pdgId_ = lAbsHit.pdgId;
				Double_t parentEng = lAbsHit.parentKE;
				if (jAbsHit == 0) {
					if (eleEWgt_ > 0)
						eleCentroid_ = eleCentroid_ * (1.0 / eleEWgt_);
					if (muEWgt_ > 0)
						muCentroid_ = muCentroid_ * (1.0 / muEWgt_);
					if (neutEWgt_ > 0)
						neutCentroid_ = neutCentroid_ * (1.0 / neutEWgt_);
					if (hadEWgt_ > 0)
						hadCentroid_ = hadCentroid_ * (1.0 / hadEWgt_);
				}
				if ((abs(pdgId_) == 11) || (abs(pdgId_) == 22)) {
					eleShowerSizeE_ += TMath::Power((v1 - eleCentroid_).Mag() * parentEng, 2);
				} else if (abs(pdgId_) == 13) {
					muShowerSizeE_ += TMath::Power((v1 - muCentroid_).Mag() * parentEng, 2);
				} else if (pdgId_ == 2112) {
					neutShowerSizeE_ += TMath::Power((v1 - neutCentroid_).Mag() * parentEng,
							2);
				} else if ((abs(pdgId_) != 111) && (abs(pdgId_) != 310)
						&& (pdgId_ != -2212)) {
					hadShowerSizeE_ += TMath::Power((v1 - hadCentroid_).Mag() * parentEng, 2);
				}
			}
			(eleEWgt_ > 0) ?
					lSec.eleShowerSize(TMath::Sqrt(eleShowerSizeE_) / eleEWgt_) :
					lSec.eleShowerSize(0);
			(muEWgt_ > 0) ?
					lSec.muShowerSize(TMath::Sqrt(muShowerSizeE_) / muEWgt_) :
					lSec.muShowerSize(0);
			(hadEWgt_ > 0) ?
					lSec.hadronShowerSize(TMath::Sqrt(hadShowerSizeE_) / hadEWgt_) :
					lSec.hadronShowerSize(0);
			(neutEWgt_ > 0) ?
					lSec.neutronShowerSize(TMath::Sqrt(neutShowerSizeE_) / neutEWgt_) :
					lSec.neutronShowerSize(0);

		}
		ssvec_.push_back(lSec);

		if (debug) {
			(*detector_)[i].report((i == 1));
		}

		(*detector_)[i].resetCounters();
	}
	if (debug) {
		G4cout << " -- Number of truth particles = " << targetvec_.size() << G4endl<< " -- Number of simhits = " << hitvec_.size() << G4endl
		<< " -- Number of sampling sections = " << ssvec_.size() << G4endl;

	}

	tree_->Fill();

	//reset vectors
	targetvec_.clear();
	hitvec_.clear();
	ssvec_.clear();
	hadronvec_.clear();
	trackids.clear();
}
