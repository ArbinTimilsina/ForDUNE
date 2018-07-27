////////////////////////////////////////////////////////////////////////
// Class:        MatchingWithSingleCRT
// Module Type:  analyzer
// File:         MatchingWithSingleCRT_module.cc
// Author:       Arbin Timilsina, arbint@bnl.gov
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// C++ includes
#include <cmath>
#include <random>

// MuonCounter Geometry includes
#include "dune/Geometry/MuonCounter35Alg.h"

// Header file
#include "MatchingWithSingleCRT_module.h"

#include "../MCTruthInformation.h"

using namespace std;


namespace CRTMatching
{
    /////////////////////////////////////////////////////////
    //Constructor
    /////////////////////////////////////////////////////////
    MatchingWithSingleCRT::MatchingWithSingleCRT(fhicl::ParameterSet const& parameterSet)
	: EDAnalyzer(parameterSet)
    {
	// Read in the parameters from the .fcl file
	this->reconfigure(parameterSet);
    }


    /////////////////////////////////////////////////////////
    //Destructor
    /////////////////////////////////////////////////////////
    MatchingWithSingleCRT::~MatchingWithSingleCRT() {}

    /////////////////////////////////////////////////////////
    //Reads parameters form the .fcl file
    /////////////////////////////////////////////////////////
    void MatchingWithSingleCRT::reconfigure(fhicl::ParameterSet const& p)
    {

	fTrackModuleLabel = p.get< std::string >("TrackModuleLabel");
	fFlashModuleLabel = p.get< std::string >("FlashModuleLabel");

	fSelectedPDG = p.get< int >("PDGcode");
	fCRTsConfiguration = p.get< int >("CRTsConfiguration");

	return;
    }

    /////////////////////////////////////////////////////////
    //Executes once at the beginning of the job
    /////////////////////////////////////////////////////////
    void MatchingWithSingleCRT::beginJob()
    {
	nEvents = 0;
	nHitsCRT = 0, nHitsCRT_F = 0, nHitsCRT_B = 0;
	nTotalRecoTracks = 0;

	for (unsigned int c = 0; c <= 3; c++)
	    {
		nPrimaryMuons[c] = 0;
		nPrimaryMuonsWithOneHit[c] = 0;
		nPrimaryMuonsWithTwoHits[c] = 0;

		nConsideredRecoTracks[c] = 0;

		nPrimaryMatchedRecoTracksOneHit[c] = 0;
		nAllCRTMatchedRecoTracksOneHit[c] = 0;
		nGoodCRTMatchedRecoTracksOneHit[c] = 0;

		nPrimaryMatchedRecoTracksTwoHits[c] = 0;
		nAllCRTMatchedRecoTracksTwoHits[c] = 0;
		nGoodCRTMatchedRecoTracksTwoHits[c] = 0;
	    }

	// Open a basic log file, will overwrite a pre-existing one
	logFile.open("MatchingWithSingleCRT.log");

	// Get local time
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	logFile << "MatchingWithSingleCRT_module log file, " << asctime(timeinfo) << endl;

	//Load the CRT positions from a text file
	char counterfile[500];
	if(fCRTsConfiguration == 0)
	    {
		strcpy(counterfile, "/nashome/a/arbint/DuneApp/SpaceChargeEffects/CRTMatching/CRTsDefault.txt");
		logFile << "Default CRT geometry loaded." << endl;
	    }
	else if(fCRTsConfiguration == 1)
	    {
		strcpy(counterfile, "/nashome/a/arbint/DuneApp/SpaceChargeEffects/CRTMatching/CRTsNew.txt");
		logFile << "New CRT geometry loaded." << endl;
	    }
	else
	    {
		logFile << "ERROT: No CRT geometry loaded." << endl;
	    }

	int countersLoaded = geoMuonCounter->loadMuonCounterGeometry(counterfile, counterGeometry);

	if(!countersLoaded)
	    {
		logFile << "ERROR: CRT geometry failed to load." << endl;
	    }
	else
	    {
		logFile << "CRT geometry loaded: " << counterfile << endl;
	    }

	// Access ART's TFileService, which will handle creating and writing histograms and n-tuples
	art::ServiceHandle<art::TFileService> fileServiceHandle;

	//Histograms
	hStatistics = fileServiceHandle->make<TH1D>("hStatistics", "Job Statistics", 100, 0, 100);

	hRecoTrackEnergy = fileServiceHandle->make<TH2D>("hRecoTrackEnergy", "", 400, 0, 20.0, 400, 0, 20.0);
	hRecoTrackLengthYZ = fileServiceHandle->make<TH1D>("hRecoTrackLengthYZ", "", 160, 0, 800.0);
	hRecoTrackLengthYZVsEnergy = fileServiceHandle->make<TH2D>("hRecoTrackLengthYZVsEnergy", "", 400, 0, 20.0, 160, 0, 800.0);

	hHitsVsEnergy = fileServiceHandle->make<TH2D>("hHitsVsEnergy", "", 400, 0, 20.0, 500, 0, 5000.0);

	hPE = fileServiceHandle->make<TH1D>("hPE", "", 100, 0.0, 500.0);
	hFlashPeakTime = fileServiceHandle->make<TH1D>("hFlashPeakTime", "", 840, -4200.0, 4200.0);

	hDeltaT = fileServiceHandle->make<TH1D>("hDeltaT", "", 200, -1.0, 1.0);
	hDeltaY = fileServiceHandle->make<TH1D>("hDeltaY", "", 200, 0.0, 200.0);
	hDeltaX = fileServiceHandle->make<TH1D>("hDeltaX", "", 200, 0.0, 200.0);

	for (unsigned int c = 0; c <= 3; c++)
	    {
		hRecoTracksCoverageYX[c] = fileServiceHandle->make<TH2D>(Form("hRecoTracksCoverageYX_%u", c), "", 770 / dY, -385.0, 385.0, 610 / dY, 0.0, 610.0);
		hRecoTracksCoverageYZ[c] = fileServiceHandle->make<TH2D>(Form("hRecoTracksCoverageYZ_%u", c), "", 700 / dY, -1.0, 699.0, 610 / dY, 0.0, 610.0);

		hPrimaryMuonsWithOneHit[c] = fileServiceHandle->make<TH1D>(Form("hPrimaryMuonsWithOneHit_%u", c), "", 20, 0, 20.0);
		hPrimaryMatchedRecoTracksOneHit[c] = fileServiceHandle->make<TH1D>(Form("hPrimaryMatchedRecoTrackOneHit_%u", c), "", 20, 0, 20.0);
		hAllCRTMatchedRecoTracksOneHit[c] = fileServiceHandle->make<TH1D>(Form("hAllCRTMatchedRecoTrackOneHit_%u", c), "", 20, 0, 20.0);
		hGoodCRTMatchedRecoTracksOneHit[c] = fileServiceHandle->make<TH1D>(Form("hGoodCRTMatchedRecoTrackOneHit_%u", c), "", 20, 0, 20.0);
		hCoverageYXOneHit[c] = fileServiceHandle->make<TH2D>(Form("hCoverageYXOneHit_%u", c), "", 770 / dY, -385.0, 385.0, 610 / dY, 0.0, 610.0);
		hCoverageYZOneHit[c] = fileServiceHandle->make<TH2D>(Form("hCoverageYZOneHit_%u", c), "", 700 / dY, -1.0, 699.0, 610 / dY, 0.0, 610.0);

		hPrimaryMuonsWithTwoHits[c] = fileServiceHandle->make<TH1D>(Form("hPrimaryMuonsWithTwoHits_%u", c), "", 20, 0, 20.0);
		hPrimaryMatchedRecoTracksTwoHits[c] = fileServiceHandle->make<TH1D>(Form("hPrimaryMatchedRecoTrackTwoHits_%u", c), "", 20, 0, 20.0);
		hAllCRTMatchedRecoTracksTwoHits[c] = fileServiceHandle->make<TH1D>(Form("hAllCRTMatchedRecoTrackTwoHits_%u", c), "", 20, 0, 20.0);
		hGoodCRTMatchedRecoTracksTwoHits[c] = fileServiceHandle->make<TH1D>(Form("hGoodCRTMatchedRecoTrackTwoHits_%u", c), "", 20, 0, 20.0);
		hCoverageYXTwoHits[c] = fileServiceHandle->make<TH2D>(Form("hCoverageYXTwoHits_%u", c), "", 770 / dY, -385.0, 385.0, 610 / dY, 0.0, 610.0);
		hCoverageYZTwoHits[c] = fileServiceHandle->make<TH2D>(Form("hCoverageYZTwoHits_%u", c), "", 700 / dY, -1.0, 699.0, 610 / dY, 0.0, 610.0);
	    }
    }

    /////////////////////////////////////////////////////////
    //Executes once at the end of the job
    /////////////////////////////////////////////////////////
    void MatchingWithSingleCRT::endJob()
    {
	hStatistics->SetBinContent(1, nEvents);

	hStatistics->SetBinContent(11, nPrimaryMuons[0]);
	hStatistics->SetBinContent(12, nPrimaryMuons[1]);
	hStatistics->SetBinContent(13, nPrimaryMuons[2]);
	hStatistics->SetBinContent(14, nPrimaryMuons[3]);

	hStatistics->SetBinContent(15, nPrimaryMuonsWithOneHit[0]);
	hStatistics->SetBinContent(16, nPrimaryMuonsWithOneHit[1]);
	hStatistics->SetBinContent(17, nPrimaryMuonsWithOneHit[2]);
	hStatistics->SetBinContent(18, nPrimaryMuonsWithOneHit[3]);
	hStatistics->SetBinContent(19, nPrimaryMuonsWithTwoHits[0]);
	hStatistics->SetBinContent(20, nPrimaryMuonsWithTwoHits[1]);
	hStatistics->SetBinContent(21, nPrimaryMuonsWithTwoHits[2]);
	hStatistics->SetBinContent(22, nPrimaryMuonsWithTwoHits[3]);

	hStatistics->SetBinContent(23, nHitsCRT);
	hStatistics->SetBinContent(24, nHitsCRT_F);
	hStatistics->SetBinContent(25, nHitsCRT_B);

	hStatistics->SetBinContent(26, nTotalRecoTracks);

	hStatistics->SetBinContent(31, nConsideredRecoTracks[0]);
	hStatistics->SetBinContent(32, nConsideredRecoTracks[1]);
	hStatistics->SetBinContent(33, nConsideredRecoTracks[2]);
	hStatistics->SetBinContent(34, nConsideredRecoTracks[3]);

	hStatistics->SetBinContent(41, nPrimaryMatchedRecoTracksOneHit[0]);
	hStatistics->SetBinContent(42, nPrimaryMatchedRecoTracksOneHit[1]);
	hStatistics->SetBinContent(43, nPrimaryMatchedRecoTracksOneHit[2]);
	hStatistics->SetBinContent(44, nPrimaryMatchedRecoTracksOneHit[3]);
	hStatistics->SetBinContent(45, nPrimaryMatchedRecoTracksTwoHits[0]);
	hStatistics->SetBinContent(46, nPrimaryMatchedRecoTracksTwoHits[1]);
	hStatistics->SetBinContent(47, nPrimaryMatchedRecoTracksTwoHits[2]);
	hStatistics->SetBinContent(48, nPrimaryMatchedRecoTracksTwoHits[3]);

	hStatistics->SetBinContent(51, nAllCRTMatchedRecoTracksOneHit[0]);
	hStatistics->SetBinContent(52, nAllCRTMatchedRecoTracksOneHit[1]);
	hStatistics->SetBinContent(53, nAllCRTMatchedRecoTracksOneHit[2]);
	hStatistics->SetBinContent(54, nAllCRTMatchedRecoTracksOneHit[3]);
	hStatistics->SetBinContent(55, nAllCRTMatchedRecoTracksTwoHits[0]);
	hStatistics->SetBinContent(56, nAllCRTMatchedRecoTracksTwoHits[1]);
	hStatistics->SetBinContent(57, nAllCRTMatchedRecoTracksTwoHits[2]);
	hStatistics->SetBinContent(58, nAllCRTMatchedRecoTracksTwoHits[3]);

	hStatistics->SetBinContent(61, nGoodCRTMatchedRecoTracksOneHit[0]);
	hStatistics->SetBinContent(62, nGoodCRTMatchedRecoTracksOneHit[1]);
	hStatistics->SetBinContent(63, nGoodCRTMatchedRecoTracksOneHit[2]);
	hStatistics->SetBinContent(64, nGoodCRTMatchedRecoTracksOneHit[3]);
	hStatistics->SetBinContent(65, nGoodCRTMatchedRecoTracksTwoHits[0]);
	hStatistics->SetBinContent(66, nGoodCRTMatchedRecoTracksTwoHits[1]);
	hStatistics->SetBinContent(67, nGoodCRTMatchedRecoTracksTwoHits[2]);
	hStatistics->SetBinContent(68, nGoodCRTMatchedRecoTracksTwoHits[3]);

	cout << endl << endl;
	cout << "Primary Muons: " << nPrimaryMuons[0] << endl;
	cout << "with one hit: " << nPrimaryMuonsWithOneHit[0] << endl;
	cout << "with two hits: " << nPrimaryMuonsWithTwoHits[0] << endl;
	cout << "Total CRT hits: " << nHitsCRT << ", hits front: " << nHitsCRT_F << ", hits back: " << nHitsCRT_B << endl;
	cout << "Track Reconstruction Efficiency (one hit): " << 100.0 * nPrimaryMatchedRecoTracksOneHit[0] / nPrimaryMuonsWithOneHit[0] << " %" << endl;
	cout << "Matching Purity (one hit): " << 100.0 * nGoodCRTMatchedRecoTracksOneHit[0] / nAllCRTMatchedRecoTracksOneHit[0] << " %" << endl;
	cout << "Matching Efficiency (one hit): " << 100.0 * nGoodCRTMatchedRecoTracksOneHit[0] / nPrimaryMuonsWithOneHit[0] << " %" << endl;
	cout << "Track Reconstruction Efficiency (two hits): " << 100.0 * nPrimaryMatchedRecoTracksTwoHits[0] / nPrimaryMuonsWithTwoHits[0] << " %" << endl;
	cout << "Matching Purity (two hits): " << 100.0 * nGoodCRTMatchedRecoTracksTwoHits[0] / nAllCRTMatchedRecoTracksTwoHits[0] << " %" << endl;
	cout << "Matching Efficiency (two hits): " << 100.0 * nGoodCRTMatchedRecoTracksTwoHits[0] / nPrimaryMuonsWithTwoHits[0] << " %" << endl;
	cout << endl << endl;

	logFile << "EndJob called, all done!!!" << endl;
    }

    /////////////////////////////////////////////////////////
    //Executes once per event
    /////////////////////////////////////////////////////////
    void MatchingWithSingleCRT::analyze( const art::Event & event )
    {
	nEvents++;

	// The basics
	int fEvent  = event.id().event();
	int fRun    = event.run();
	int fSubRun = event.subRun();
	logFile << "Event " << fEvent << ", run " << fRun << ", subrun " << fSubRun << endl;

	verbo = (fRun == 1);

	//Clear the vectors for each event
	hitsCRT.clear();

	//Detector clocks service
	auto const* detectorClocksService = lar::providerFrom<detinfo::DetectorClocksService>();

	//Detector properties service
	auto const* detectorPropertiesService = lar::providerFrom<detinfo::DetectorPropertiesService>();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Store primary hits
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(verbo)
	    {
		logFile << endl;
		logFile << "Beginning Primary and CRT hit info!!" << endl;
	    }

	//Map of trackID and energy to see how primary were matched to reco tracks
	map<int, double> primaryMatchedRecoTrackOneHit[4];
	map<int, double> primaryMatchedRecoTrackTwoHits[4];

	//Get the MC truth and CRT hit information
	int tempId = 0;
	std::vector< art::Handle< std::vector<simb::MCTruth> > > allMCTruthList;
	event.getManyByType(allMCTruthList);
	for(size_t mcl = 0; mcl < allMCTruthList.size(); ++mcl)
	    {
		art::Handle< std::vector<simb::MCTruth> > allMCTruthListHandle = allMCTruthList[mcl];
		for(size_t mch = 0; mch < allMCTruthListHandle->size(); ++mch)
		    {
			art::Ptr<simb::MCTruth> mcTruth(allMCTruthListHandle, mch);
			auto particleOrigin = mcTruth->Origin();

			for(Int_t iParticle = 0; iParticle < mcTruth->NParticles(); iParticle++)
			    {
				const simb::MCParticle& particle(mcTruth->GetParticle(iParticle));

				int primaryPDG = particle.PdgCode();
				double primaryEnergy = particle.E();
				double primaryT0 = particle.T();

				//Make sure it's a muon
				if (fabs(primaryPDG) != fabs(fSelectedPDG))
				    {
					continue;
				    }

				// Match this primary muon to a GEANT track and assign GEANT track id
				int primaryTrackId = -9999;
				const sim::ParticleList& geantList = particleInventory->ParticleList();
				for (const auto& PartPair : geantList)
				    {
					const simb::MCParticle& geantParticle = *(PartPair.second);
					if((primaryPDG == geantParticle.PdgCode()) &&
					   (fabs(particle.Px() - geantParticle.Px()) < 0.0001) &&
					   (fabs(particle.Py() - geantParticle.Py()) < 0.0001) &&
					   (fabs(particle.Pz() - geantParticle.Pz()) < 0.0001))
					    {
						primaryTrackId = geantParticle.TrackId();
						break;
					    }
				    }
				//Skip the particle if it can't be assigned to a GEANT track
				if(primaryTrackId == -9999)
				    {
					continue;
				    }

				// Check if the primary muon falls with the readout window
				const simb::MCParticle* geantParticle = particleInventory->TrackIdToParticle_P(primaryTrackId);
				bool inReadOutWindow = MCTruthInformation::InReadoutWindow(geantParticle);
				if (!inReadOutWindow)
				    {
					continue;
				    }

				//Find where this primary muon originated from
				string primaryOrigin;
				if (particleOrigin == simb::kCosmicRay)
				    {
					primaryOrigin = "Cosmic";
				    }
				else if (particleOrigin == simb::kSingleParticle)
				    {
					// Look for Beam particles: Beam + Muon halo
					if ( particle.Process() == "primary")
					    {
						primaryOrigin = "Beam";
					    }
					else if( particle.Process() == "primaryBackground")
					    {
						primaryOrigin = "Halo";
					    }
				    }
				else
				    {
					continue;
				    }

				nPrimaryMuons[0]++;
				if(primaryOrigin == "Cosmic")
				    {
					nPrimaryMuons[1]++;
				    }
				if(primaryOrigin == "Beam")
				    {
					nPrimaryMuons[2]++;
				    }
				if(primaryOrigin == "Halo")
				    {
					nPrimaryMuons[3]++;
				    }

				const TLorentzVector& primaryPositionStart = particle.Position(0);
				const TLorentzVector& primaryMomentumStart = particle.Momentum(0);
				TVector3 primaryStart(primaryPositionStart.X(), primaryPositionStart.Y(), primaryPositionStart.Z());

				double primaryAngleYZ = setAngle(atan2(primaryMomentumStart.Y(), primaryMomentumStart.Z()));
				if(verbo)
				    {
					logFile << endl;
					logFile << "Start position of primary muon (x, y, z, t): (" << primaryPositionStart.X() << ", "
						<< primaryPositionStart.Y() << ", " << primaryPositionStart.Z() << ", " << primaryT0 << ")" << endl;
					logFile << "TrackId: " << primaryTrackId << endl;
					logFile << "Energy: " << primaryEnergy << endl;
					logFile << "AngleYZ: " << primaryAngleYZ << endl;
				    }

				// Get CRT Hits
				vector< vector<double> > crtHitCounters;
				unsigned int primaryCountersHit = geoMuonCounter->testTrackInAllCounters(primaryTrackId,
													 primaryStart, primaryMomentumStart.Vect(), counterGeometry, crtHitCounters);

				if(primaryCountersHit != crtHitCounters.size())
				    {
					logFile << "ERROR: size of primary hit counters vector is not the same as number of primary hit counters." << endl;
				    }

				if(verbo)
				    {
					logFile << "Primary Hit Counter Size: " << crtHitCounters.size() << endl;
				    }

				if(crtHitCounters.size() == 1)
				    {
					nPrimaryMuonsWithOneHit[0]++;
					hPrimaryMuonsWithOneHit[0]->Fill(primaryEnergy);
					primaryMatchedRecoTrackOneHit[0].insert(make_pair(primaryTrackId, primaryEnergy));
					if(primaryOrigin == "Cosmic")
					    {
						nPrimaryMuonsWithOneHit[1]++;
						hPrimaryMuonsWithOneHit[1]->Fill(primaryEnergy);
						primaryMatchedRecoTrackOneHit[1].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
					if(primaryOrigin == "Beam")
					    {
						nPrimaryMuonsWithOneHit[2]++;
						hPrimaryMuonsWithOneHit[2]->Fill(primaryEnergy);
						primaryMatchedRecoTrackOneHit[2].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
					if(primaryOrigin == "Halo")
					    {
						nPrimaryMuonsWithOneHit[3]++;
						hPrimaryMuonsWithOneHit[3]->Fill(primaryEnergy);
						primaryMatchedRecoTrackOneHit[3].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
				    }
				if(crtHitCounters.size() == 2)
				    {
					nPrimaryMuonsWithTwoHits[0]++;
					hPrimaryMuonsWithTwoHits[0]->Fill(primaryEnergy);
					primaryMatchedRecoTrackTwoHits[0].insert(make_pair(primaryTrackId, primaryEnergy));
					if(primaryOrigin == "Cosmic")
					    {
						nPrimaryMuonsWithTwoHits[1]++;
						hPrimaryMuonsWithTwoHits[1]->Fill(primaryEnergy);
						primaryMatchedRecoTrackTwoHits[1].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
					if(primaryOrigin == "Beam")
					    {
						nPrimaryMuonsWithTwoHits[2]++;
						hPrimaryMuonsWithTwoHits[2]->Fill(primaryEnergy);
						primaryMatchedRecoTrackTwoHits[2].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
					if(primaryOrigin == "Halo")
					    {
						nPrimaryMuonsWithTwoHits[3]++;
						hPrimaryMuonsWithTwoHits[3]->Fill(primaryEnergy);
						primaryMatchedRecoTrackTwoHits[3].insert(make_pair(primaryTrackId, primaryEnergy));
					    }
				    }

				// Smear the position of hits with CRT resolutions
				if(verbo && crtHitCounters.size() > 0)
				    {
					logFile << "CRT hit position and T0 before smearing:" << endl;
					for(unsigned int hc = 0; hc < crtHitCounters.size(); hc++)
					    {
						logFile << "Intersection point for primary hit " << hc << " (x, y, z, t) = (" << crtHitCounters[hc][3] << ", "
							<< crtHitCounters[hc][4] << ", " << crtHitCounters[hc][5] << ", " << primaryT0 << ")" << endl;
					    }
				    }

				for(unsigned int hc = 0; hc < crtHitCounters.size(); hc++)
				    {
					for(unsigned int nd = 3; nd < crtHitCounters[hc].size(); nd++)
					    {
						uniform_real_distribution<double> xyzDistribution(-CRT_XYZ_RESOLUTION, CRT_XYZ_RESOLUTION);
						default_random_engine xyzGenerator;

						if((crtHitCounters[hc][0] <= 19 && nd != 5))
						    {
							crtHitCounters[hc][nd] = crtHitCounters[hc][nd] + xyzDistribution(xyzGenerator);
						    }
					    }

					TVector3 hitPosition(crtHitCounters[hc][3], crtHitCounters[hc][4], crtHitCounters[hc][5]);
					double distanceTravelled = (primaryStart - hitPosition).Mag() * 1E-2; //convert cm to m
					double velocity = (particle.P() / particle.E()) * TMath::C();
					double hitsT0 = primaryPositionStart.T() + ((distanceTravelled / (velocity)) * 1E9); //convert to ns

					hits tempHitsCRT;
					tempHitsCRT.tempId = tempId;

					tempHitsCRT.hitPositionX = crtHitCounters[hc][3];
					tempHitsCRT.hitPositionY = crtHitCounters[hc][4];
					tempHitsCRT.hitPositionZ = crtHitCounters[hc][5];

					uniform_real_distribution<double> tDistribution(-CRT_T_RESOLUTION, CRT_T_RESOLUTION);
					default_random_engine tGenerator;
					hitsT0 = hitsT0 + tDistribution(tGenerator);

					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					// TPCG4Time2Tick: Given G4 time returns TPC time-tick
					// { return (G4ToElecTime(g4time) - (TriggerTime() + TriggerOffsetTPC())) / fTPCClock.TickPeriod(); }
					// G4ToElecTime: Given Geant4 time [ns], returns relative time [us] w.r.t. electronics time T0
					// {return g4_time * 1.e-3 - fG4RefTime; }
					/////////////////////////////////////////////////////////////////////////////////////////////////////////
					tempHitsCRT.hitT0 = detectorClocksService->TPCG4Time2Tick(hitsT0) / detectorClocksService->TPCClock().Frequency(); //in us

					if(crtHitCounters[hc][0] <= 7)
					    {
						tempHitsCRT.hitFront = true;
						tempHitsCRT.hitBack = false;
						nHitsCRT_F++;
					    }
					else
					    {
						tempHitsCRT.hitFront = false;
						tempHitsCRT.hitBack = true;
						nHitsCRT_B++;
					    }

					tempHitsCRT.primaryTrackId = primaryTrackId;
					tempHitsCRT.primaryOrigin = primaryOrigin;
					tempHitsCRT.primaryEnergy = primaryEnergy;
					tempHitsCRT.primaryAngleYZ = primaryAngleYZ;

					hitsCRT.push_back(tempHitsCRT);
					tempId++;
					nHitsCRT++;
				    }//done with CRT counters
			    }// end loop over particles
		    }//MCTruthListHandle
	    }//MCTruthList

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(verbo)
	    {
		logFile << endl;
		logFile << "Done with CRT hit info!!" << endl;
		logFile << "Total CRT hits: " << hitsCRT.size() << endl << endl;
		if(verbo && hitsCRT.size() > 0)
		    {
			logFile << "CRT hit position and To after smearing:" << endl;
			for(unsigned int hc = 0; hc < hitsCRT.size(); hc++)
			    {
				logFile << "Intersection point for hit " << hc + 1 << " (x, y, z, t) = (" << hitsCRT[hc].hitPositionX << ", "
					<< hitsCRT[hc].hitPositionY << ", " << hitsCRT[hc].hitPositionZ << ", " << hitsCRT[hc].hitT0 << ")" << endl;
			    }
		    }

	    }
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Reconstructed tracks
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//TrackList handle
	art::Handle< std::vector<recob::Track> > trackListHandle;
	std::vector<art::Ptr<recob::Track> > trackList;
	if (event.getByLabel(fTrackModuleLabel, trackListHandle))
	    {
		art::fill_ptr_vector(trackList, trackListHandle);
	    }
	int nTracksReco = trackList.size();

	//FlashList handle
	art::Handle< std::vector<recob::OpFlash> > flashListHandle;
	std::vector<art::Ptr<recob::OpFlash> > flashList;
	if (event.getByLabel(fFlashModuleLabel, flashListHandle))
	    {
		art::fill_ptr_vector(flashList, flashListHandle);
	    }

	int nFlashes = flashList.size();
	for (int iFlash = 0; iFlash < nFlashes; iFlash++)
	    {
		double totalPE = flashList[iFlash]->TotalPE(); //The amount of light, in units of photo electrons, observed in this hit
		hPE->Fill(totalPE);

		double flashTime = flashList[iFlash]->Time(); //The time the hit occurred
		hFlashPeakTime->Fill(flashTime);
	    }

	if(verbo)
	    {
		logFile << endl;
		logFile << "Begining reco track info!!" << endl;
		logFile << "Total number of reco tracks: " << nTracksReco << endl << endl;
	    }

	//Find the associations between tracks and hits
	const art::FindManyP<recob::Hit> hitsFromTrack(trackListHandle, event, fTrackModuleLabel);

	vector<tracksPair> tracksPairF;
	vector<tracksPair> tracksPairB;
	for(int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack)
	    {
		nTotalRecoTracks++;
		std::vector< art::Ptr<recob::Hit> > allHits =  hitsFromTrack.at(iRecoTrack);

		unordered_map<int, double> trkIDE;
		for (auto const & hit : allHits)
		    {
			for (auto const & ide : backTracker->HitToTrackIDEs(hit)) // IDE: Ionization energy from a Geant4 track
			    {
				trkIDE[ide.trackID] += ide.energy; // Sum energy from the particle with Geant4 supplied trackID [MeV]
			    }
		    }

		int bestTrackId = 0;
		double totalTrackEnergy = 0.0, maxEnergy = 0.0;
		for (auto const & contrib : trkIDE)
		    {
			totalTrackEnergy += contrib.second; // Sum total energy in these hits
			if (contrib.second > maxEnergy) // Find track ID corresponding to max energy
			    {
				maxEnergy = contrib.second;
				bestTrackId = contrib.first;
			    }
		    }

		totalTrackEnergy *= 1E-3; //MeV to GeV

		if(bestTrackId < 0)  //-ve id means this is EM activity caused by track with the same but positive ID
		    {
			bestTrackId = -bestTrackId;
		    }

		if(primaryMatchedRecoTrackOneHit[0].find(bestTrackId) != primaryMatchedRecoTrackOneHit[0].end())
		    {
			hRecoTrackEnergy->Fill(totalTrackEnergy, primaryMatchedRecoTrackOneHit[0][bestTrackId]);
		    }

		int nTrajectoryPoints = trackList[iRecoTrack]->NumberTrajectoryPoints();
		int nHits = allHits.size();
		hHitsVsEnergy->Fill(totalTrackEnergy, nHits);

		// Don't use the very end points of the tracks (trackList[iRecoTrack]->End()) in case of scatter or distortion
		int nSkip = 2;
		int firstPoint = 0;
		int lastPoint = nTrajectoryPoints - nSkip; //End point is (numberTrajectoryPoints - 1)
		int firstHit = nHits - nSkip; // Last entry is first 'hit'
		int lastHit = 0;

		if(nHits < nSkip)
		    {
			continue;
		    }

		TVector3 trackStartPosition = trackList[iRecoTrack]->LocationAtPoint(firstPoint);;
		TVector3 trackEndPosition  = trackList[iRecoTrack]->LocationAtPoint(lastPoint);

		if(trackList[iRecoTrack]->End().Z() < trackList[iRecoTrack]->Vertex().Z()) //flip
		    {
			firstHit = 0;
			lastHit = nHits - nSkip;
			trackStartPosition = trackList[iRecoTrack]->LocationAtPoint(lastPoint);
			trackEndPosition = trackList[iRecoTrack]->LocationAtPoint(firstPoint);
		    }

		double trackStartPositionX = trackStartPosition.X();
		double trackStartPositionY = trackStartPosition.Y();
		double trackStartPositionZ = trackStartPosition.Z();
		double trackEndPositionX = trackEndPosition.X();
		double trackEndPositionY = trackEndPosition.Y();
		double trackEndPositionZ = trackEndPosition.Z();

		double trackLengthYZ = sqrt(pow(trackEndPositionY - trackStartPositionY, 2) + pow(trackEndPositionZ - trackStartPositionZ, 2));
		hRecoTrackLengthYZ->Fill(trackLengthYZ);
		hRecoTrackLengthYZVsEnergy->Fill(totalTrackEnergy, trackLengthYZ);

		if((trackLengthYZ < MIN_TRACK_LENGTH_YZ) || (totalTrackEnergy < MIN_TRACK_ENERGY))
		    {
			continue;
		    }

		//Find where this track originated from
		string trackOrigin;
		auto origin = particleInventory->TrackIdToMCTruth_P(bestTrackId)->Origin();
		if((fabs((particleInventory->TrackIdToParticle_P(bestTrackId))->PdgCode()) == fabs(fSelectedPDG)))
		    {
			if (origin == simb::kSingleParticle)
			    {
				// Look for Beam particles: Beam + Muon halo
				if(particleInventory->TrackIdToParticle_P(bestTrackId)->Process() == "primary")
				    {
					trackOrigin = "Halo";
				    }
				else
				    {
					trackOrigin = "Beam";
				    }
			    }
			else if (origin == simb::kCosmicRay)
			    {
				trackOrigin = "Cosmic";
			    }
			else
			    {
				trackOrigin = "Unknown";
			    }
		    }

		nConsideredRecoTracks[0]++;
		if(trackOrigin == "Cosmic")
		    {
			nConsideredRecoTracks[1]++;
		    }
		if(trackOrigin == "Beam")
		    {
			nConsideredRecoTracks[2]++;
		    }
		if(trackOrigin == "Halo")
		    {
			nConsideredRecoTracks[3]++;
		    }

		//Coverage of all reco tracks
		double trackSlopeYZ = (trackEndPositionY - trackStartPositionY) / (trackEndPositionZ - trackStartPositionZ);
		double yzIntercept = trackStartPositionY - (trackSlopeYZ * trackStartPositionZ);

		double trackSlopeYX = (trackEndPositionY - trackStartPositionY) / (trackEndPositionX - trackStartPositionX);
		double yxIntercept = trackStartPositionY - (trackSlopeYX * trackStartPositionX);

		double yHere = trackStartPositionY;
		double yLast = trackEndPositionY;

		if(yHere > yLast)
		    {
			yHere = trackEndPositionY;
			yLast = trackStartPositionY;
		    }

		while(yHere < yLast)
		    {
			double xHere = (yHere - yxIntercept) / trackSlopeYX;
			double zHere = (yHere - yzIntercept) / trackSlopeYZ;

			hRecoTracksCoverageYX[0]->Fill(xHere, yHere);
			hRecoTracksCoverageYZ[0]->Fill(zHere, yHere);

			if(trackOrigin == "Cosmic")
			    {
				hRecoTracksCoverageYX[1]->Fill(xHere, yHere);
				hRecoTracksCoverageYZ[1]->Fill(zHere, yHere);
			    }
			if(trackOrigin == "Beam")
			    {
				hRecoTracksCoverageYX[2]->Fill(xHere, yHere);
				hRecoTracksCoverageYZ[2]->Fill(zHere, yHere);
			    }
			if(trackOrigin == "Halo")
			    {
				hRecoTracksCoverageYX[3]->Fill(xHere, yHere);
				hRecoTracksCoverageYZ[3]->Fill(zHere, yHere);
			    }

			yHere += dY;
		    }

		double trackAngleYZ = setAngle(atan2(trackEndPositionY - trackStartPositionY, trackEndPositionZ - trackStartPositionZ));
		if(verbo)
		    {
			logFile << "Considered Reconstruced Track info: " << endl;
			logFile << "TrackId: " << bestTrackId << endl;
			logFile << "Length: " << trackList[iRecoTrack]->Length() << endl;
			logFile << "Energy: " << totalTrackEnergy << endl;
			logFile << "Start (x,y,z): " << trackStartPositionX << ", " << trackStartPositionY << ", " << trackStartPositionZ << endl;
			logFile << "End (x,y,z): " << trackEndPositionX << ", " << trackEndPositionY << ", " << trackEndPositionZ << endl;
			logFile << "AngleYZ: " << trackAngleYZ << endl;
			logFile << "Track origin: " << trackOrigin << endl << endl;
		    }

		if(primaryMatchedRecoTrackOneHit[0].find(bestTrackId) != primaryMatchedRecoTrackOneHit[0].end())
		    {
			nPrimaryMatchedRecoTracksOneHit[0]++;
			hPrimaryMatchedRecoTracksOneHit[0]->Fill(primaryMatchedRecoTrackOneHit[0][bestTrackId]);
		    }
		if((primaryMatchedRecoTrackOneHit[1].find(bestTrackId) != primaryMatchedRecoTrackOneHit[1].end()) && (trackOrigin == "Cosmic"))
		    {
			nPrimaryMatchedRecoTracksOneHit[1]++;
			hPrimaryMatchedRecoTracksOneHit[1]->Fill(primaryMatchedRecoTrackOneHit[1][bestTrackId]);
		    }
		if((primaryMatchedRecoTrackOneHit[2].find(bestTrackId) != primaryMatchedRecoTrackOneHit[2].end()) && (trackOrigin == "Beam"))
		    {
			nPrimaryMatchedRecoTracksOneHit[2]++;
			hPrimaryMatchedRecoTracksOneHit[2]->Fill(primaryMatchedRecoTrackOneHit[2][bestTrackId]);
		    }
		if((primaryMatchedRecoTrackOneHit[3].find(bestTrackId) != primaryMatchedRecoTrackOneHit[3].end()) && (trackOrigin == "Halo"))
		    {
			nPrimaryMatchedRecoTracksOneHit[3]++;
			hPrimaryMatchedRecoTracksOneHit[3]->Fill(primaryMatchedRecoTrackOneHit[3][bestTrackId]);
		    }

		if(primaryMatchedRecoTrackTwoHits[0].find(bestTrackId) != primaryMatchedRecoTrackTwoHits[0].end())
		    {
			nPrimaryMatchedRecoTracksTwoHits[0]++;
			hPrimaryMatchedRecoTracksTwoHits[0]->Fill(primaryMatchedRecoTrackTwoHits[0][bestTrackId]);
		    }

		if((primaryMatchedRecoTrackTwoHits[1].find(bestTrackId) != primaryMatchedRecoTrackTwoHits[1].end()) && (trackOrigin == "Cosmic"))
		    {
			nPrimaryMatchedRecoTracksTwoHits[1]++;
			hPrimaryMatchedRecoTracksTwoHits[1]->Fill(primaryMatchedRecoTrackTwoHits[1][bestTrackId]);
		    }
		if((primaryMatchedRecoTrackTwoHits[2].find(bestTrackId) != primaryMatchedRecoTrackTwoHits[2].end()) && (trackOrigin == "Beam"))
		    {
			nPrimaryMatchedRecoTracksTwoHits[2]++;
			hPrimaryMatchedRecoTracksTwoHits[2]->Fill(primaryMatchedRecoTrackTwoHits[2][bestTrackId]);
		    }
		if((primaryMatchedRecoTrackTwoHits[3].find(bestTrackId) != primaryMatchedRecoTrackTwoHits[3].end()) && (trackOrigin == "Halo"))
		    {
			nPrimaryMatchedRecoTracksTwoHits[3]++;
			hPrimaryMatchedRecoTracksTwoHits[3]->Fill(primaryMatchedRecoTrackTwoHits[3][bestTrackId]);
		    }

		for(unsigned int iHitsCRT = 0; iHitsCRT < hitsCRT.size(); iHitsCRT++)
		    {
			// First constrain in Y-Z plane
			double predictedHitPositionY = (trackSlopeYZ * hitsCRT[iHitsCRT].hitPositionZ) + yzIntercept;

			double deltaY = fabs(hitsCRT[iHitsCRT].hitPositionY - predictedHitPositionY);
			hDeltaY->Fill(deltaY);

			if(deltaY > MAX_DELTA_Y)
			    {
				continue;
			    }

			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// GetXTicksOffset - TriggerOffset: The offsets for each plane
			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			double g4TicksStart = (hitsCRT[iHitsCRT].hitT0 * detectorClocksService->TPCClock().Frequency()) +
			    detectorPropertiesService->GetXTicksOffset(allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat) -
			    detectorPropertiesService->TriggerOffset();
			double g4TicksEnd = (hitsCRT[iHitsCRT].hitT0 * detectorClocksService->TPCClock().Frequency()) +
			    detectorPropertiesService->GetXTicksOffset(allHits[lastHit]->WireID().Plane, allHits[lastHit]->WireID().TPC, allHits[lastHit]->WireID().Cryostat) -
			    detectorPropertiesService->TriggerOffset();

			double xOffsetStart = detectorPropertiesService->ConvertTicksToX(g4TicksStart,
											 allHits[firstHit]->WireID().Plane, allHits[firstHit]->WireID().TPC, allHits[firstHit]->WireID().Cryostat);
			double xOffsetEnd = detectorPropertiesService->ConvertTicksToX(g4TicksEnd,
										       allHits[lastHit]->WireID().Plane, allHits[lastHit]->WireID().TPC, allHits[lastHit]->WireID().Cryostat);

			// Move the track start x-position to correspond to T0
			double trackNewStartPositionX = trackStartPositionX - xOffsetStart;
			double trackNewEndPositionX = trackEndPositionX - xOffsetEnd;

			//Make sure the new start and end track positions are inside the TPC
			art::ServiceHandle<geo::Geometry> geom;
			double trackStartPositionArray[3] = {trackNewStartPositionX, trackStartPositionY, trackStartPositionZ};
			double trackEndPositionArray[3] = {trackNewEndPositionX, trackEndPositionY, trackEndPositionZ};
			geo::TPCID tpcIdStart = geom->FindTPCAtPosition(trackStartPositionArray);
			geo::TPCID tpcIdEnd = geom->FindTPCAtPosition(trackEndPositionArray);

			if(!tpcIdStart.isValid || !tpcIdEnd.isValid)
			    {
				continue;
			    }

			double trackNewSlopeXZ = (trackNewEndPositionX - trackNewStartPositionX) / (trackEndPositionZ - trackStartPositionZ);
			double xIntercept = trackNewStartPositionX - (trackNewSlopeXZ * trackStartPositionZ);
			double predictedHitPositionX = (trackNewSlopeXZ * hitsCRT[iHitsCRT].hitPositionZ) + xIntercept;
			double deltaX = fabs(hitsCRT[iHitsCRT].hitPositionX - predictedHitPositionX);

			tracksPair tempTracksPair;
			tempTracksPair.index.first = iRecoTrack;
			tempTracksPair.index.second = hitsCRT[iHitsCRT].tempId;
			tempTracksPair.trackId.first = bestTrackId;
			tempTracksPair.trackId.second = hitsCRT[iHitsCRT].primaryTrackId;
			tempTracksPair.hitFB.first = hitsCRT[iHitsCRT].hitFront;
			tempTracksPair.hitFB.second = hitsCRT[iHitsCRT].hitBack;
			tempTracksPair.T0.first = allHits[firstHit]->PeakTime() / detectorClocksService->TPCClock().Frequency();
			tempTracksPair.T0.second = hitsCRT[iHitsCRT].hitT0;
			tempTracksPair.AngleYZ.first = trackAngleYZ;
			tempTracksPair.AngleYZ.second = hitsCRT[iHitsCRT].primaryAngleYZ;
			tempTracksPair.energy.first = totalTrackEnergy;
			tempTracksPair.energy.second = hitsCRT[iHitsCRT].primaryEnergy;
			tempTracksPair.origin.first = trackOrigin;
			tempTracksPair.origin.second = hitsCRT[iHitsCRT].primaryOrigin;

			tempTracksPair.slopeYX = (trackEndPositionY - trackStartPositionY) / (trackNewEndPositionX - trackNewStartPositionX);
			tempTracksPair.slopeYZ = trackSlopeYZ;
			tempTracksPair.startPositionArray[0] = trackStartPositionArray[0];
			tempTracksPair.startPositionArray[1] = trackStartPositionArray[1];
			tempTracksPair.startPositionArray[2] = trackStartPositionArray[2];
			tempTracksPair.endPositionArray[0] = trackEndPositionArray[0];
			tempTracksPair.endPositionArray[1] = trackEndPositionArray[1];
			tempTracksPair.endPositionArray[2] = trackEndPositionArray[2];

			tempTracksPair.deltaX = deltaX;

			if(hitsCRT[iHitsCRT].hitFront)
			    {
				tracksPairF.push_back(tempTracksPair);
			    }
			if(hitsCRT[iHitsCRT].hitBack)
			    {
				tracksPairB.push_back(tempTracksPair);
			    }
		    }//iHitsCRT
	    }//iRecoTrack

	//Sort pair by ascending order of deltaR
	sort(tracksPairF.begin(), tracksPairF.end(), sortByR());
	sort(tracksPairB.begin(), tracksPairB.end(), sortByR());

	//Require 1 to 1 matching- and save as unique pair
	vector<tracksPair> allUniqueTracksPair;
	while ( tracksPairF.size())
	    {
		allUniqueTracksPair.push_back(tracksPairF.front());
		tracksPairF.erase(remove_if(tracksPairF.begin(), tracksPairF.end(), removePairIndex(tracksPairF.front())),
				  tracksPairF.end());
	    }
	while ( tracksPairB.size())
	    {
		allUniqueTracksPair.push_back(tracksPairB.front());
            tracksPairB.erase(remove_if(tracksPairB.begin(), tracksPairB.end(), removePairIndex(tracksPairB.front())),
                              tracksPairB.end());
        }

    //Sort pair by ascending order of track index- will be used not to double count
    sort(allUniqueTracksPair.begin(), allUniqueTracksPair.end(), sortByIndex());

    //Keep track of tracks with two hits (front and back)
    unordered_map<int, int> twoHitsTracks;
    for(unsigned int u = 0; u < allUniqueTracksPair.size(); u++)
        {
            twoHitsTracks[allUniqueTracksPair[u].index.first]++;
        }

    //Access information
    if(verbo)
        {
            logFile << endl;
            logFile << "Done with reco track info!!" << endl;
            logFile << "Begining matching info!!" << endl << endl;
            logFile << "Unique pair size: " << allUniqueTracksPair.size() << endl;
        }

    for(unsigned int u = 0; u < allUniqueTracksPair.size(); u++)
        {
            //Don't double count; skip the next pair- it's the same track matched to differnet CRT plane
            bool oneHitTrack = true;
            if(twoHitsTracks[allUniqueTracksPair[u].index.first] == 2)
                {
                    oneHitTrack = false;
                    u++;
                }

            double deltaX = allUniqueTracksPair[u].deltaX;
            hDeltaX->Fill(deltaX);

            if(deltaX > MAX_DELTA_X)
                {
                    continue;
                }

            double minTimeDifference = 9999.99;
            for (int iFlash = 0; iFlash < nFlashes; iFlash++)
                {
                    if (flashList[iFlash]->TotalPE() < MIN_PE)
                        {
                            continue;
                        }

                    double flashTime = flashList[iFlash]->Time();
                    double timeDifference = allUniqueTracksPair[u].T0.second - flashTime;

                    if(fabs(timeDifference) < fabs(minTimeDifference))
                        {
                            minTimeDifference = timeDifference;
                        }
                }
            hDeltaT->Fill(minTimeDifference);

            if((minTimeDifference > MAX_DELTA_T) || (minTimeDifference < MIN_DELTA_T) )
                {
                    continue;
                }

            if(oneHitTrack)
                {
                    nAllCRTMatchedRecoTracksOneHit[0]++;
                    hAllCRTMatchedRecoTracksOneHit[0]->Fill(allUniqueTracksPair[u].energy.second);
                    if(allUniqueTracksPair[u].origin.first == "Cosmic")
                        {
                            nAllCRTMatchedRecoTracksOneHit[1]++;
                            hAllCRTMatchedRecoTracksOneHit[1]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                    if(allUniqueTracksPair[u].origin.first == "Beam")
                        {
                            nAllCRTMatchedRecoTracksOneHit[2]++;
                            hAllCRTMatchedRecoTracksOneHit[2]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                    if(allUniqueTracksPair[u].origin.first == "Halo")
                        {
                            nAllCRTMatchedRecoTracksOneHit[3]++;
                            hAllCRTMatchedRecoTracksOneHit[3]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                }
            else
                {
                    nAllCRTMatchedRecoTracksTwoHits[0]++;
                    hAllCRTMatchedRecoTracksTwoHits[0]->Fill(allUniqueTracksPair[u].energy.second);
                    if(allUniqueTracksPair[u].origin.first == "Cosmic")
                        {
                            nAllCRTMatchedRecoTracksTwoHits[1]++;
                            hAllCRTMatchedRecoTracksTwoHits[1]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                    if(allUniqueTracksPair[u].origin.first == "Beam")
                        {
                            nAllCRTMatchedRecoTracksTwoHits[2]++;
                            hAllCRTMatchedRecoTracksTwoHits[2]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                    if(allUniqueTracksPair[u].origin.first == "Halo")
                        {
                            nAllCRTMatchedRecoTracksTwoHits[3]++;
                            hAllCRTMatchedRecoTracksTwoHits[3]->Fill(allUniqueTracksPair[u].energy.second);
                        }
                }
            //Look at the coverage of the matched tracks
            double yxIntercept = allUniqueTracksPair[u].startPositionArray[1] - (allUniqueTracksPair[u].slopeYX * allUniqueTracksPair[u].startPositionArray[0]);
            double yzIntercept = allUniqueTracksPair[u].startPositionArray[1] - (allUniqueTracksPair[u].slopeYZ * allUniqueTracksPair[u].startPositionArray[2]);

            double yHere = allUniqueTracksPair[u].startPositionArray[1];
            double yLast = allUniqueTracksPair[u].endPositionArray[1];

            if(yHere > yLast)
                {
                    yHere = allUniqueTracksPair[u].endPositionArray[1];
                    yLast = allUniqueTracksPair[u].startPositionArray[1];
                }

            while(yHere < yLast)
                {
                    double xHere = (yHere - yxIntercept) / allUniqueTracksPair[u].slopeYX;
                    double zHere = (yHere - yzIntercept) / allUniqueTracksPair[u].slopeYZ;

                    if(oneHitTrack)
                        {
                            hCoverageYXOneHit[0]->Fill(xHere, yHere);
                            hCoverageYZOneHit[0]->Fill(zHere, yHere);
                            if(allUniqueTracksPair[u].origin.first == "Cosmic")
                                {
                                    hCoverageYXOneHit[1]->Fill(xHere, yHere);
                                    hCoverageYZOneHit[1]->Fill(zHere, yHere);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Beam")
                                {
                                    hCoverageYXOneHit[2]->Fill(xHere, yHere);
                                    hCoverageYZOneHit[2]->Fill(zHere, yHere);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Halo")
                                {
                                    hCoverageYXOneHit[3]->Fill(xHere, yHere);
                                    hCoverageYZOneHit[3]->Fill(zHere, yHere);
                                }
                        }
                    else
                        {
                            hCoverageYXTwoHits[0]->Fill(xHere, yHere);
                            hCoverageYZTwoHits[0]->Fill(zHere, yHere);
                            if(allUniqueTracksPair[u].origin.first == "Cosmic")
                                {
                                    hCoverageYXTwoHits[1]->Fill(xHere, yHere);
                                    hCoverageYZTwoHits[1]->Fill(zHere, yHere);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Beam")
                                {
                                    hCoverageYXTwoHits[2]->Fill(xHere, yHere);
                                    hCoverageYZTwoHits[2]->Fill(zHere, yHere);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Halo")
                                {
                                    hCoverageYXTwoHits[3]->Fill(xHere, yHere);
                                    hCoverageYZTwoHits[3]->Fill(zHere, yHere);
                                }
                        }
                    yHere += dY;
                }

            if(verbo)
                {
                    logFile << "All CRT Matched info (track, primary): " << endl;
                    logFile << "TrackId: (" << allUniqueTracksPair[u].trackId.first << ", " << allUniqueTracksPair[u].trackId.second << ")" << endl;
                    logFile << "AngleYZ: (" << allUniqueTracksPair[u].AngleYZ.first << ", " << allUniqueTracksPair[u].AngleYZ.second << ")" << endl;
                    logFile << "Energy: (" << allUniqueTracksPair[u].energy.first << ", " << allUniqueTracksPair[u].energy.second << ")" << endl << endl;
                }

            if(oneHitTrack)
                {
                    if(allUniqueTracksPair[u].trackId.first == allUniqueTracksPair[u].trackId.second)
                        {
                            nGoodCRTMatchedRecoTracksOneHit[0]++;
                            hGoodCRTMatchedRecoTracksOneHit[0]->Fill(allUniqueTracksPair[u].energy.second);
                            if(allUniqueTracksPair[u].origin.first == "Cosmic")
                                {
                                    nGoodCRTMatchedRecoTracksOneHit[1]++;
                                    hGoodCRTMatchedRecoTracksOneHit[1]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Beam")
                                {
                                    nGoodCRTMatchedRecoTracksOneHit[2]++;
                                    hGoodCRTMatchedRecoTracksOneHit[2]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Halo")
                                {
                                    nGoodCRTMatchedRecoTracksOneHit[3]++;
                                    hGoodCRTMatchedRecoTracksOneHit[3]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                        }
                }
            else
                {
                    if((allUniqueTracksPair[u].trackId.first == allUniqueTracksPair[u].trackId.second) &&
                            (allUniqueTracksPair[u + 1].trackId.first == allUniqueTracksPair[u + 1].trackId.second))
                        {
                            nGoodCRTMatchedRecoTracksTwoHits[0]++;
                            hGoodCRTMatchedRecoTracksTwoHits[0]->Fill(allUniqueTracksPair[u].energy.second);
                            if(allUniqueTracksPair[u].origin.first == "Cosmic")
                                {
                                    nGoodCRTMatchedRecoTracksTwoHits[1]++;
                                    hGoodCRTMatchedRecoTracksTwoHits[1]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Beam")
                                {
                                    nGoodCRTMatchedRecoTracksTwoHits[2]++;
                                    hGoodCRTMatchedRecoTracksTwoHits[2]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                            if(allUniqueTracksPair[u].origin.first == "Halo")
                                {
                                    nGoodCRTMatchedRecoTracksTwoHits[3]++;
                                    hGoodCRTMatchedRecoTracksTwoHits[3]->Fill(allUniqueTracksPair[u].energy.second);
                                }
                        }
                }
        }
}//analyze
}// namespace CRTMatching

