////////////////////////////////////////////////////////////////////////
// Class:        MatchingWithTwoCRTs
// Module Type:  analyzer
// File:         MatchingWithTwoCRTs_module.cc
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
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// C++ includes
#include <cmath>
#include <random>

// MuonCounter Geometry includes
#include "dune/Geometry/MuonCounter35Alg.h"

// Header file
#include "MatchingWithTwoCRTs_module.h"

#include "../MCTruthInformation.h"

using namespace std;

namespace CRTMatching
{
    /////////////////////////////////////////////////////////
    //Constructor
    /////////////////////////////////////////////////////////
    MatchingWithTwoCRTs::MatchingWithTwoCRTs(fhicl::ParameterSet const& parameterSet)
	: EDAnalyzer(parameterSet)
    {
	// Read in the parameters from the .fcl file
	this->reconfigure(parameterSet);
    }


    /////////////////////////////////////////////////////////
    //Destructor
    /////////////////////////////////////////////////////////
    MatchingWithTwoCRTs::~MatchingWithTwoCRTs() {}

    /////////////////////////////////////////////////////////
    //Reads parameters form the .fcl file
    /////////////////////////////////////////////////////////
    void MatchingWithTwoCRTs::reconfigure(fhicl::ParameterSet const& p)
    {

	fTrackModuleLabel = p.get< string >("TrackModuleLabel");
	fSelectedPDG = p.get< int >("PDGcode");

	return;
    }

    /////////////////////////////////////////////////////////
    //Executes once at the beginning of the job
    /////////////////////////////////////////////////////////
    void MatchingWithTwoCRTs::beginJob()
    {
	nEvents = 0, nPrimaryMuons = 0, nPrimaryMuonsWithTwoHits = 0;
	nPrimaryHits_F = 0, nPrimaryHits_B = 0, nPrimaryHits_T = 0;
	nPrimaryHits_FB = 0, nPrimaryHits_FT = 0, nPrimaryHits_BT = 0;
	nCombinatorialTracks_FB = 0, nCombinatorialTracks_FT = 0, nCombinatorialTracks_BT = 0;
	nTotalRecoTracks = 0, nConsideredRecoTracks = 0, nPrimaryMatchedRecoTracks = 0;
	nAllCRTMatchedRecoTracks = 0, nGoodCRTMatchedRecoTracks = 0;

	// Open a basic log file, will overwrite a pre-existing one
	logFile.open("MatchingWithTwoCRTs.log");

	// Get local time
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	logFile << "MatchingWithTwoCRTs_module log file, " << asctime(timeinfo) << endl;

	// Load the CRT positions from a text file
	char counterfile[] = "/nashome/a/arbint/DuneApp/SpaceChargeEffects/CRTMatching/CRTs.txt";
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
	hStatistics = fileServiceHandle->make<TH1D>("hStatistics", "Job Statistics", 50, 0, 50);

	hPrimaryAngleYZ_FB = fileServiceHandle->make<TH1D>("hPrimaryAngleYZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hPrimaryAngleYZ_FT = fileServiceHandle->make<TH1D>("hPrimaryAngleYZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hPrimaryAngleYZ_BT = fileServiceHandle->make<TH1D>("hPrimaryAngleYZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hPrimaryAngleXZ_FB = fileServiceHandle->make<TH1D>("hPrimaryAngleXZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hPrimaryAngleXZ_FT = fileServiceHandle->make<TH1D>("hPrimaryAngleXZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hPrimaryAngleXZ_BT = fileServiceHandle->make<TH1D>("hPrimaryAngleXZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);

	hCombinatorialAngleYZ_FB = fileServiceHandle->make<TH1D>("hCombinatorialAngleYZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialAngleYZ_FT = fileServiceHandle->make<TH1D>("hCombinatorialAngleYZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialAngleYZ_BT = fileServiceHandle->make<TH1D>("hCombinatorialAngleYZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialAngleXZ_FB = fileServiceHandle->make<TH1D>("hCombinatorialAngleXZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialAngleXZ_FT = fileServiceHandle->make<TH1D>("hCombinatorialAngleXZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialAngleXZ_BT = fileServiceHandle->make<TH1D>("hCombinatorialAngleXZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);

	hCombinatorialPrimaryAngleYZ_FB = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleYZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialPrimaryAngleYZ_FT = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleYZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialPrimaryAngleYZ_BT = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleYZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialPrimaryAngleXZ_FB = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleXZ_FB", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialPrimaryAngleXZ_FT = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleXZ_FT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);
	hCombinatorialPrimaryAngleXZ_BT = fileServiceHandle->make<TH1D>("hCombinatorialPrimaryAngleXZ_BT", "", BINS_ANGLE, MIN_BINS_ANGLE, MAX_BINS_ANGLE);

	hRecoTrackLengthYZ = fileServiceHandle->make<TH1D>("hRecoTrackLengthYZ", "", 160, 0, 800.0);
	hRecoTrackLengthYZVsEnergy = fileServiceHandle->make<TH2D>("hRecoTrackLengthYZVsEnergy", "", 100, 0, 5.0, 160, 0, 800.0);

	hDeltaSlopeYZ = fileServiceHandle->make<TH1D>("hDeltaSlopeYZ", "", 100, 0.0, 1.0);
	hDeltaSlopeXZ = fileServiceHandle->make<TH1D>("hDeltaSlopeXZ", "", 100, 0.0, 1.0);
	hDeltaY1 = fileServiceHandle->make<TH1D>("hDeltaY1", "", 100, 0.0, 100.0);
	hDeltaY2 = fileServiceHandle->make<TH1D>("hDeltaY2", "", 100, 0.0, 100.0);
	hDeltaX = fileServiceHandle->make<TH1D>("hDeltaX", "", 500, 0.0, 500.0);

	hConsideredRecoTrack_RecoEnergy = fileServiceHandle->make<TH1D>("hConsideredRecoTrack_RecoEnergy", "", 40, 0, 20.0);
	hPrimaryMatchedRecoTrack_RecoEnergy = fileServiceHandle->make<TH1D>("hPrimaryMatchedRecoTrack_RecoEnergy", "", 40, 0, 20.0);
	hPrimaryMatchedRecoTrack_TrueEnergy = fileServiceHandle->make<TH1D>("hPrimaryMatchedRecoTrack_TrueEnergy", "", 40, 0, 20.0);
	hAllCRTMatchedRecoTrack_RecoEnergy = fileServiceHandle->make<TH1D>("hAllCRTMatchedRecoTrack_RecoEnergy", "", 40, 0, 20.0);
	hAllCRTMatchedRecoTrack_TrueEnergy = fileServiceHandle->make<TH1D>("hAllCRTMatchedRecoTrack_TrueEnergy", "", 40, 0, 20.0);
	hGoodCRTMatchedRecoTrack_RecoEnergy = fileServiceHandle->make<TH1D>("hGoodCRTMatchedRecoTrack_RecoEnergy", "", 40, 0, 20.0);
	hGoodCRTMatchedRecoTrack_TrueEnergy = fileServiceHandle->make<TH1D>("hGoodCRTMatchedRecoTrack_TrueEnergy", "", 40, 0, 20.0);
    }

    /////////////////////////////////////////////////////////
    //Executes once at the end of the job
    /////////////////////////////////////////////////////////
    void MatchingWithTwoCRTs::endJob()
    {
	hStatistics->SetBinContent(1, nEvents);

	hStatistics->SetBinContent(11, nPrimaryMuons);
	hStatistics->SetBinContent(12, nPrimaryMuonsWithTwoHits);
	hStatistics->SetBinContent(13, nPrimaryHits_F);
	hStatistics->SetBinContent(14, nPrimaryHits_B);
	hStatistics->SetBinContent(15, nPrimaryHits_T);
	hStatistics->SetBinContent(16, nPrimaryHits_FB);
	hStatistics->SetBinContent(17, nPrimaryHits_FT);
	hStatistics->SetBinContent(18, nPrimaryHits_BT);
	hStatistics->SetBinContent(19, nCombinatorialTracks_FB);
	hStatistics->SetBinContent(20, nCombinatorialTracks_FT);
	hStatistics->SetBinContent(21, nCombinatorialTracks_BT);

	hStatistics->SetBinContent(31, nTotalRecoTracks);
	hStatistics->SetBinContent(32, nConsideredRecoTracks);
	hStatistics->SetBinContent(33, nPrimaryMatchedRecoTracks);
	hStatistics->SetBinContent(34, nAllCRTMatchedRecoTracks);
	hStatistics->SetBinContent(35, nGoodCRTMatchedRecoTracks);

	cout << endl << endl;
	cout << "Primary Muons: " << nPrimaryMuons << endl;
	cout << "Primary Muons with hit: " << nPrimaryMuonsWithTwoHits << endl;
	cout << "CRT hits: " << nPrimaryHits_F + nPrimaryHits_B + nPrimaryHits_T << endl;
	cout << "Track Reconstruction Efficiency: " << 100.0 * nPrimaryMatchedRecoTracks / nPrimaryMuonsWithTwoHits << " %" << endl;
	cout << "Matching Purity: " << 100.0 * nGoodCRTMatchedRecoTracks / nAllCRTMatchedRecoTracks << " %" << endl;
	cout << "Matching Efficiency: " << 100.0 * nGoodCRTMatchedRecoTracks / nPrimaryMatchedRecoTracks << " %" << endl;
	cout << "Total Efficiency: " << 100.0 * nGoodCRTMatchedRecoTracks / nConsideredRecoTracks << " %" << endl;
	cout << endl << endl;

	logFile << "endJob called, all done!!!" << endl;
    }

    /////////////////////////////////////////////////////////
    //Executes once per event
    /////////////////////////////////////////////////////////
    void MatchingWithTwoCRTs::analyze( const art::Event & event )
    {
	nEvents++;

	//Clear the vectors for each event
	primaryHits_F.clear();
	primaryHits_B.clear();
	primaryHits_T.clear();

	// The basics
	int fEvent  = event.id().event();
	int fRun    = event.run();
	int fSubRun = event.subRun();
	logFile << "Event " << fEvent << ", run " << fRun << ", subrun " << fSubRun << endl << endl;

	verbo = (fEvent == 1);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Store primary hits
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(verbo)
	    {
		logFile << "Beginning CRT hit info!!" << endl;
	    }
	//Get the MC truth information
	const sim::ParticleList& plist = particleInventory->ParticleList();

	int nGEANTparticles = plist.size();
	if(verbo)
	    {
		logFile << "Number of geant particles: " << nGEANTparticles << endl;
	    }

	// Loop over all the particles, find the primary muon, and get CRT hits
	int tempId = 0;

	// Map of trackID  to see how primary were matched to reco tracks
	map<int, double> primaryMatchedRecoTrack;
	for (const auto& PartPair : plist)
	    {
		const simb::MCParticle& particle = *(PartPair.second);

		int primaryTrackId = particle.TrackId();
		int primaryPDG = particle.PdgCode();
		double primaryEnergy = particle.E();

		// Locate the primary muon in readout window
		bool inReadOutWindow = MCTruthInformation::InReadoutWindow(particle);
		if ( particle.Process() == "primary"  &&  fabs(primaryPDG) == fabs(fSelectedPDG) && inReadOutWindow)
		    {
			nPrimaryMuons++;

			const TLorentzVector& primaryPositionStart = particle.Position(0);
			const TLorentzVector& primaryMomentumStart = particle.Momentum(0);
			TVector3 primaryStart(primaryPositionStart.X(), primaryPositionStart.Y(), primaryPositionStart.Z());

			double primaryAngleYZ = setAngle(atan2(primaryMomentumStart.Y(), primaryMomentumStart.Z()));
			double primaryAngleXZ = setAngle(atan2(primaryMomentumStart.X(), primaryMomentumStart.Z()));

			if(verbo)
			    {
				logFile << endl;
				logFile << "Start position of primary muon (x, y, z): (" << primaryPositionStart.X() << ", "
					<< primaryPositionStart.Y() << ", " << primaryPositionStart.Z() << ")" << endl;
				logFile << "TrackId of primary muon: " << primaryTrackId << endl;
				logFile << "Energy of primary muon: " << primaryEnergy << endl;
				logFile << "Angle YZ of primary muon: " << primaryAngleYZ << endl;
				logFile << "Angle XZ of primary muon: " << primaryAngleXZ << endl << endl;
			    }

			// Get CRT Hits
			vector< vector<double> > primaryHitCounters;
			unsigned int primaryCountersHit = geoMuonCounter->testTrackInAllCounters(primaryTrackId,
												 primaryStart, primaryMomentumStart.Vect(), counterGeometry, primaryHitCounters);

			if(primaryCountersHit != primaryHitCounters.size())
			    {
				logFile << "ERROR: size of primary hit counters vector is not the same as number of primary hit counters." << endl;
			    }

			if(verbo)
			    {
				logFile << "Primary Hit Counter Size: " << primaryHitCounters.size() << endl;
			    }


			if(primaryHitCounters.size() == 2)
			    {
				nPrimaryMuonsWithTwoHits++;
				primaryMatchedRecoTrack.insert(make_pair(primaryTrackId, primaryEnergy));
			    }

			//Smear the position of hits with CRT resolutions
			if(verbo && primaryHitCounters.size() > 0)
			    {
				logFile << "CRT hit position before smearing:" << endl;
				for(unsigned int hc = 0; hc < primaryHitCounters.size(); hc++)
				    {
					logFile << "Intersection point for primary hit " << hc << " (x, y, z) = (" << primaryHitCounters[hc][3] << ", "
						<< primaryHitCounters[hc][4] << ", " << primaryHitCounters[hc][5] << ")" << endl;
				    }
			    }
			bool primaryFrontHit = false, primaryBackHit = false, primaryTopHit = false;
			for(unsigned int hc = 0; hc < primaryHitCounters.size(); hc++)
			    {
				for(unsigned int nd = 3; nd < primaryHitCounters[hc].size(); nd++)
				    {
					normal_distribution<double> distribution(primaryHitCounters[hc][nd], CRT_RESOLUTION);
					default_random_engine generator;

					if((primaryHitCounters[hc][0] <= 19 && nd != 5) || (primaryHitCounters[hc][0] > 19 && nd != 4))
					    {
						primaryHitCounters[hc][nd] = distribution(generator);
					    }
				    }

				hits tempPrimaryHits;
				tempPrimaryHits.tempId = tempId;

				tempPrimaryHits.hitPositionX = primaryHitCounters[hc][3];
				tempPrimaryHits.hitPositionY = primaryHitCounters[hc][4];
				tempPrimaryHits.hitPositionZ = primaryHitCounters[hc][5];

				tempPrimaryHits.primaryTrackId = primaryTrackId;
				tempPrimaryHits.primaryEnergy = primaryEnergy;
				tempPrimaryHits.primaryAngleYZ = primaryAngleYZ;
				tempPrimaryHits.primaryAngleXZ = primaryAngleXZ;

				if(primaryHitCounters[hc][0] <= 9)
				    {
					primaryHits_F.push_back(tempPrimaryHits);
					primaryFrontHit = true;
					nPrimaryHits_F++;
				    }
				if((primaryHitCounters[hc][0] >= 10) && (primaryHitCounters[hc][0] <= 19))
				    {
					primaryHits_B.push_back(tempPrimaryHits);
					primaryBackHit = true;
					nPrimaryHits_B++;
				    }
				if((primaryHitCounters[hc][0] >= 20) && (primaryHitCounters[hc][0] <= 29))
				    {
					primaryHits_T.push_back(tempPrimaryHits);
					primaryTopHit = true;
					nPrimaryHits_T++;
				    }
				tempId++;
			    }

			if(primaryFrontHit && primaryBackHit)
			    {
				nPrimaryHits_FB++;
				hPrimaryAngleYZ_FB->Fill(primaryAngleYZ);
				hPrimaryAngleXZ_FB->Fill(primaryAngleXZ);
			    }
			if(primaryFrontHit && primaryTopHit)
			    {
				nPrimaryHits_FT++;
				hPrimaryAngleYZ_FT->Fill(primaryAngleYZ);
				hPrimaryAngleXZ_FT->Fill(primaryAngleXZ);
			    }
			if(primaryBackHit && primaryTopHit)
			    {
				nPrimaryHits_BT++;
				hPrimaryAngleYZ_BT->Fill(primaryAngleYZ);
				hPrimaryAngleXZ_BT->Fill(primaryAngleXZ);
			    }

			if(verbo && primaryHitCounters.size() > 0)
			    {
				logFile << "CRT hit position after smearing:" << endl;
				for(unsigned int hc = 0; hc < primaryHitCounters.size(); hc++)
				    {
					logFile << "Intersection point for hit " << hc << " (x, y, z) = (" << primaryHitCounters[hc][3] << ", "
						<< primaryHitCounters[hc][4] << ", " << primaryHitCounters[hc][5] << ")" << endl;
				    }
			    }
			if(verbo && primaryHitCounters.size() == 2)
			    {
				double hitAngleYZ = setAngle(atan2(primaryHitCounters[0][4] - primaryHitCounters[1][4], primaryHitCounters[0][5] - primaryHitCounters[1][5]));
				double hitAngleXZ = setAngle(atan2(primaryHitCounters[0][3] - primaryHitCounters[1][3], primaryHitCounters[0][5] - primaryHitCounters[1][5]));
				logFile << "Angle YZ of CRT hits: " << hitAngleYZ << endl;
				logFile << "Angle XZ of CRT hits: " << hitAngleXZ << endl;
			    }
		    }// done with primary muon
	    } // end loop over particles
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(verbo)
	    {
		logFile << endl;
		logFile << "Done with CRT hit info!!" << endl << endl << endl << endl;
		logFile << "Beginning  combinatorial track info!!" << endl << endl;;
	    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Combinatorial tracks
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vector<combinatorialTracks> allCombinatorialTracks;
	for(unsigned int f = 0; f < primaryHits_F.size(); f++)
	    {
		for(unsigned int b = 0; b < primaryHits_B.size(); b++)
		    {
			combinatorialTracks tempCombinatorialTrack;
			tempCombinatorialTrack.tempId.first = primaryHits_F[f].tempId;
			tempCombinatorialTrack.tempId.second = primaryHits_B[b].tempId;

			tempCombinatorialTrack.hitPositionX.first = primaryHits_F[f].hitPositionX;
			tempCombinatorialTrack.hitPositionX.second = primaryHits_B[b].hitPositionX;
			tempCombinatorialTrack.hitPositionY.first = primaryHits_F[f].hitPositionY;
			tempCombinatorialTrack.hitPositionY.second = primaryHits_B[b].hitPositionY;
			tempCombinatorialTrack.hitPositionZ.first = primaryHits_F[f].hitPositionZ;
			tempCombinatorialTrack.hitPositionZ.second = primaryHits_B[b].hitPositionZ;

			tempCombinatorialTrack.primaryTrackId.first = primaryHits_F[f].primaryTrackId;
			tempCombinatorialTrack.primaryTrackId.second = primaryHits_B[b].primaryTrackId;
			tempCombinatorialTrack.primaryEnergy.first = primaryHits_F[f].primaryEnergy;
			tempCombinatorialTrack.primaryEnergy.second = primaryHits_B[b].primaryEnergy;
			tempCombinatorialTrack.primaryAngleYZ.first = primaryHits_F[f].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleYZ.second = primaryHits_B[b].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleXZ.first = primaryHits_F[f].primaryAngleXZ;
			tempCombinatorialTrack.primaryAngleXZ.second = primaryHits_B[b].primaryAngleXZ;


			tempCombinatorialTrack.combinatorialTrack_FB = true;
			tempCombinatorialTrack.combinatorialTrack_FT = false;
			tempCombinatorialTrack.combinatorialTrack_BT = false;

			allCombinatorialTracks.push_back(tempCombinatorialTrack);
			nCombinatorialTracks_FB++;
		    }
	    }

	for(unsigned int f = 0; f < primaryHits_F.size(); f++)
	    {
		for(unsigned int t = 0; t < primaryHits_T.size(); t++)
		    {
			combinatorialTracks tempCombinatorialTrack;
			tempCombinatorialTrack.tempId.first = primaryHits_F[f].tempId;
			tempCombinatorialTrack.tempId.second = primaryHits_T[t].tempId;

			tempCombinatorialTrack.hitPositionX.first = primaryHits_F[f].hitPositionX;
			tempCombinatorialTrack.hitPositionX.second = primaryHits_T[t].hitPositionX;
			tempCombinatorialTrack.hitPositionY.first = primaryHits_F[f].hitPositionY;
			tempCombinatorialTrack.hitPositionY.second = primaryHits_T[t].hitPositionY;
			tempCombinatorialTrack.hitPositionZ.first = primaryHits_F[f].hitPositionZ;
			tempCombinatorialTrack.hitPositionZ.second = primaryHits_T[t].hitPositionZ;

			tempCombinatorialTrack.primaryTrackId.first = primaryHits_F[f].primaryTrackId;
			tempCombinatorialTrack.primaryTrackId.second = primaryHits_T[t].primaryTrackId;
			tempCombinatorialTrack.primaryEnergy.first = primaryHits_F[f].primaryEnergy;
			tempCombinatorialTrack.primaryEnergy.second = primaryHits_T[t].primaryEnergy;
			tempCombinatorialTrack.primaryAngleYZ.first = primaryHits_F[f].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleYZ.second = primaryHits_T[t].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleXZ.first = primaryHits_F[f].primaryAngleXZ;
			tempCombinatorialTrack.primaryAngleXZ.second = primaryHits_T[t].primaryAngleXZ;

			tempCombinatorialTrack.combinatorialTrack_FB = false;
			tempCombinatorialTrack.combinatorialTrack_FT = true;
			tempCombinatorialTrack.combinatorialTrack_BT = false;

			allCombinatorialTracks.push_back(tempCombinatorialTrack);
			nCombinatorialTracks_FT++;
		    }
	    }
	for(unsigned int b = 0; b < primaryHits_B.size(); b++)
	    {
		for(unsigned int t = 0; t < primaryHits_T.size(); t++)
		    {
			combinatorialTracks tempCombinatorialTrack;
			tempCombinatorialTrack.tempId.first = primaryHits_B[b].tempId;
			tempCombinatorialTrack.tempId.second = primaryHits_T[t].tempId;

			tempCombinatorialTrack.hitPositionX.first = primaryHits_B[b].hitPositionX;
			tempCombinatorialTrack.hitPositionX.second = primaryHits_T[t].hitPositionX;
			tempCombinatorialTrack.hitPositionY.first = primaryHits_B[b].hitPositionY;
			tempCombinatorialTrack.hitPositionY.second = primaryHits_T[t].hitPositionY;
			tempCombinatorialTrack.hitPositionZ.first = primaryHits_B[b].hitPositionZ;
			tempCombinatorialTrack.hitPositionZ.second = primaryHits_T[t].hitPositionZ;

			tempCombinatorialTrack.primaryTrackId.first = primaryHits_B[b].primaryTrackId;
			tempCombinatorialTrack.primaryTrackId.second = primaryHits_T[t].primaryTrackId;
			tempCombinatorialTrack.primaryEnergy.first = primaryHits_B[b].primaryEnergy;
			tempCombinatorialTrack.primaryEnergy.second = primaryHits_T[t].primaryEnergy;
			tempCombinatorialTrack.primaryAngleYZ.first = primaryHits_B[b].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleYZ.second = primaryHits_T[t].primaryAngleYZ;
			tempCombinatorialTrack.primaryAngleXZ.first = primaryHits_B[b].primaryAngleXZ;
			tempCombinatorialTrack.primaryAngleXZ.second = primaryHits_T[t].primaryAngleXZ;

			tempCombinatorialTrack.combinatorialTrack_FB = false;
			tempCombinatorialTrack.combinatorialTrack_FT = false;
			tempCombinatorialTrack.combinatorialTrack_BT = true;

			allCombinatorialTracks.push_back(tempCombinatorialTrack);
			nCombinatorialTracks_BT++;
		    }
	    }

	if(verbo)
	    {
		logFile << "Total combinatorial tracks: " << allCombinatorialTracks.size()  << endl;
	    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Cross-check Combinatorial tracks
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(unsigned int iCombinatorialTrack = 0; iCombinatorialTrack < allCombinatorialTracks.size(); iCombinatorialTrack++)
	    {

		double combinatorialTrackAngleYZ = setAngle(atan2(allCombinatorialTracks[iCombinatorialTrack].hitPositionY.first -
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionY.second,
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first -
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second));
		double combinatorialTrackAngleXZ = setAngle(atan2(allCombinatorialTracks[iCombinatorialTrack].hitPositionX.first -
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionX.second,
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first -
								  allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second));

		bool matchToPrimary = (allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.first == allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.second);

		if(allCombinatorialTracks[iCombinatorialTrack].combinatorialTrack_FB)
		    {
			hCombinatorialAngleYZ_FB->Fill(combinatorialTrackAngleYZ);
			hCombinatorialAngleXZ_FB->Fill(combinatorialTrackAngleXZ);

			if(matchToPrimary)
			    {
				hCombinatorialPrimaryAngleYZ_FB->Fill(combinatorialTrackAngleYZ);
				hCombinatorialPrimaryAngleXZ_FB->Fill(combinatorialTrackAngleXZ);
				if(verbo)
				    {
					logFile << "Combinatorial Primary Angle YZ Front-Back: " << combinatorialTrackAngleYZ << endl;
					logFile << "Combinatorial Primary Angle XZ Front-Back: " << combinatorialTrackAngleXZ << endl;
					logFile << "Combinatorial Primary Front-Back TrackId (Energy): " << allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.first <<
                                            " (" << allCombinatorialTracks[iCombinatorialTrack].primaryEnergy.first << ")" << endl;
				    }
			    }
		    }
		if(allCombinatorialTracks[iCombinatorialTrack].combinatorialTrack_FT)
		    {
			hCombinatorialAngleYZ_FT->Fill(combinatorialTrackAngleYZ);
			hCombinatorialAngleXZ_FT->Fill(combinatorialTrackAngleXZ);

			if(matchToPrimary)
			    {
				hCombinatorialPrimaryAngleYZ_FT->Fill(combinatorialTrackAngleYZ);
				hCombinatorialPrimaryAngleXZ_FT->Fill(combinatorialTrackAngleXZ);
				if(verbo)
				    {
					logFile << "Combinatorial Primary Angle YZ Front-Top: " << combinatorialTrackAngleYZ << endl;
					logFile << "Combinatorial Primary Angle XZ Front-Top: " << combinatorialTrackAngleXZ << endl;
					logFile << "Combinatorial Primary Front-Top TrackId (Energy): " << allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.first <<
                                            " (" << allCombinatorialTracks[iCombinatorialTrack].primaryEnergy.first << ")" << endl;
				    }
			    }
		    }

		if(allCombinatorialTracks[iCombinatorialTrack].combinatorialTrack_BT)
		    {
			hCombinatorialAngleYZ_FT->Fill(combinatorialTrackAngleYZ);
			hCombinatorialAngleXZ_FT->Fill(combinatorialTrackAngleXZ);

			if(matchToPrimary)
			    {
				hCombinatorialPrimaryAngleYZ_BT->Fill(combinatorialTrackAngleYZ);
				hCombinatorialPrimaryAngleXZ_BT->Fill(combinatorialTrackAngleXZ);
				if(verbo)
				    {
					logFile << "Combinatorial Primary Angle YZ Back-Top: " << combinatorialTrackAngleYZ << endl;
					logFile << "Combinatorial Primary Angle XZ Back-Top: " << combinatorialTrackAngleXZ << endl;
					logFile << "Combinatorial Primary Back-Top TrackId (Energy): " << allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.first <<
                                            " (" << allCombinatorialTracks[iCombinatorialTrack].primaryEnergy.first << ")" << endl;
				    }
			    }
		    }
	    }
	if(verbo)
	    {
		logFile << endl;
		logFile << "Done with combinatorial track info!!" << endl << endl << endl << endl;
		logFile << "Beginning  reconstructed track info!!" << endl << endl;;
	    }
   
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Reconstructed tracks
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        art::Handle< vector<recob::Track> > trackListHandle;
	vector<art::Ptr<recob::Track> > trackList;
	if(event.getByLabel(fTrackModuleLabel, trackListHandle))
	    {
		art::fill_ptr_vector(trackList, trackListHandle);
	    }
	int nTracksReco = trackList.size();
	if(verbo)
	    {
		logFile << endl;
		logFile << "Total number of combinatorial tracks: " << allCombinatorialTracks.size() << endl;
		logFile << "Total number of reco tracks: " << nTracksReco << endl;
	    }
    
	// Find the associations between tracks and hits
	art::FindManyP<recob::Hit> hitsFromTrack(trackListHandle, event, fTrackModuleLabel);
	vector<tracksPair> allTracksPair;
	for(int iRecoTrack = 0; iRecoTrack < nTracksReco; ++iRecoTrack)
	    {
		nTotalRecoTracks++;
		unordered_map<int, double> trkIDE;
		for (auto const & hit : hitsFromTrack.at(iRecoTrack))
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

		TVector3 trackStartPosition = trackList[iRecoTrack]->Vertex();
		double trackStartPositionX = trackStartPosition.X();
		double trackStartPositionY = trackStartPosition.Y();
		double trackStartPositionZ = trackStartPosition.Z();

		// Don't use the very end points of the tracks (trackList[iRecoTrack]->End()) in case of scatter or distortion
		int nTrajectoryPoints = trackList[iRecoTrack]->NumberTrajectoryPoints();
		int lastPoint = nTrajectoryPoints - 2;
		TVector3 trackEndPosition = trackList[iRecoTrack]->LocationAtPoint(lastPoint);
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

		nConsideredRecoTracks++;
		hConsideredRecoTrack_RecoEnergy->Fill(totalTrackEnergy);

		double trackSlopeYZ = (trackEndPositionY - trackStartPositionY) / (trackEndPositionZ - trackStartPositionZ);
		double yIntercept = trackStartPositionY - (trackSlopeYZ * trackStartPositionZ);
		double trackSlopeXZ = (trackEndPositionX - trackStartPositionX) / (trackEndPositionZ - trackStartPositionZ);
		double xIntercept = trackStartPositionX - (trackSlopeXZ * trackStartPositionZ);
		if(verbo)
		    {
			logFile << "Considered Reconstruced Track info: " << endl;
			logFile << "TrackId: " << bestTrackId << endl;
			logFile << "Length: " << trackList[iRecoTrack]->Length() << endl;
			logFile << "Energy: " << totalTrackEnergy << endl;
			logFile << "Start (x,y,z): " << trackStartPositionX << ", " << trackStartPositionY << ", " << trackStartPositionZ << endl;
			logFile << "End (x,y,z): " << trackEndPositionX << ", " << trackEndPositionY << ", " << trackEndPositionZ << endl;
			logFile << "trackSlopeYZ: " << trackSlopeYZ << endl;
			logFile << "trackSlopeXZ: " << trackSlopeXZ << endl << endl;
		    }


		if(primaryMatchedRecoTrack.find(bestTrackId) != primaryMatchedRecoTrack.end())
		    {
			nPrimaryMatchedRecoTracks++;
			hPrimaryMatchedRecoTrack_RecoEnergy->Fill(totalTrackEnergy);
			hPrimaryMatchedRecoTrack_TrueEnergy->Fill(primaryMatchedRecoTrack[bestTrackId]);
			if(verbo)
			    {
				logFile << "Primary Matched Reconstruced Track info: " << endl;
				logFile << "TrackId: " << bestTrackId << endl;
				logFile << "Length: " << trackList[iRecoTrack]->Length() << endl;
				logFile << "Energy: " << totalTrackEnergy << endl;
				logFile << "Start (x,y,z): " << trackStartPositionX << ", " << trackStartPositionY << ", " << trackStartPositionZ << endl;
				logFile << "End (x,y,z): " << trackEndPositionX << ", " << trackEndPositionY << ", " << trackEndPositionZ << endl;
				logFile << "trackSlopeYZ: " << trackSlopeYZ << endl;
				logFile << "trackSlopeXZ: " << trackSlopeXZ << endl << endl;
			    }
		    }

		for(unsigned int iCombinatorialTrack = 0; iCombinatorialTrack < allCombinatorialTracks.size(); iCombinatorialTrack++)
		    {
			int tempId1 = allCombinatorialTracks[iCombinatorialTrack].tempId.first;
			int tempId2 = allCombinatorialTracks[iCombinatorialTrack].tempId.second;

			double X1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionX.first;
			double X2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionX.second;
			double Y1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionY.first;
			double Y2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionY.second;
			double Z1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first;
			double Z2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second;


			if(((trackEndPositionZ > trackStartPositionZ) && (Z2 < Z1)) || ((trackEndPositionZ < trackStartPositionZ) && (Z2 > Z1)))
			    {
				tempId1 = allCombinatorialTracks[iCombinatorialTrack].tempId.second;
				tempId2 = allCombinatorialTracks[iCombinatorialTrack].tempId.first;

				X1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionX.second;
				X2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionX.first;
				Y1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionY.second;
				Y2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionY.first;
				Z1 = allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second;
				Z2 = allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first;
			    }
		   
                    double combinatorialSlopeXZ = (X2 - X1) / (Z2 - Z1);
                    double combinatorialSlopeYZ = (Y2 - Y1) / (Z2 - Z1);
                    if(((trackSlopeYZ < 0) && (combinatorialSlopeYZ > 0)) || ((trackSlopeYZ > 0) && (combinatorialSlopeYZ < 0)))
                        {
                            continue;
                        }

                    double deltaSlopeYZ = fabs(trackSlopeYZ - combinatorialSlopeYZ);
                    hDeltaSlopeYZ->Fill(deltaSlopeYZ);

                    if(deltaSlopeYZ > MAX_DELTA_SLOPE_YZ)
                        {
                            continue;
                        }

                    double predictedHitPositionY1 = (trackSlopeYZ * allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first) + yIntercept;
                    double deltaY1 = fabs(allCombinatorialTracks[iCombinatorialTrack].hitPositionY.first - predictedHitPositionY1);
                    hDeltaY1->Fill(deltaY1);

                    double predictedHitPositionY2 = (trackSlopeYZ * allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second) + yIntercept;
                    double deltaY2 = fabs(allCombinatorialTracks[iCombinatorialTrack].hitPositionY.second - predictedHitPositionY2);
                    hDeltaY2->Fill(deltaY2);

                    if((deltaY1 > MAX_DELTA_Y) || (deltaY2 > MAX_DELTA_Y))
                        {
                            continue;
                        }

                    if(((trackSlopeXZ < 0) && (combinatorialSlopeXZ > 0)) || ((trackSlopeXZ > 0) && (combinatorialSlopeXZ < 0)))
                        {
                            continue;
                        }

                    double deltaSlopeXZ = fabs(trackSlopeXZ - combinatorialSlopeXZ);
                    hDeltaSlopeXZ->Fill(deltaSlopeXZ);

                    if(deltaSlopeXZ > MAX_DELTA_SLOPE_XZ)
                        {
                            continue;
                        }

                    double predictedHitPositionX1 = (trackSlopeXZ * allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.first) + xIntercept;
                    double deltaX1 = fabs(allCombinatorialTracks[iCombinatorialTrack].hitPositionX.first - predictedHitPositionX1);

                    double predictedHitPositionX2 = (trackSlopeXZ * allCombinatorialTracks[iCombinatorialTrack].hitPositionZ.second) + xIntercept;
                    double deltaX2 = fabs(allCombinatorialTracks[iCombinatorialTrack].hitPositionX.second - predictedHitPositionX2);

                    double deltaX = fabs(deltaX1 - deltaX2);

                    tracksPair tempTracksPair;
                    tempTracksPair.tempId.first = tempId1;
                    tempTracksPair.tempId.second = tempId2;
                    tempTracksPair.index.first = iRecoTrack;
                    tempTracksPair.index.second = iCombinatorialTrack;
                    tempTracksPair.slopeYZ.first = trackSlopeYZ;
                    tempTracksPair.slopeYZ.second = combinatorialSlopeYZ;
                    tempTracksPair.slopeXZ.first = trackSlopeXZ;
                    tempTracksPair.slopeXZ.second = combinatorialSlopeXZ;

                    tempTracksPair.trackId[0] = bestTrackId;
                    tempTracksPair.trackId[1] = allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.first;
                    tempTracksPair.trackId[2] = allCombinatorialTracks[iCombinatorialTrack].primaryTrackId.second;
                    tempTracksPair.energy[0] = totalTrackEnergy;
                    tempTracksPair.energy[1] = allCombinatorialTracks[iCombinatorialTrack].primaryEnergy.first;
                    tempTracksPair.energy[2] = allCombinatorialTracks[iCombinatorialTrack].primaryEnergy.second;

                    tempTracksPair.deltaX = deltaX;

                    allTracksPair.push_back(tempTracksPair);
                }
        }//iRecoTrack

    //Sort pair by ascending order of deltaR
    sort(allTracksPair.begin(), allTracksPair.end(), sortPair());

    //Require 1 to 1 matching- and save as unique pair
    vector<tracksPair> allUniqueTracksPair;
    while ( allTracksPair.size())
        {
            allUniqueTracksPair.push_back(allTracksPair.front());
            allTracksPair.erase(remove_if(allTracksPair.begin(), allTracksPair.end(), removePairIndex(allTracksPair.front())),
                                allTracksPair.end());
        }

    //Access information
    if(verbo)
        {
            logFile << "Unique pair size: " << allUniqueTracksPair.size() << endl;
        }
    for(unsigned int u = 0; u < allUniqueTracksPair.size(); u++)
        {
            double deltaX = allUniqueTracksPair[u].deltaX;
            hDeltaX->Fill(deltaX);

            if(deltaX > MIN_DELTA_X)
                {
                    continue;
                }

            nAllCRTMatchedRecoTracks++;
            hAllCRTMatchedRecoTrack_RecoEnergy->Fill(allUniqueTracksPair[u].energy[0]);
            //hAllCRTMatchedRecoTrack_TrueEnergy->Fill();

            if((allUniqueTracksPair[u].trackId[0] == allUniqueTracksPair[u].trackId[1]) && (allUniqueTracksPair[u].trackId[1] == allUniqueTracksPair[u].trackId[2]))
                {
                    nGoodCRTMatchedRecoTracks++;
                    hGoodCRTMatchedRecoTrack_RecoEnergy->Fill(allUniqueTracksPair[u].energy[0]);
                    hGoodCRTMatchedRecoTrack_TrueEnergy->Fill(allUniqueTracksPair[u].energy[1]);

                    if(verbo)
                        {
                            logFile << "Matching info: " << endl;
                            logFile << "ID of track: " << allUniqueTracksPair[u].trackId[0] << endl;
                            logFile << "IDs hits: " << allUniqueTracksPair[u].trackId[1] << ", " << allUniqueTracksPair[u].trackId[2] << endl;
                            logFile << "Energy of track: " << allUniqueTracksPair[u].energy[0] << endl;
                            logFile << "Energy of hits: " << allUniqueTracksPair[u].energy[1] << ", " << allUniqueTracksPair[u].energy[2] << endl;
                            logFile << "SlopeYZ of track: " << allUniqueTracksPair[u].slopeYZ.first << endl;
                            logFile << "SlopeYZ of combinatorial: " << allUniqueTracksPair[u].slopeYZ.second << endl;
                            logFile << "SlopeXZ of track: " << allUniqueTracksPair[u].slopeXZ.first << endl;
                            logFile << "SlopeXZ of combinatorial: " << allUniqueTracksPair[u].slopeXZ.second << endl;
                            logFile << "deltaX: " << deltaX << endl << endl;
                        }
                }
        }
    }//analyze
}// namespace CRTMatching

