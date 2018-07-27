////////////////////////////////////////////////////////////////////////
// Class:        MatchingWithSingleCRT
// Module Type:  analyzer
// File:         MatchingWithSingleCRT_module.h
// Author:       Arbin Timilsina, arbint@bnl.gov
////////////////////////////////////////////////////////////////////////

#ifndef MatchingWithSingleCRT_Module_H
#define MatchingWithSingleCRT_Module_H

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"

#include "larcore/Geometry/Geometry.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"

const double RAD_TO_DEG = 180.0 / M_PI;

const double CRT_XYZ_RESOLUTION = 2.5; //cm
const double CRT_T_RESOLUTION = 16.0; //ns
//The smallest possible increase of time the clock model allows is resolution; hence for 62.5 MHz clock, resolution is 16 ns

// Selection cuts
const double MIN_TRACK_LENGTH_YZ = 40.0; //cm
const double MIN_TRACK_ENERGY = 0.1; //GeV
const double MAX_DELTA_Y = 40.0; //cm
const double MAX_DELTA_X = 40.0;//
const double MAX_DELTA_T = 0.05; //us
const double MIN_DELTA_T = -0.20; //us
const double MIN_PE = 100.0;

const int dY = 2; //cm

typedef struct
{
    int tempId;

    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;
    double hitT0;

    bool hitFront;
    bool hitBack;

    int primaryTrackId;
    std::string primaryOrigin;
    double primaryEnergy;
    double primaryAngleYZ;

} hits;

typedef struct
{
    std::pair<int, int> index;
    std::pair<int, int> trackId;
    std::pair<bool, bool> hitFB;
    std::pair<double, double> T0;
    std::pair<double, double> AngleYZ;
    std::pair<double, double> energy;
    std::pair<std::string, std::string> origin;

    double slopeYX;
    double slopeYZ;
    double startPositionArray[3];
    double endPositionArray[3];

    double deltaX;
} tracksPair;

struct sortByR
{
    bool operator()(const tracksPair &pair1, const tracksPair &pair2)
    {
        return (pair1.deltaX < pair2.deltaX);
    }
};

struct sortByIndex
{
    bool operator()(const tracksPair &pair1, const tracksPair &pair2)
    {
        return (pair1.index.first < pair2.index.first);
    }
};

struct removePairIndex
{
    const tracksPair tracksPair1;
removePairIndex(const tracksPair &tracksPair0)
: tracksPair1(tracksPair0)
    {}

    bool operator() (const tracksPair &tracksPair2)
    {
        return ((tracksPair1.index.first == tracksPair2.index.first) || (tracksPair1.index.second == tracksPair2.index.second));
    }
};

namespace CRTMatching
{
    class MatchingWithSingleCRT : public art::EDAnalyzer
	{
	public:
	    /////////////////////////////////////////////////////////
	    // Constructor and destructor
	    /////////////////////////////////////////////////////////
	    explicit MatchingWithSingleCRT(fhicl::ParameterSet const& pset);
	    virtual ~MatchingWithSingleCRT();

	    /////////////////////////////////////////////////////////
	    //Required functions
	    /////////////////////////////////////////////////////////
	    //Called once per event
	    void analyze(const art::Event& evt);

	    /////////////////////////////////////////////////////////
	    //Selected optional functions
	    /////////////////////////////////////////////////////////
	    // Called once at the beginning of the job
	    void beginJob();

	    // Called once at the end of the job
	    void endJob();

	    // Called in to read .fcl file parameters
	    void reconfigure(fhicl::ParameterSet const& pset);

	    double setAngle(double angle)
	    {
		if(angle < 0)
		    {
			angle += M_PI;
		    }
		angle *= RAD_TO_DEG;
		return angle;
	    }

	private:

	    // Log file
	    std::ofstream logFile;

	    // The parameters that will be read from the .fcl file
	    std::string fTrackModuleLabel;
	    std::string fFlashModuleLabel;
	    int fSelectedPDG;
	    int fCRTsConfiguration;

	    // Muon Counter geometry
	    geo::MuonCounter35Alg *geoMuonCounter;
	    std::vector< std::vector<double> > counterGeometry;

	    //BackTracker
	    art::ServiceHandle<cheat::BackTrackerService> backTracker;

	    //ParticleInventory
	    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

	    // Pointers to services
	    art::ServiceHandle<geo::Geometry> fGeom;

	    //For statistics
	    unsigned int nEvents;
	    unsigned int nPrimaryMuons[4], nPrimaryMuonsWithOneHit[4], nPrimaryMuonsWithTwoHits[4];
	    unsigned int nHitsCRT, nHitsCRT_F, nHitsCRT_B;
	    unsigned int nTotalRecoTracks;
	    unsigned int nConsideredRecoTracks[4];
	    unsigned int nPrimaryMatchedRecoTracksOneHit[4], nAllCRTMatchedRecoTracksOneHit[4], nGoodCRTMatchedRecoTracksOneHit[4];
	    unsigned int nPrimaryMatchedRecoTracksTwoHits[4], nAllCRTMatchedRecoTracksTwoHits[4], nGoodCRTMatchedRecoTracksTwoHits[4];

	    //For hits
	    std::vector<hits> hitsCRT;

	    bool verbo = false;

	    //Histograms
	    TH1D *hStatistics;

	    TH2D *hRecoTrackEnergy;
	    TH1D *hRecoTrackLengthYZ;
	    TH2D *hRecoTrackLengthYZVsEnergy;

	    TH2D *hHitsVsEnergy;

	    TH1D *hPE;
	    TH1D *hFlashPeakTime;

	    TH1D *hDeltaT;
	    TH1D *hDeltaY;
	    TH1D *hDeltaX;

	    TH2D *hRecoTracksCoverageYX[4];
	    TH2D *hRecoTracksCoverageYZ[4];

	    TH1D *hPrimaryMuonsWithOneHit[4];
	    TH1D *hPrimaryMatchedRecoTracksOneHit[4];
	    TH1D *hAllCRTMatchedRecoTracksOneHit[4];
	    TH1D *hGoodCRTMatchedRecoTracksOneHit[4];
	    TH2D *hCoverageYXOneHit[4];
	    TH2D *hCoverageYZOneHit[4];

	    TH1D *hPrimaryMuonsWithTwoHits[4];
	    TH1D *hPrimaryMatchedRecoTracksTwoHits[4];
	    TH1D *hAllCRTMatchedRecoTracksTwoHits[4];
	    TH1D *hGoodCRTMatchedRecoTracksTwoHits[4];
	    TH2D *hCoverageYXTwoHits[4];
	    TH2D *hCoverageYZTwoHits[4];
	}; // class MatchingWithSingleCRT

}// namespace CRTMatching

DEFINE_ART_MODULE(CRTMatching::MatchingWithSingleCRT)

#endif // MatchingWithSingleCRT_Module

