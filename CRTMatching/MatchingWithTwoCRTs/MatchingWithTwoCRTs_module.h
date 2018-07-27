////////////////////////////////////////////////////////////////////////
// Class:        MatchingWithTwoCRTs
// Module Type:  analyzer
// File:         MatchingWithTwoCRTs_module.h
// Author:       Arbin Timilsina, arbint@bnl.gov
////////////////////////////////////////////////////////////////////////

#ifndef MatchingWithTwoCRTs_Module_H
#define MatchingWithTwoCRTs_Module_H

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"

#include "larcore/Geometry/Geometry.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"

const double RAD_TO_DEG = 180.0 / M_PI;
const double CRT_RESOLUTION = 2.5; //cm

// Angle binning
const int BINS_ANGLE = 180;
const double MIN_BINS_ANGLE = 0.0;
const double MAX_BINS_ANGLE = 180.0;

// Selection cuts
const double MIN_TRACK_LENGTH_YZ = 40.0;
const double MIN_TRACK_ENERGY = 0.10;

const double MAX_DELTA_SLOPE_YZ = 0.1;
const double MAX_DELTA_Y = 25.0;
const double MAX_DELTA_SLOPE_XZ = 0.1;
const double MIN_DELTA_X = 35.0;

typedef struct
{
    int tempId;

    double hitPositionX;
    double hitPositionY;
    double hitPositionZ;

    int primaryTrackId;
    double primaryEnergy;
    double primaryAngleYZ;
    double primaryAngleXZ;
} hits;

typedef struct
{
    std::pair<int, int> tempId;

    std::pair<double, double> hitPositionX;
    std::pair<double, double> hitPositionY;
    std::pair<double, double> hitPositionZ;

    std::pair<int, int> primaryTrackId;
    std::pair<double, double> primaryEnergy;
    std::pair<double, double> primaryAngleYZ;
    std::pair<double, double> primaryAngleXZ;

    bool combinatorialTrack_FB;
    bool combinatorialTrack_FT;
    bool combinatorialTrack_BT;
} combinatorialTracks;

typedef struct
{
    std::pair<int, int> tempId;
    std::pair<int, int> index;
    std::pair<double, double> slopeYZ;
    std::pair<double, double> slopeXZ;

    int trackId[3];
    double energy[3];

    double deltaX;
} tracksPair;

struct sortPair
{
    bool operator()(const tracksPair &pair1, const tracksPair &pair2)
    {
        return (pair1.deltaX < pair2.deltaX);
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
        return (tracksPair1.index.first == tracksPair2.index.first ||
                tracksPair1.tempId.first == tracksPair2.tempId.first ||
                tracksPair1.tempId.second == tracksPair2.tempId.second);
    }
};

namespace CRTMatching
{
    class MatchingWithTwoCRTs : public art::EDAnalyzer
	{
	public:
	    /////////////////////////////////////////////////////////
	    // Constructor and destructor
	    /////////////////////////////////////////////////////////
	    explicit MatchingWithTwoCRTs(fhicl::ParameterSet const& pset);
	    virtual ~MatchingWithTwoCRTs();

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
	    int fSelectedPDG;

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
	    unsigned int nEvents, nPrimaryMuons, nPrimaryMuonsWithTwoHits;
	    unsigned int nPrimaryHits_F, nPrimaryHits_B, nPrimaryHits_T;
	    unsigned int nPrimaryHits_FB, nPrimaryHits_FT, nPrimaryHits_BT;
	    unsigned int nCombinatorialTracks_FB, nCombinatorialTracks_FT, nCombinatorialTracks_BT;
	    unsigned int nTotalRecoTracks, nConsideredRecoTracks, nPrimaryMatchedRecoTracks;
	    unsigned int nAllCRTMatchedRecoTracks, nGoodCRTMatchedRecoTracks;

	    //For hits
	    std::vector<hits> primaryHits_F;
	    std::vector<hits> primaryHits_B;
	    std::vector<hits> primaryHits_T;

	    bool verbo = false;

	    //Histograms
	    TH1D *hStatistics;

	    //////////////////////////////////////////
	    //For sanity check
	    TH1D *hPrimaryAngleYZ_FB;
	    TH1D *hPrimaryAngleYZ_FT;
	    TH1D *hPrimaryAngleYZ_BT;
	    TH1D *hPrimaryAngleXZ_FB;
	    TH1D *hPrimaryAngleXZ_FT;
	    TH1D *hPrimaryAngleXZ_BT;

	    TH1D *hCombinatorialAngleYZ_FB;
	    TH1D *hCombinatorialAngleYZ_FT;
	    TH1D *hCombinatorialAngleYZ_BT;
	    TH1D *hCombinatorialAngleXZ_FB;
	    TH1D *hCombinatorialAngleXZ_FT;
	    TH1D *hCombinatorialAngleXZ_BT;

	    TH1D *hCombinatorialPrimaryAngleYZ_FB;
	    TH1D *hCombinatorialPrimaryAngleYZ_FT;
	    TH1D *hCombinatorialPrimaryAngleYZ_BT;
	    TH1D *hCombinatorialPrimaryAngleXZ_FB;
	    TH1D *hCombinatorialPrimaryAngleXZ_FT;
	    TH1D *hCombinatorialPrimaryAngleXZ_BT;
	    //////////////////////////////////////////

	    TH1D *hRecoTrackLengthYZ;
	    TH2D *hRecoTrackLengthYZVsEnergy;

	    TH1D *hDeltaSlopeYZ;
	    TH1D *hDeltaSlopeXZ;

	    TH1D *hDeltaY1;
	    TH1D *hDeltaY2;
	    TH1D *hDeltaX;

	    TH1D *hConsideredRecoTrack_RecoEnergy;
	    TH1D *hPrimaryMatchedRecoTrack_RecoEnergy;
	    TH1D *hPrimaryMatchedRecoTrack_TrueEnergy;
	    TH1D *hAllCRTMatchedRecoTrack_RecoEnergy;
	    TH1D *hAllCRTMatchedRecoTrack_TrueEnergy;
	    TH1D *hGoodCRTMatchedRecoTrack_RecoEnergy;
	    TH1D *hGoodCRTMatchedRecoTrack_TrueEnergy;
	}; // class MatchingWithTwoCRTs

}// namespace CRTMatching

DEFINE_ART_MODULE(CRTMatching::MatchingWithTwoCRTs)

#endif // MatchingWithTwoCRTs_Module

