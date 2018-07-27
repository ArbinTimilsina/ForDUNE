namespace MCTruthInformation
{
    bool InReadoutWindow( const simb::MCParticle* particle)
    {
	art::ServiceHandle<geo::Geometry> geom;
	auto const* detectorClocksService = lar::providerFrom<detinfo::DetectorClocksService>();
	auto const* detectorPropertiesService = lar::providerFrom<detinfo::DetectorPropertiesService>();
	double xDriftVelocity = detectorPropertiesService->DriftVelocity() * 1e-3; //Got cm/us, converted to cm/ns
	double windowSize = detectorPropertiesService->NumberTimeSamples() * detectorClocksService->TPCClock().TickPeriod() * 1e3; // Got in us, converted to ns

	// Look at each MC hit and see if it falls in readout window
	bool beenInVolume = false;
	for(unsigned int iTrajectoryPoint = 0; iTrajectoryPoint < particle->NumberTrajectoryPoints(); ++iTrajectoryPoint)
	    {
		const TLorentzVector& tmpPosition = particle->Position(iTrajectoryPoint);
		double tmpPosArray[3] = {tmpPosition[0], tmpPosition[1], tmpPosition[2]};

		// Check if hit is in TPC...
		geo::TPCID tpcid = geom->FindTPCAtPosition(tmpPosArray);
		if (tpcid.isValid)
		    {
			// Check if hit is within drift window...
			const geo::CryostatGeo& cryo = geom->Cryostat(tpcid.Cryostat);
			const geo::TPCGeo& tpc  = cryo.TPC(tpcid.TPC);
			double xPositionInPlane = tpc.PlaneLocation(0)[0];
			double driftTimeCorrection = fabs( tmpPosition[0] - xPositionInPlane ) / xDriftVelocity;
			double timeAtPlane = particle->T() + driftTimeCorrection;

			if ( (timeAtPlane < detectorPropertiesService->TriggerOffset()) || (timeAtPlane > (detectorPropertiesService->TriggerOffset() + windowSize)))
			    {
				continue;
			    }

                    if(!beenInVolume )
                        {
                            beenInVolume = true;
                        }
                } // TPC.valid
        } // NumberTrajectoryPoints()
    return beenInVolume;
} // InReadoutWindow
}//namespace MCTruthInformation
