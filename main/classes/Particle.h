#include "PassingFailing.h"
#include "InvariantMass.h"

//Holder for 3 PassingFailing class
class Particle
{
private:
	int method = 0;
	//if 1 -> sideband by histogram
	//if 2 -> sideband by fitting

public:
	const char *particleName = "Muon";

	PassingFailing TrackerPass{&this->method, &this->particleName, "Tracker Passing"};
	PassingFailing All        {&this->method, &this->particleName, "Tracker All"};

	InvariantMass MassTracker{&this->method, &this->particleName, &this->TrackerPass, &this->All};

	void setMethod(int method)
	{
		this->method = method;
		this->TrackerPass.prepareMethod();
		this->All        .prepareMethod();

		this->MassTracker.defineNumbers(240, 2.8, 3.4);
		this->MassTracker.createAllMassHistograms();
	}

	void doFit()
	{
		this->MassTracker.doFit();
	}

	void updateSelectionParameters()
	{
		this->MassTracker.updateMassAll();
	}

	void subtractSigHistograms()
	{
		this->TrackerPass.subtractSigHistograms();
		this->All        .subtractSigHistograms();
	}

	void writeMassHistograms(bool writehPass, bool writehAll)
	{
		if (writehPass == true)
		{
			TrackerPass.hMass->Write();
		}
		
		if (writehAll == true)
		{
			All.hMass->Write();
		}
	}

	void consistencyDebugCout()
	{
		cout << endl;
		cout << "Checking histograms number inconsistency (should be 0)" << endl;
		this->TrackerPass.debugCout();
		this->All        .debugCout();
	}
};