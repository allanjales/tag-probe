#include "LocalReconstruction.h"

//Holder for 3 Type class
class Particle
{
private:
	int method = 0;
	//if 1 -> sideband by histogram
	//if 2 -> sideband by fitting

	void prepareMethod()
	{
		this->Tracker.prepareMethod();
	}

public:
	const char *particleName = "Muon";

	LocalReconstruction Tracker{&this->method, &this->particleName, "Tracker"};

	void setMethod(int method)
	{
		this->method = method;
		this->prepareMethod();
	}

	void doFit()
	{
		this->Tracker.doFit();
	}

	void updateSelectionParameters()
	{
		this->Tracker.updateSelectionParameters();
	}

	void subtractSigHistograms()
	{
		this->Tracker.subtractSigHistograms();
	}

	void writehMass()
	{
		this->Tracker.writehMass();
	}

	void massDebugCout()
	{
		this->Tracker.massDebugCout();
	}

	void consistencyDebugCout()
	{
		this->Tracker.consistencyDebugCout();
	}
};