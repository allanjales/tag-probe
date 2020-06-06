#include "PassingFailing.h"

//Holder for 3 PassingFailing class
class Particle
{
private:
	int method = 0;
	//if 1 -> sideband by histogram
	//if 2 -> sideband by fitting

	void prepareMethod()
	{
		this->Pass.prepareMethod();
		this->Fail.prepareMethod();
		this->Both.prepareMethod();
	}

public:
	const char *particleName = "Muon";

	PassingFailing Pass{&this->method, &this->particleName, "Passing"};
	PassingFailing Fail{&this->method, &this->particleName, "Failing"};
	PassingFailing Both{&this->method, &this->particleName, "All"};

	void setMethod(int method)
	{
		this->method = method;
		this->prepareMethod();
	}

	void doFit()
	{
		this->Pass.Mass.doFit();
		this->Fail.Mass.doFit();
		this->Both.Mass.doFit();
	}

	void updateSelectionParameters()
	{
		this->Pass.updateSelectorParameters();
		this->Fail.updateSelectorParameters();
		this->Both.updateSelectorParameters();
	}

	void subtractSigHistograms()
	{
		this->Pass.subtractSigHistograms();
		this->Fail.subtractSigHistograms();
		this->Both.subtractSigHistograms();
	}

	void massDebugCout()
	{
		this->Pass.Mass.debugCout();
		this->Fail.Mass.debugCout();
		this->Both.Mass.debugCout();
	}

	void consistencyDebugCout()
	{
		cout << endl;
		cout << "Checking histograms number inconsistency (should be 0)" << endl;
		this->Pass.debugCout();
		this->Fail.debugCout();
		this->Both.debugCout();
	}
};