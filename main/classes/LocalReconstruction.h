#include "PassingFailing.h"

//Holder for 3 PassingFailing class
class LocalReconstruction
{
private:
	int *method;
	const char **particleName;

public:
	const char *particleReconstruction = NULL;

	PassingFailing Pass{this->method, this->particleName, &this->particleReconstruction, "Passing"};
	PassingFailing Fail{this->method, this->particleName, &this->particleReconstruction, "Failing"};
	PassingFailing All {this->method, this->particleName, &this->particleReconstruction, "All"};

	void prepareMethod()
	{
		this->Pass.prepareMethod();
		this->Fail.prepareMethod();
		this->All .prepareMethod();
	}

	void doFit()
	{
		this->Pass.Mass.doFit();
		this->Fail.Mass.doFit();
		this->All .Mass.doFit();
	}

	void updateSelectionParameters()
	{
		this->Pass.updateSelectorParameters();
		this->Fail.updateSelectorParameters();
		this->All .updateSelectorParameters();
	}

	void subtractSigHistograms()
	{
		this->Pass.subtractSigHistograms();
		this->Fail.subtractSigHistograms();
		this->All .subtractSigHistograms();
	}

	void writehMass()
	{
		this->Pass.Mass.hMass->Write();
		this->Fail.Mass.hMass->Write();
		this->All .Mass.hMass->Write();
	}

	void massDebugCout()
	{
		this->Pass.Mass.debugCout();
		this->Fail.Mass.debugCout();
		this->All .Mass.debugCout();
	}

	void consistencyDebugCout()
	{
		cout << endl;
		cout << "Checking histograms number inconsistency (should be 0)" << endl;
		this->Pass.debugCout();
		this->Fail.debugCout();
		this->All .debugCout();
	}

	LocalReconstruction(int *method, const char **particleName, const char *particleReconstruction)
		: method(method), particleName(particleName), particleReconstruction(particleReconstruction)
	{}
};