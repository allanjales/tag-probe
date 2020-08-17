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
	const char* particleName = "Muon";

	PassingFailing Tracker    {&this->method, &this->particleName, "Tracker"};
	PassingFailing Standalone {&this->method, &this->particleName, "Standalone"};
	PassingFailing Global     {&this->method, &this->particleName, "Global"};
	PassingFailing All        {&this->method, &this->particleName, "All"};

	InvariantMass MassTracker    {&this->method, &this->particleName, &this->Tracker, 	 &this->All};
	InvariantMass MassStandalone {&this->method, &this->particleName, &this->Standalone, &this->All};
	InvariantMass MassGlobal     {&this->method, &this->particleName, &this->Global, 	 &this->All};

	void setMethod(int method)
	{
		this->method = method;
		this->Tracker    .prepareMethod();
		this->Standalone .prepareMethod();
		this->Global     .prepareMethod();
		this->All        .prepareMethod();

		this->MassTracker.defineNumbers(240, 2.8, 3.4);
		this->MassTracker.createPassingAndAllMassHistograms();

		this->MassStandalone.defineNumbers(240, 2.8, 3.4);
		this->MassStandalone.createPassingAndAllMassHistograms();

		this->MassGlobal.defineNumbers(240, 2.8, 3.4);
		this->MassGlobal.createPassingAndAllMassHistograms();
	}

	void doFit()
	{
		this->MassTracker   .doFit();
		this->MassStandalone.doFit();
		this->MassGlobal    .doFit();
	}

	void updateSelectionParameters()
	{
		this->MassTracker   .updateMassAll();
		this->MassStandalone.updateMassAll();
		this->MassGlobal    .updateMassAll();
	}

	void subtractSigHistograms()
	{
		this->Tracker   .subtractSigHistograms();
		this->Standalone.subtractSigHistograms();
		this->Global    .subtractSigHistograms();
		this->All       .subtractSigHistograms();
	}

	void writeMassHistograms(bool writehPass, bool writehAll)
	{
		if (writehPass == true)
		{
			Tracker   .hMass->Write();
			Standalone.hMass->Write();
			Global    .hMass->Write();
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
		this->Tracker   .debugCout();
		this->Standalone.debugCout();
		this->Global    .debugCout();
		this->All       .debugCout();
	}
};