#include "PtEtaPhi.h"

//Holder for 3 set of histograms for each quantity
class TagProbe{
private:
	int* method 				 = NULL;
	const char** particleName	 = NULL;
	const char** directoryToSave = NULL;
	const char** particleType    = NULL;

	InvariantMass* ObjMass = NULL;

public:
	const char* tagOrProbe = NULL;

	PtEtaPhi Pt  {this->method, this->particleName, this->directoryToSave, this->particleType, this->ObjMass, &this->tagOrProbe,
		"Pt",  "p_{t}", "GeV/c", "Transversal Momentum", 100,  0.00, 100.00, 1};
	PtEtaPhi Eta {this->method, this->particleName, this->directoryToSave, this->particleType, this->ObjMass, &this->tagOrProbe,
		"Eta", "#eta",  "", 	 "Pseudorapidity",       200, -2.50,   2.50};
	PtEtaPhi Phi {this->method, this->particleName, this->directoryToSave, this->particleType, this->ObjMass, &this->tagOrProbe,
		"Phi", "#phi",  "rad",   "Azimuthal Angle",       79, -3.15,   3.15};

	void subtractSigHistograms()
	{
		this->Pt .subtractSigHistograms();
		this->Eta.subtractSigHistograms();
		this->Phi.subtractSigHistograms();
	}

	void fillQuantitiesHistograms(double** quantities, double* InvariantMass, int* isPassing)
	{
		/*
		//Assign variables for easy visualization
		double &pt  = *quantities[0];
		double &eta = *quantities[1];
		double &phi = *quantities[2];
		*/

		this->Pt .fillQuantitiesHistograms(quantities[0], InvariantMass, isPassing);
		this->Eta.fillQuantitiesHistograms(quantities[1], InvariantMass, isPassing);
		this->Phi.fillQuantitiesHistograms(quantities[2], InvariantMass, isPassing);
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Pt .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		this->Eta.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		this->Phi.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	void consistencyDebugCout()
	{
		this->Pt .consistencyDebugCout();
		this->Eta.consistencyDebugCout();
		this->Phi.consistencyDebugCout();
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		this->Pt .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->Eta.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->Phi.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		this->Pt .createEfficiencyPlot(shouldWrite);
		this->Eta.createEfficiencyPlot(shouldWrite);
		this->Phi.createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Pt .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		this->Eta.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		this->Phi.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}



	TagProbe(int* method,
		const char** particleName,
		const char** directoryToSave,
	 	const char** particleType,
	 	InvariantMass* ObjMass,
	 	const char*  tagOrProbe)
		  : method(method),
		    particleName(particleName),
		    directoryToSave(directoryToSave),
		    particleType(particleType),
		    ObjMass(ObjMass),
		    tagOrProbe(tagOrProbe)
	{}
};