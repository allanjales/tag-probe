#include "InvariantMass.h"
#include "TagProbe.h"

//Holder for 2 TagProbe class
class Type{
private:
	int& method;
	const char*& ressonance;
	const char*& particleName;
	const char*& directoryToSave;

public:
	const char* particleType = NULL;

	InvariantMass Mass  {this->method, this->ressonance, this->particleName, this->directoryToSave, this->particleType};
	TagProbe      Tag   {this->method, this->particleName, this->directoryToSave, this->particleType, this->Mass, "Tag"};
	TagProbe      Probe {this->method, this->particleName, this->directoryToSave, this->particleType, this->Mass, "Probe"};

	void defineMassHistogramNumbers(double xMin, double xMax, int nBins, int decimals = 3)
	{
		Mass.defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
	}
	
	void doFit()
	{
		this->Mass.doFit();
	}

	void updateMassValuesAll()
	{
		this->Mass.updateMassValuesAll();
	}

	void createMassCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Mass.createCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}

	void writeMassHistogramsOnFile(bool writehPass, bool writehAll)
	{
		this->Mass.writeMassHistogramsOnFile(writehPass, writehAll);
	}

	void subtractSigHistograms()
	{
		this->Tag  .subtractSigHistograms();
		this->Probe.subtractSigHistograms();
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Tag  .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		this->Probe.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	void consistencyDebugCout()
	{
		this->Tag  .consistencyDebugCout();
		this->Probe.consistencyDebugCout();
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		this->Tag  .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->Probe.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		this->Tag  .createEfficiencyPlot(shouldWrite);
		this->Probe.createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Tag  .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		this->Probe.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}



	void fillMassHistograms(double& InvariantMass, int& isPassing)
	{
		this->Mass.fillMassHistograms(InvariantMass, isPassing);
	}

	void fillQuantitiesHistograms(double** quantities, int& isPassing)
	{
		/*
		//Assign variables for easy visualization
		double &ProbeMuon_Pt            = *quantities[0];
		double &ProbeMuon_Eta           = *quantities[1];
		double &ProbeMuon_Phi           = *quantities[2];
		double &TagMuon_Pt              = *quantities[3];
		double &TagMuon_Eta             = *quantities[4];
		double &TagMuon_Phi             = *quantities[5];
		double &InvariantMass           = *quantities[6];
		*/

		this->Tag  .fillQuantitiesHistograms(&quantities[3], *quantities[6], isPassing);
		this->Probe.fillQuantitiesHistograms(quantities,     *quantities[6], isPassing);
	}


	Type(int& method,
		const char*& ressonance,
		const char*& particleName,
		const char*& directoryToSave,
	 	const char*  particleType)
		  : method(method),
		    ressonance(ressonance),
		    particleName(particleName),
		    directoryToSave(directoryToSave),
		    particleType(particleType)
	{}
};