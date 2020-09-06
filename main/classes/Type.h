#include "InvariantMass.h"
#include "TagProbe.h"

//Holder for 2 TagProbe class
class Type{
private:
	int& method;
	const char*& ressonance;
	const char*& particleName;
	const char*& canvasWatermark;
	const char*& directoryToSave;

public:
	const char* particleType = NULL;

	InvariantMass Mass  {method, ressonance, particleName, canvasWatermark, directoryToSave, particleType};
	TagProbe      Tag   {method, ressonance, particleName, canvasWatermark, directoryToSave, particleType, Mass, "Tag"};
	TagProbe      Probe {method, ressonance, particleName, canvasWatermark, directoryToSave, particleType, Mass, "Probe"};

	void defineMassHistogramNumbers(double xMin, double xMax, int nBins, int decimals = 3)
	{
		Mass.defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
	}
	
	void doFit()
	{
		Mass.doFit();
	}

	void updateMassValuesAll()
	{
		Mass.updateMassValuesAll();
	}

	void createMassCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSavePNG = false)
	{
		Mass.createCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}

	void writeMassHistogramsOnFile(bool writehPass, bool writehAll)
	{
		Mass.writeMassHistogramsOnFile(writehPass, writehAll);
	}

	void subtractSigHistograms()
	{
		Tag  .subtractSigHistograms();
		Probe.subtractSigHistograms();
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		Tag  .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		Probe.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	void consistencyDebugCout()
	{
		Tag  .consistencyDebugCout();
		Probe.consistencyDebugCout();
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		Tag  .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		Probe.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		Tag  .createEfficiencyPlot(shouldWrite);
		Probe.createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		Tag  .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		Probe.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}



	void fillMassHistograms(double& InvariantMass, int& isPassing)
	{
		Mass.fillMassHistograms(InvariantMass, isPassing);
	}

	void fillQuantitiesHistograms(double** quantities, int& isPassing, bool storeInSignalHistogram = false)
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

		Tag  .fillQuantitiesHistograms(&quantities[3], *quantities[6], isPassing, storeInSignalHistogram);
		Probe.fillQuantitiesHistograms(quantities,     *quantities[6], isPassing, storeInSignalHistogram);
	}


	Type(int& method,
		const char*& ressonance,
		const char*& particleName,
		const char*& canvasWatermark,
		const char*& directoryToSave,
	 	const char*  particleType)
		  : method(method),
		    ressonance(ressonance),
		    particleName(particleName),
		    canvasWatermark(canvasWatermark),
		    directoryToSave(directoryToSave),
		    particleType(particleType)
	{}
};