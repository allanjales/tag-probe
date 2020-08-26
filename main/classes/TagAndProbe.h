#include "Type.h"

//Holder for 3 Type class
class TagAndProbe
{
public:
	int method = 1;	// 1 | 2
	const char* ressonance = "Jpsi"; // "Jpsi" | "Upsilon"
	const char* particleName = "Muon";
	const char* directoryToSave = "../result/";

	//Variables for computing each type
	bool doTracker    = true;
	bool doStandalone = true;
	bool doGlobal     = true;

	Type Tracker    {this->method, this->ressonance, this->particleName, this->directoryToSave, "Tracker"};
	Type Standalone {this->method, this->ressonance, this->particleName, this->directoryToSave, "Standalone"};
	Type Global     {this->method, this->ressonance, this->particleName, this->directoryToSave, "Global"};
	
	void defineMassHistogramNumbers(double xMin, double xMax, int nBins, int decimals = 3)
	{
		if (doTracker)
			this->Tracker   .defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
		if (doStandalone)
			this->Standalone.defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
		if (doGlobal)
			this->Global    .defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
	}

	void doFit()
	{
		if (doTracker)
			this->Tracker   .doFit();
		if (doStandalone)
			this->Standalone.doFit();
		if (doGlobal)
			this->Global    .doFit();
	}

	void updateMassValuesAll()
	{
		if (doTracker)
			this->Tracker   .updateMassValuesAll();
		if (doStandalone)
			this->Standalone.updateMassValuesAll();
		if (doGlobal)
			this->Global    .updateMassValuesAll();
	}

	void createMassCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSavePNG = false)
	{
		if (doTracker)
			this->Tracker   .createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
		if (doStandalone)
			this->Standalone.createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
		if (doGlobal)
			this->Global    .createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}

	void subtractSigHistograms()
	{
		if (doTracker)
			this->Tracker   .subtractSigHistograms();
		if (doStandalone)
			this->Standalone.subtractSigHistograms();
		if (doGlobal)
			this->Global    .subtractSigHistograms();
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		if (doTracker)
			this->Tracker   .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		if (doStandalone)
			this->Standalone.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		if (doGlobal)
			this->Global    .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	void consistencyDebugCout()
	{
		cout << "\nChecking histograms number inconsistency (should be 0)" << endl;
		cout << "* total - (background + alpha*signal)" << endl;
		if (doTracker)
			this->Tracker   .consistencyDebugCout();
		if (doStandalone)
			this->Standalone.consistencyDebugCout();
		if (doGlobal)
			this->Global    .consistencyDebugCout();
	}

	void writeMassHistogramsOnFile(bool writehPass, bool writehAll)
	{
		if (doTracker)
			this->Tracker   .writeMassHistogramsOnFile(writehPass, writehAll);
		if (doStandalone)
			this->Standalone.writeMassHistogramsOnFile(writehPass, writehAll);
		if (doGlobal)
			this->Global    .writeMassHistogramsOnFile(writehPass, writehAll);
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		if (doTracker)
			this->Tracker   .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		if (doStandalone)
			this->Standalone.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		if (doGlobal)
			this->Global    .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		if (doTracker)
			this->Tracker   .createEfficiencyPlot(shouldWrite);
		if (doStandalone)
			this->Standalone.createEfficiencyPlot(shouldWrite);
		if (doGlobal)
			this->Global    .createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		if (doTracker)
			this->Tracker   .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		if (doStandalone)
			this->Standalone.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		if (doGlobal)
			this->Global    .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}



	void fillMassHistograms(double** quantities, int** types)
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
		int &PassingProbeTrackingMuon   = *types[0];
		int &PassingProbeStandAloneMuon = *types[1];
		int &PassingProbeGlobalMuon     = *types[2];
		*/

		if (doTracker)
			this->Tracker   .fillMassHistograms(*quantities[6], *types[0]);
		if (doStandalone)
			this->Standalone.fillMassHistograms(*quantities[6], *types[1]);
		if (doGlobal)
			this->Global    .fillMassHistograms(*quantities[6], *types[2]);
	}

	void fillQuantitiesHistograms(double** quantities, int** types)
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
		int &PassingProbeTrackingMuon   = *types[0];
		int &PassingProbeStandAloneMuon = *types[1];
		int &PassingProbeGlobalMuon     = *types[2];
		*/

		if (doTracker)
			this->Tracker   .fillQuantitiesHistograms(quantities, *types[0]);
		if (doStandalone)
			this->Standalone.fillQuantitiesHistograms(quantities, *types[1]);
		if (doGlobal)
			this->Global    .fillQuantitiesHistograms(quantities, *types[2]);
	}
	


	TagAndProbe()
	{}

	TagAndProbe(const char* ressonance)
			: ressonance(ressonance)
	{}
};