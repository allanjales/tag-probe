#include "Type.h"

//Holder for 3 Type class
class TagAndProbe
{
public:
	int method = 1;	// 1 | 2
	const char* ressonance = "Jpsi"; // "Jpsi" | "Upsilon"
	const char* particleName = "Muon";
	const char* directoryToSave = "../result/";

	Type Tracker    {this->method, this->ressonance, this->particleName, this->directoryToSave, "Tracker"};
	Type Standalone {this->method, this->ressonance, this->particleName, this->directoryToSave, "Standalone"};
	Type Global     {this->method, this->ressonance, this->particleName, this->directoryToSave, "Global"};
	
	void defineMassHistogramNumbers(double xMin, double xMax, int nBins, int decimals = 3)
	{
		this->Tracker   .defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
		this->Standalone.defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
		this->Global    .defineMassHistogramNumbers(xMin, xMax, nBins, decimals);
	}

	void doFit()
	{
		this->Tracker   .doFit();
		this->Standalone.doFit();
		this->Global    .doFit();
	}

	void updateMassValuesAll()
	{
		this->Tracker   .updateMassValuesAll();
		this->Standalone.updateMassValuesAll();
		this->Global    .updateMassValuesAll();
	}

	void createMassCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Tracker   .createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
		this->Standalone.createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
		this->Global    .createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}

	void subtractSigHistograms()
	{
		this->Tracker   .subtractSigHistograms();
		this->Standalone.subtractSigHistograms();
		this->Global    .subtractSigHistograms();
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Tracker   .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		this->Standalone.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
		this->Global    .createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	void consistencyDebugCout()
	{
		cout << "\nChecking histograms number inconsistency (should be 0)" << endl;
		cout << "* total - (background + signal)" << endl;
		this->Tracker   .consistencyDebugCout();
		this->Standalone.consistencyDebugCout();
		this->Global    .consistencyDebugCout();
	}

	void writeMassHistogramsOnFile(bool writehPass, bool writehAll)
	{
		this->Tracker   .writeMassHistogramsOnFile(writehPass, writehAll);
		this->Standalone.writeMassHistogramsOnFile(writehPass, writehAll);
		this->Global    .writeMassHistogramsOnFile(writehPass, writehAll);
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		this->Tracker   .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->Standalone.writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->Global    .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		this->Tracker   .createEfficiencyPlot(shouldWrite);
		this->Standalone.createEfficiencyPlot(shouldWrite);
		this->Global    .createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Tracker   .createEfficiencyCanvas(shouldWrite, shouldSavePNG);
		this->Standalone.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
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

		this->Tracker   .fillMassHistograms(*quantities[6], *types[0]);
		this->Standalone.fillMassHistograms(*quantities[6], *types[1]);
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

		this->Tracker   .fillQuantitiesHistograms(quantities, *types[0]);
		this->Standalone.fillQuantitiesHistograms(quantities, *types[1]);
		this->Global    .fillQuantitiesHistograms(quantities, *types[2]);
	}
	


	TagAndProbe()
	{}

	TagAndProbe(const char* ressonance)
			: ressonance(ressonance)
	{}
};