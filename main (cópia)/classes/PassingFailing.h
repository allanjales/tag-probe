#include "TagProbe.h"

//Holder for 2 TagProbe class
class PassingFailing{
private:
	int* method 				= NULL;
	const char** particleName 	= NULL;

public:
	const char* passingOrFailing = NULL;

	TH1D* hMass = NULL;

	//In sigmas
	double signalRegion 	= 3.0;
	double sidebandRegion	= 6.0;

	//This will be changed by InvariantMass object
	double M_JPSI = 0;
	double W_JPSI = 0;
	double subtractionFactor = 1.;

	TagProbe Tag  {this->method, &this->subtractionFactor, this->particleName, &this->passingOrFailing, "Tag"};
	TagProbe Probe{this->method, &this->subtractionFactor, this->particleName, &this->passingOrFailing, "Probe"};

	void prepareMethod()
	{
		this->defineDefaultHistogramsTexts();
		this->defineDefaultHistogramsNumbers();
		this->createSigBackHistograms();
		this->createBackHistograms();
	}

	void defineDefaultHistogramsTexts()
	{
		this->Tag  .defineDefaultHistogramsTexts();
		this->Probe.defineDefaultHistogramsTexts();
	}

	void defineDefaultHistogramsNumbers()
	{
		this->Tag  .defineDefaultHistogramsNumbers();
		this->Probe.defineDefaultHistogramsNumbers();
	}

	void fillSigBackHistograms(double TagPtValue, double TagEtaValue, double TagPhiValue, double ProbePtValue, double ProbeEtaValue, double ProbePhiValue)
	{
		this->Tag  .fillSigBackHistograms(TagPtValue, 	TagEtaValue, 	TagPhiValue);
		this->Probe.fillSigBackHistograms(ProbePtValue, ProbeEtaValue, ProbePhiValue);
	}

	void fillBackHistograms(double TagPtValue, double TagEtaValue, double TagPhiValue, double ProbePtValue, double ProbeEtaValue, double ProbePhiValue)
	{
		this->Tag  .fillBackHistograms(TagPtValue, 	TagEtaValue, 	TagPhiValue);
		this->Probe.fillBackHistograms(ProbePtValue, ProbeEtaValue, ProbePhiValue);
	}

	void createSigBackHistograms()
	{
		this->Tag  .createSigBackHistograms();
		this->Probe.createSigBackHistograms();
	}

	void createBackHistograms()
	{
		this->Tag  .createBackHistograms();
		this->Probe.createBackHistograms();
	}

	void subtractSigHistograms()
	{
		this->Tag  .subtractSigHistograms();
		this->Probe.subtractSigHistograms();
	}

	void createDividedCanvas(bool shouldWrite = false, const char* directoryToSave = "../result/", bool shouldSavePNG = false)
	{
		this->Tag  .createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Probe.createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
	}

	void write(bool hSigBack, bool hSig, bool hBack)
	{
		this->Tag  .write(hSigBack, hSig, hBack);
		this->Probe.write(hSigBack, hSig, hBack);
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		this->Tag  .createEfficiencyPlot(shouldWrite);
		this->Probe.createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, const char* directoryToSave = "../result/", bool shouldSavePNG = true)
	{
		this->Tag  .createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Probe.createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
	}

	TBox* createTBox(double Ymax, int index = 0)
	{
		//index = -1 -> left region
		//index = 0 -> signal region
		//index = 1 -> right region

		double dx1, dx2 = 0;

		switch(index)
		{
			case -1:
				dx1 = -sidebandRegion;
				dx2 = -signalRegion;
				break;
			case 0:
				dx1 = -signalRegion;
				dx2 = +signalRegion;
				break;
			case 1:
				dx1 = +sidebandRegion;
				dx2 = +signalRegion;
				break;
		}

		double x1 = M_JPSI + W_JPSI * dx1;
		double x2 = M_JPSI + W_JPSI * dx2;

		TBox* region = new TBox(x1, 0., x2, Ymax);

		return region;
	}

	bool isInSignalRegion(double InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) < W_JPSI * signalRegion)
			return true;

		return false;
	}

	bool isInSidebandRegion(double InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) > W_JPSI * signalRegion && fabs(InvariantMass - M_JPSI) < W_JPSI * sidebandRegion)
			return true;

		return false;
	}

	void fillHistograms(double InvariantMass, double TagMuon_Pt, double TagMuon_Eta, double TagMuon_Phi, double ProbeMuon_Pt, double ProbeMuon_Eta, double ProbeMuon_Phi)
	{
		//If is inside signal region
		if (this->isInSignalRegion(InvariantMass))
		{
			this->fillSigBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}

		//If is inside sideband region
		if (this->isInSidebandRegion(InvariantMass))
		{
			this->fillBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}

	void debugCout()
	{
		this->Tag  .debugCout();
		this->Probe.debugCout();
	}

	PassingFailing(int* method, const char** particleName, const char* passingOrFailing)
		: method(method), particleName(particleName), passingOrFailing(passingOrFailing)
	{}
};