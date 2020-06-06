#include "InvariantMass.h"
#include "TagProbe.h"

//Holder for 2 TagProbe class
class PassingFailing{
private:
	int *method;
	const char **particleName;

protected:
	double subtractionFactor = 1.0;

public:
	const char *PassingOrFailing;

	InvariantMass Mass{this->method, &this->subtractionFactor, this->particleName, &this->PassingOrFailing};
	TagProbe Tag  {this->method, &this->subtractionFactor, this->particleName, &this->PassingOrFailing, "Tag"};
	TagProbe Probe{this->method, &this->subtractionFactor, this->particleName, &this->PassingOrFailing, "Probe"};

	void prepareMethod()
	{
		this->defineHistogramsTexts();
		this->defineHistogramsNumbers();
		this->createSigBackHistograms();
		this->createBackHistograms();
		this->Mass.createMassHistogram();
	}

	void defineHistogramsTexts()
	{
		this->Tag  .defineHistogramsTexts();
		this->Probe.defineHistogramsTexts();
	}

	void defineHistogramsNumbers()
	{
		this->Tag  .defineHistogramsNumbers();
		this->Probe.defineHistogramsNumbers();
		this->Mass.defineNumbers(240, 2.8, 3.4);
	}

	void fillSigBackHistograms(Double_t TagPtValue, Double_t TagEtaValue, Double_t TagPhiValue, Double_t ProbePtValue, Double_t ProbeEtaValue, Double_t ProbePhiValue)
	{
		this->Tag  .fillSigBackHistograms(TagPtValue, 	TagEtaValue, 	TagPhiValue);
		this->Probe.fillSigBackHistograms(ProbePtValue, ProbeEtaValue, ProbePhiValue);
	}

	void fillBackHistograms(Double_t TagPtValue, Double_t TagEtaValue, Double_t TagPhiValue, Double_t ProbePtValue, Double_t ProbeEtaValue, Double_t ProbePhiValue)
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

	void createDividedCanvas(bool shouldWrite = false, bool shouldSave = false)
	{
		this->Tag  .createDividedCanvas(shouldWrite, shouldSave);
		this->Probe.createDividedCanvas(shouldWrite, shouldSave);
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

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSave = true)
	{
		this->Tag  .createEfficiencyCanvas(shouldWrite, shouldSave);
		this->Probe.createEfficiencyCanvas(shouldWrite, shouldSave);
	}

	void fillHistograms(double InvariantMass, double TagMuon_Pt, double TagMuon_Eta, double TagMuon_Phi, double ProbeMuon_Pt, double ProbeMuon_Eta, double ProbeMuon_Phi)
	{
		//If is inside signal region
		if (this->Mass.isInSignalRegion(InvariantMass))
		{
			this->fillSigBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}

		//If is inside sideband region
		if (this->Mass.isInSidebandRegion(InvariantMass))
		{
			this->fillBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}

	void updateSelectorParameters()
	{
		this->Mass.updateMassParameters();
	}

	void debugCout()
	{
		this->Tag  .debugCout();
		this->Probe.debugCout();
	}

	PassingFailing(int *method, const char **particleName, const char *PassingOrFailing)
		: method(method), particleName(particleName), PassingOrFailing(PassingOrFailing)
	{}
};