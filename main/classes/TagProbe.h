#include "Histograms.h"

//Holder for 3 Histograms Quantities
class TagProbe{
private:
	int *method;
	double *subtractionFactor;
	const char **particleName;
	const char **PassingOrFailing;

public:
	const char *tagOrProbe			= NULL;

	Histograms Pt {this->method, this->subtractionFactor, this->particleName, this->PassingOrFailing, &this->tagOrProbe};
	Histograms Eta{this->method, this->subtractionFactor, this->particleName, this->PassingOrFailing, &this->tagOrProbe};
	Histograms Phi{this->method, this->subtractionFactor, this->particleName, this->PassingOrFailing, &this->tagOrProbe};

	void defineHistogramsTexts()
	{
		this->Pt .defineTexts("Pt",  "p_{t}",	"GeV/c", 	"Transversal Momentum");
		this->Eta.defineTexts("Eta", "#eta", 	"", 		"Pseudorapidity");
		this->Phi.defineTexts("Phi", "#phi", 	"rad", 		"Azimuthal Angle");
	}

	void defineHistogramsNumbers()
	{
		this->Pt .defineNumbers(100,	 0., 	100., 1);
		this->Eta.defineNumbers(200, 	-2.5, 	2.5);
		this->Phi.defineNumbers(79, 	-3.15, 	3.15);
	}

	void fillSigBackHistograms(Double_t PtValue, Double_t EtaValue, Double_t PhiValue)
	{
		this->Pt .hSigBack->Fill(PtValue);
		this->Eta.hSigBack->Fill(EtaValue);
		this->Phi.hSigBack->Fill(PhiValue);
	}

	void fillBackHistograms(Double_t PtValue, Double_t EtaValue, Double_t PhiValue)
	{
		this->Pt .hBack->Fill(PtValue);
		this->Eta.hBack->Fill(EtaValue);
		this->Phi.hBack->Fill(PhiValue);
	}

	void fillSigHistograms(Double_t PtValue, Double_t EtaValue, Double_t PhiValue)
	{
		this->Pt .hSig->Fill(PtValue);
		this->Eta.hSig->Fill(EtaValue);
		this->Phi.hSig->Fill(PhiValue);
	}

	void createSigBackHistograms()
	{
		this->Pt .createSigBackHistogram();
		this->Eta.createSigBackHistogram();
		this->Phi.createSigBackHistogram();
	}

	void createBackHistograms()
	{
		this->Pt .createBackHistogram();
		this->Eta.createBackHistogram();
		this->Phi.createBackHistogram();
	}
	
	void subtractSigHistograms()
	{
		this->Pt .subtractSigHistogram();
		this->Eta.subtractSigHistogram();
		this->Phi.subtractSigHistogram();
	}

	void createDividedCanvas(bool shouldWrite = false, bool shouldSave = false)
	{
		this->Pt .createDividedCanvas(shouldWrite, shouldSave);
		this->Eta.createDividedCanvas(shouldWrite, shouldSave);
		this->Phi.createDividedCanvas(shouldWrite, shouldSave);
	}

	void write(bool hSigBack, bool hSig, bool hBack)
	{
		if (hSigBack == true)
		{
			this->Pt .hSigBack->Write();
			this->Eta.hSigBack->Write();
			this->Phi.hSigBack->Write();
		}

		if (hSig == true)
		{
			this->Pt .hSig->Write();
			this->Eta.hSig->Write();
			this->Phi.hSig->Write();
		}

		if (hBack == true)
		{
			this->Pt .hBack->Write();
			this->Eta.hBack->Write();
			this->Phi.hBack->Write();
		}
	}

	void createEfficiencyPlot(bool shouldWrite = false)
	{
		this->Pt .createEfficiencyPlot(shouldWrite);
		this->Eta.createEfficiencyPlot(shouldWrite);
		this->Phi.createEfficiencyPlot(shouldWrite);
	}

	void createEfficiencyCanvas(bool shouldWrite = false, bool shouldSave = false)
	{
		this->Pt .createEfficiencyCanvas(shouldWrite, shouldSave);
		this->Eta.createEfficiencyCanvas(shouldWrite, shouldSave);
		this->Phi.createEfficiencyCanvas(shouldWrite, shouldSave);
	}

	void debugCout()
	{
		this->Pt .debugCout();
		this->Eta.debugCout();
		this->Phi.debugCout();
	}

	TagProbe(int *method, double *subtractionFactor, const char **particleName, const char **PassingOrFailing, const char *tagOrProbe)
		: method(method), subtractionFactor(subtractionFactor), particleName(particleName), PassingOrFailing(PassingOrFailing), tagOrProbe(tagOrProbe)
	{}
};