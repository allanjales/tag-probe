#include "PtEtaPhi.h"

//Holder for 3 set of histograms for each quantity
class TagProbe{
private:
	int* method					  = NULL;
	double* subtractionFactor	  = NULL;
	const char** particleName	  = NULL;
	const char** passingOrFailing = NULL;

public:
	const char* tagOrProbe = NULL;

	PtEtaPhi Pt {this->method, this->subtractionFactor, this->particleName, this->passingOrFailing, &this->tagOrProbe};
	PtEtaPhi Eta{this->method, this->subtractionFactor, this->particleName, this->passingOrFailing, &this->tagOrProbe};
	PtEtaPhi Phi{this->method, this->subtractionFactor, this->particleName, this->passingOrFailing, &this->tagOrProbe};

	void defineDefaultHistogramsTexts()
	{
		this->Pt .defineTexts("Pt",  "p_{t}",	"GeV/c", 	"Transversal Momentum");
		this->Eta.defineTexts("Eta", "#eta", 	"", 		"Pseudorapidity");
		this->Phi.defineTexts("Phi", "#phi", 	"rad", 		"Azimuthal Angle");
	}

	void defineDefaultHistogramsNumbers()
	{
		this->Pt .defineNumbers(100,	 0., 	100., 1);
		this->Eta.defineNumbers(200, 	-2.5, 	2.5);
		this->Phi.defineNumbers(79, 	-3.15, 	3.15);
	}

	void fillSigBackHistograms(double PtValue, double EtaValue, double PhiValue)
	{
		this->Pt .hSigBack->Fill(PtValue);
		this->Eta.hSigBack->Fill(EtaValue);
		this->Phi.hSigBack->Fill(PhiValue);
	}

	void fillBackHistograms(double PtValue, double EtaValue, double PhiValue)
	{
		this->Pt .hBack->Fill(PtValue);
		this->Eta.hBack->Fill(EtaValue);
		this->Phi.hBack->Fill(PhiValue);
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

	void createDividedCanvas(bool shouldWrite = false, const char* directoryToSave = "../result/", bool shouldSavePNG = false)
	{
		this->Pt .createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Eta.createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Phi.createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
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

	void createEfficiencyCanvas(bool shouldWrite = false, const char* directoryToSave = "../result/", bool shouldSavePNG = false)
	{
		this->Pt .createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Eta.createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		this->Phi.createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
	}

	void debugCout()
	{
		this->Pt .debugCout();
		this->Eta.debugCout();
		this->Phi.debugCout();
	}

	TagProbe(int* method, double* subtractionFactor, const char** particleName, const char** passingOrFailing, const char* tagOrProbe)
		: method(method), subtractionFactor(subtractionFactor), particleName(particleName), passingOrFailing(passingOrFailing), tagOrProbe(tagOrProbe)
	{}
};