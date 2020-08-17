#include "PassingFailing.h"

//Holder for 3 set of histograms for each quantity
class PtEtaPhi{
private:
	int* method 			     = NULL;
	const char** particleName    = NULL;
	const char** directoryToSave = NULL;
	const char** particleType    = NULL;
	const char** tagOrProbe      = NULL;

	InvariantMass* ObjMass = NULL;

public:
	//About histograms
	const char* quantityName 		 = NULL;
	const char* xAxisName			 = NULL;
	const char* quantityUnit		 = NULL;
	const char* extendedQuantityName = NULL;

	int 	nBins;
	double 	xMin;
	double	xMax;
	int 	decimals = 3;

	TEfficiency* pEff 	= NULL;

	PassingFailing Pass {this->method, this->particleName, this->directoryToSave, this->particleType, this->ObjMass, this->tagOrProbe,
		"Passing", &this->quantityName, &this->xAxisName, &this->quantityUnit, &this->extendedQuantityName,
		&this->nBins, &this->xMin, &this->xMax, &this->decimals};
	PassingFailing All  {this->method, this->particleName, this->directoryToSave, this->particleType, this->ObjMass, this->tagOrProbe,
		"All", &this->quantityName, &this->xAxisName, &this->quantityUnit, &this->extendedQuantityName,
		&this->nBins, &this->xMin, &this->xMax, &this->decimals};

	void subtractSigHistograms()
	{
		this->Pass.subtractSigHistogram();
		this->All .subtractSigHistogram();
	}

	void fillQuantitiesHistograms(double* quantity, double* InvariantMass, int* isPassing)
	{
		if (*isPassing)
			this->Pass.fillQuantitiesHistograms(quantity, InvariantMass);
		this->All.fillQuantitiesHistograms(quantity, InvariantMass);
	}

	void createQuantitiesCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		this->Pass.createDividedCanvas(shouldWrite, shouldSavePNG);
		this->All .createDividedCanvas(shouldWrite, shouldSavePNG);
	}

	//Creates a efficiency plot with histograms
	TEfficiency* createEfficiencyPlot(bool shouldWrite = false)
	{
		//References
		TH1D* &hPass  = this->Pass.hSigBack;
		TH1D* &hTotal = this->All .hSigBack;

		string pName 	= string(*particleName) + " " + string(*particleType) + " " + string(*tagOrProbe) + " " + string(quantityName) + " Efficiency";
		string pTitle 	= string(extendedQuantityName) + " Efficiency (" + string(*particleType) + " " + string(*tagOrProbe) + ")";

		//Set Y axis title for efficiency plot
		hTotal->GetYaxis()->SetTitle("Efficiency");

		//Check if are valid and consistent histograms
		if(TEfficiency::CheckConsistency(*hPass, *hTotal))
		{
			//Fills histogram
			this->pEff = new TEfficiency(*hPass, *hTotal);
			this->pEff->SetName(pName.data());
		}

		//Set plot config
		this->pEff->SetTitle(pTitle.data());
		this->pEff->SetLineWidth(2);
		this->pEff->SetLineColor(kBlack);
		this->pEff->SetMarkerStyle(21);
		this->pEff->SetMarkerSize(0.5);
		this->pEff->SetMarkerColor(kBlack);

		//Writes in file
		if (shouldWrite == true)
		{
			this->pEff->Write("",TObject::kOverwrite);
		}

		return pEff;
	}

	//Creates canvas for efficiency plots
	TCanvas* createEfficiencyCanvas(bool shouldWrite = false, bool shouldSavePNG = false)
	{
		//Supress canvas
		gROOT->SetBatch(0);

		string canvasName 	= string(*particleName) + " " + string(*tagOrProbe) + " " + string(quantityName) + " " + string(*particleType) + " Efficiency" ;
		string canvasTitle 	= string(extendedQuantityName) + " Efficiency (" + string(*particleType) + " " + string(*tagOrProbe) + ")";
		string saveAs 		= string(*directoryToSave) + string("Efficiency_") + string(*tagOrProbe) + "_" + string(quantityName) + "_" + string(*particleType) +".png";

		//Draw on canvas
		TCanvas* c1 = new TCanvas(canvasName.data(), canvasTitle.data(), 800, 600);
		c1->SetRightMargin(0.05);
		this->pEff->Draw();
		gPad->Update();

		//Set range in y axis
		auto graph = this->pEff->GetPaintedGraph(); 
		graph->SetMinimum(0.0);
		graph->SetMaximum(1.2);
		gPad->Update();

		//Add legend
		TLegend* l = new TLegend(0.75,0.82,0.92,0.88);
		l->SetTextSize(0.04);
		l->AddEntry(pEff, "Data", "lp");
		l->Draw();

		//CMS Open Data
		TLatex* txCOD = new TLatex();
		txCOD->SetTextSize(0.04);
		txCOD->SetTextAlign(12);
		txCOD->SetTextFont(42);
		txCOD->SetNDC(kTRUE);
		txCOD->DrawLatex(0.14,0.85,Form("#bf{CMS Open Data}"));

		//Set x range
		if (strcmp(quantityName, "Pt") == 0)
		{
			this->pEff->GetPaintedGraph()->GetHistogram()->GetXaxis()->SetRange(0.,40.);
		}

		//Writes in file
		if (shouldWrite == true)
		{
			c1->Write("",TObject::kOverwrite);
		}

		//If should save
		if (shouldSavePNG == true)
		{
			//Saves as image
			c1->SaveAs(saveAs.data());
		} 

		return c1;
	}

	void consistencyDebugCout()
	{
		this->Pass.consistencyDebugCout();
		this->All .consistencyDebugCout();
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		this->Pass .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
		this->All  .writeQuantitiesHistogramsOnFile(hSigBack, hSig, hBack);
	}



	PtEtaPhi(int* method,
		const char** particleName,
		const char** directoryToSave,
	 	const char** particleType,
	 	InvariantMass* ObjMass,
	 	const char** tagOrProbe,
		const char*  quantityName,
		const char*  xAxisName,
		const char*  quantityUnit,
		const char*  extendedQuantityName,
		int	    	 nBins,
		double	 	 xMin,
		double	 	 xMax,
		int	    	 decimals = 3)
		  : method(method),
		    particleName(particleName),
		    directoryToSave(directoryToSave),
		    particleType(particleType),
		    ObjMass(ObjMass),
		    tagOrProbe(tagOrProbe),
			quantityName(quantityName),
			xAxisName(xAxisName),
			quantityUnit(quantityUnit),
			extendedQuantityName(extendedQuantityName),
			nBins(nBins),
			xMin(xMin),
			xMax(xMax),
			decimals(decimals)
	{}
};