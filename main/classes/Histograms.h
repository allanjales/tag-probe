//Pt, Eta, Phi histograms
class Histograms{
private:
	int *method;
	double *subtractionFactor;
	const char **particleName;
	const char **PassingOrFailing;
	const char **tagOrProbe;

	void createHistogram(TH1D* &histo, const char *histoName)
	{
		//Define parameters
		string hName 		= string(*PassingOrFailing) + string(*tagOrProbe) + string(*particleName) + "_" + string(quantityName) + string(histoName);
		string hTitle 		= string(extendedQuantityName) + " (" + string(*PassingOrFailing) + " " + string(*tagOrProbe) + ")";
		string xAxisTitle 	= xAxisName;
		string yAxisTitleForm;
		if (strcmp(quantityUnit, "") == 0)
		{
			yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f)";
		}
		else
		{
			xAxisTitle += " (" + string(quantityUnit) + ")";
			yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f " + string(quantityUnit) + ")";
		}

		//Create histogram
		histo = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		histo->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), histo->GetBinWidth(0)));
		histo->GetXaxis()->SetTitle(xAxisTitle.data());
	}

public:
	const char *quantityName 			= NULL;
	const char *extendedQuantityName	= NULL;
	const char *quantityUnit			= NULL;
	const char *xAxisName				= NULL;

	Int_t 		nBins;
	int 		decimals = 3;
	Double_t 	xMin;
	Double_t	xMax;

	double signalRegionEnd 		= 0;
	double sidebandRegionEnd	= 0;

	TH1D *hSigBack 		= NULL;
	TH1D *hSig 			= NULL;
	TH1D *hBack 		= NULL;
	TEfficiency* pEff 	= NULL;

	void defineTexts(const char *quantityName, const char *xAxisName, const char *quantityUnit,const char *extendedQuantityName)
	{
		this->quantityName 			= quantityName;
		this->extendedQuantityName 	= extendedQuantityName;
		this->xAxisName 			= xAxisName;
		this->quantityUnit			= quantityUnit;
	}

	void defineNumbers(Int_t nBins, Double_t xMin, Double_t xMax, int decimals = 3)
	{
		this->nBins 	= nBins;
		this->xMin 		= xMin;
		this->xMax 		= xMax;
		this->decimals 	= decimals;
	}

	void createSigBackHistogram()
	{
		this->createHistogram(hSigBack, "SigBack");
	}

	void createBackHistogram()
	{
		this->createHistogram(hBack, "Back");
	}

	void subtractSigHistogram()
	{
		if (this->hSig == NULL)
		{
			this->createHistogram(hSig, "Sig");

			this->hSig->Add(hSigBack,1);
			this->hSig->Add(hBack,-*this->subtractionFactor);
		}
		else
		{
			cout << "WARNING! Sig Histogram already exists! Could not Subtract" << endl;
		}
	}

	TCanvas *createDividedCanvas(bool shouldWrite = false, bool shouldSave = true)
	{
		string canvasName 	= string(*PassingOrFailing) + string(*tagOrProbe) + string(*particleName) + "_" + string(quantityName);
		string titleLeft 	= string(extendedQuantityName) + " (" + string(*PassingOrFailing) + " " + string(*tagOrProbe) + ")";
		string titleRight 	= string(extendedQuantityName) + " of Signal (" + string(*PassingOrFailing) + " " + string(*tagOrProbe) + ")";
		string saveAs 		= "../result/" + string(quantityName) + string(*PassingOrFailing) + string(*tagOrProbe) + ".png";

		//Create canvas and divide it
		TCanvas *c1 = new TCanvas(canvasName.data(), titleLeft.data(), 1200, 600);
		c1->Divide(2,1);

		//Select canvas part and set margin
		c1->cd(1);
		c1->cd(1)->SetMargin(0.14, 0.03, 0.11, 0.07);

		//Draws Main histogram
		hSigBack->SetMinimum(0);
		hSigBack->SetLineWidth(2);		//Line Width
		hSigBack->Draw();

		//Draws Background histogram
		hBack->SetLineColor(kBlue); 	//Line Color
		hBack->SetLineStyle(kDashed);	//Line Style
		hBack->SetLineWidth(2);			//Line Width
		hBack->Draw("same");

		//Draws Signal histogram
		hSig->SetLineColor(kMagenta); 	//Line Color
		hSig->SetLineWidth(2);			//Line Width
		hSig->Draw("same");

		//Get Y range of draw
		gPad->Update();
		Double_t Ymax = gPad->GetFrame()->GetY2();

		/*
		//Set fill color
		hSig->SetFillColor(kMagenta);
		hSigBack->SetFillColor(kBlue);
		hBack->SetFillColor(kYellow);
		*/

		//Add legend
		TLegend *l1_1 = new TLegend(0.65,0.75,0.92,0.90);
		l1_1->SetTextSize(0.04);
		l1_1->AddEntry(hSigBack,	"All",			"lp");
		l1_1->AddEntry(hSig,		"Signal",		"l");
		l1_1->AddEntry(hBack,		"Background",	"l");
		l1_1->Draw();
		
		//Draws text information
		TLatex *tx1_1 = new TLatex();
		tx1_1->SetTextSize(0.04);
		tx1_1->SetTextFont(42);
		tx1_1->SetNDC(kTRUE);
		if (strcmp(quantityName, "Pt") == 0) 
		{
			tx1_1->SetTextAlign(12);	//Align left, center
			tx1_1->DrawLatex(0.48,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
			tx1_1->DrawLatex(0.48,0.45,Form("%g entries (signal)",		hSig->GetEntries()));
			tx1_1->DrawLatex(0.48,0.40,Form("%g entries (background)",	hBack->GetEntries()));
			tx1_1->DrawLatex(0.25,0.87,Form("#bf{CMS Open Data}"));
		}
		else
		{
			tx1_1->SetTextAlign(22);	//Align center, center
			tx1_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
			tx1_1->DrawLatex(0.55,0.45,Form("%g entries (signal)",		hSig->GetEntries()));
			tx1_1->DrawLatex(0.55,0.40,Form("%g entries (background)",	hBack->GetEntries()));
			tx1_1->DrawLatex(0.32,0.87,Form("#bf{CMS Open Data}"));
		}

		//Select canvas part and set margin
		c1->cd(2);
		c1->cd(2)->SetMargin(0.14, 0.03, 0.11, 0.07);

		//Same range as comparision and Draws
		hSig->SetMinimum(0);
   		//hSig->SetMaximum(Ymax);
		hSig->SetTitle(titleRight.data());
		hSig->Draw("same");

		//Add legend
		TLegend *l1_2 = new TLegend(0.65,0.85,0.92,0.90);
		l1_2->SetTextSize(0.04);
		l1_2->AddEntry(hSig, "Signal","l");
		l1_2->Draw();

		//Draws text information
		TLatex *tx1_2 = new TLatex();
		tx1_2->SetTextSize(0.04);
		tx1_2->SetTextFont(42);
		tx1_2->SetNDC(kTRUE);
		if (strcmp(quantityName, "Pt") == 0)
		{
			tx1_2->SetTextAlign(12);	//Align left, center
			tx1_2->DrawLatex(0.48,0.50,Form("%g entries (signal)", hSig->GetEntries()));
			tx1_2->DrawLatex(0.25,0.87,Form("#bf{CMS Open Data}"));
		}
		else
		{
			tx1_2->SetTextAlign(22);	//Align center, center
			tx1_2->DrawLatex(0.55,0.5,Form("%g entries (signal)", hSig->GetEntries()));
			tx1_2->DrawLatex(0.32,0.87,Form("#bf{CMS Open Data}"));
		}

		//Not show frame with mean, std dev
		gStyle->SetOptStat(0);

		//Writes in file
		if (shouldWrite == true)
		{
			c1->Write();
		}

		//If should save
		if (shouldSave == true)
		{
			//Saves as image
			c1->SaveAs(saveAs.data());
		} 

		return c1;
	}

	//Creates a efficiency plot with histograms
	TEfficiency *createEfficiencyPlot(bool shouldWrite = false)
	{
		//References
		TH1D* &hPass  = hSig;
		TH1D* &hTotal = hSigBack;

		string pName 	= string(*PassingOrFailing) + string(*tagOrProbe) + string(*particleName) + "_" + string(quantityName) + "Efficiency";
		string pTitle 	= string(extendedQuantityName) + " Efficiency (" + string(*PassingOrFailing) + " " + string(*tagOrProbe) + ")";

		//Set Y axis title for efficiency plot
		hTotal->GetYaxis()->SetTitle("Efficiency");

		//Check if are valid and consistent histograms
		if(TEfficiency::CheckConsistency(*hPass, *hTotal))
		{
			//Fills histogram
			pEff = new TEfficiency(*hPass, *hTotal);
		}

		//Set plot config
		pEff->SetTitle(pTitle.data());
		pEff->SetLineWidth(2);
		pEff->SetLineColor(kBlack);
		pEff->SetMarkerStyle(21);
		pEff->SetMarkerSize(0.5);
		pEff->SetMarkerColor(kBlack);

		//Writes in file
		if (shouldWrite == true)
		{
			pEff->Write("",TObject::kOverwrite);
		}

		return pEff;
	}

	//Creates canvas for efficiency plots
	TCanvas *createEfficiencyCanvas(bool shouldWrite = false, bool shouldSave = false)
	{
		string canvasName 	= string(*PassingOrFailing) + string(*tagOrProbe) + string(*particleName) + "_" + string(quantityName) + "Efficiency";
		string canvasTitle 	= string(extendedQuantityName) + " Efficiency (" + string(*PassingOrFailing) + " " + string(*tagOrProbe) + ")";
		string saveAs 		= "../result/" + string(quantityName) + string(*PassingOrFailing) + string(*tagOrProbe) + "_Efficiency.png";

		//Draw on canvas
		TCanvas *c1 = new TCanvas(canvasName.data(), canvasTitle.data(), 800, 600);
		c1->SetRightMargin(0.05);
		pEff->Draw();
		gPad->Update();

		//Set range in y axis
		auto graph = pEff->GetPaintedGraph(); 
		graph->SetMinimum(0.0);
		graph->SetMaximum(1.2);
		if (strcmp(quantityName, "Phi") == 0)
		{
			graph->SetMinimum(0.5);
			graph->SetMaximum(1.0);
		} 
		gPad->Update();

		//Add legend
		TLegend *l = new TLegend(0.75,0.82,0.92,0.88);
		l->SetTextSize(0.04);
		l->AddEntry(pEff, "Data", "lp");
		l->Draw();

		//CMS Open Data
		TLatex *txCOD = new TLatex();
		txCOD->SetTextSize(0.04);
		txCOD->SetTextAlign(12);
		txCOD->SetTextFont(42);
		txCOD->SetNDC(kTRUE);
		txCOD->DrawLatex(0.14,0.85,Form("#bf{CMS Open Data}"));

		//Set x range
		if (strcmp(quantityName, "Pt") == 0)
		{
			pEff->GetPaintedGraph()->GetHistogram()->GetXaxis()->SetRange(0.,40.);
		}

		//Writes in file
		if (shouldWrite == true)
		{
			c1->Write("",TObject::kOverwrite);
		}

		//If should save
		if (shouldSave == true)
		{
			//Saves as image
			c1->SaveAs(saveAs.data());
		} 

		return c1;
	}

	string fillAfter(string text, char fillWith, int targetLength)
	{
		//Store size of text
		int textLength = strlen(text.data());

		//Fill de text in the end
		if(textLength < targetLength)
		{
			text.append(targetLength - textLength, fillWith);
		}

		return text;
	}

	void debugCout()
	{
		const int minLegendSpace = 21;

		//Set information what are shown
		string legend = "- #";
		legend += fillAfter(string(*PassingOrFailing) + " " + string(*tagOrProbe), ' ', 15);
		legend += fillAfter(string(quantityName), 		' ', 4);
		legend += "= ";

		//Show information
		cout << legend << hSigBack->GetEntries() - hSig->GetEntries() - hBack->GetEntries() << endl;
	}

	Histograms(int *method, double *subtractionFactor, const char **particleName, const char **PassingOrFailing, const char **tagOrProbe)
		: method(method), subtractionFactor(subtractionFactor), particleName(particleName), PassingOrFailing(PassingOrFailing), tagOrProbe(tagOrProbe)
	{}
};