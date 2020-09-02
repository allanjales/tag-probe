#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include <iostream>

using namespace std;

//Holder signal, background and signal+background histograms
class PassingFailing{
private:
	int& method;
	const char*& particleName;
	const char*& canvasWatermark;
	const char*& directoryToSave;
	const char*& particleType;
	const char*& tagOrProbe;
	InvariantMass& ObjMass;

	//About histograms
	const char*& quantityName;
	const char*& xAxisName;
	const char*& quantityUnit;
	const char*& extendedQuantityName;

	int& 	nBins;
	int& 	decimals;
	double& xMin;
	double&	xMax;

	void createHistogram(TH1D* &histo, const char* histoName)
	{
		//Set parameters
		string hName 		= string(particleType) + string(passingOrFailing) + string(tagOrProbe) + string(particleName) + "_" + string(quantityName) + string(histoName);
		string hTitle 		= string(extendedQuantityName) + " (" + string(passingOrFailing) + " in " + string(particleType) + " " + string(tagOrProbe) + ")";
		string xAxisTitle 	= xAxisName;
		string yAxisTitleForm = "Events";

		//Add unit if has
		if (strcmp(quantityUnit, "") != 0)
			xAxisTitle += " (" + string(quantityUnit) + ")";

		//Change title is passing
		if (strcmp(passingOrFailing, "Passing") == 0)
			hTitle = string(extendedQuantityName) + " (" + string(particleType) + " " + string(tagOrProbe) + ")";


		//Variable bin for pT
		if (strcmp(quantityName, "Pt") == 0)
		{
			double xbins[10000];
			xbins[0] = .1;
			int nbins = 0;
			double binWidth = 0.3;
			for (int i = 1; xbins[i-1] < xMax+binWidth; i++)
			{
				xbins[i] = xbins[i-1] < 0.3 ? 0.3 : xbins[i-1] *(1+binWidth);
				nbins++;
			}

			histo = new TH1D(hName.data(), hTitle.data(), nbins, xbins);
		}

		//Variable bin for eta
		else if (strcmp(quantityName, "Eta") == 0)
		{
			/*
			double xbins [] = {-2.4,-2.,-1.6,-1.2,-0.8,-0.4, 0., 0.4, 0.8, 1.2, 1.6, 2., 2.4};
			int nbins = sizeof(xbins)/sizeof(*xbins) - 1;
			*/

			double xbins[10000];
			xbins[0] = .1;
			int nbins = 0;
			double binWidth = 0.2;

			//For positive
			for (int i = 1; xbins[i-1] < xMax+binWidth; i++)
			{
				xbins[i] = xbins[i-1] < 0.2 ? 0.2 : xbins[i-1] *(1+binWidth);
				nbins++;
			}

			//Duplicate array and create another
			double rxbins[nbins*2+1];
			int entry = 0;
			for (int i = nbins; i >= 0; i--)
			{
				rxbins[entry] = -xbins[i];
				entry++;
			}
			rxbins[entry] = 0.;
			entry++;
			for (int i = 0; i <= nbins; i++)
			{
				rxbins[entry] = xbins[i];
				entry++;
			}

			/*
			//DEBUG
			cout << "---------------------------------" << endl;
			for (int i = 0; i < entry; i++)
			{
				cout << "[" << i << "] = " << rxbins[i] << endl;
			}
			*/


			
			histo = new TH1D(hName.data(), hTitle.data(), entry-1, rxbins);
		}

		//Histograms with no variable bin
		else
		{
			if (strcmp(quantityUnit, "") == 0)
			{
				yAxisTitleForm += " / (%1." + to_string(decimals) + "f)";
			}
			else
			{
				yAxisTitleForm += " / (%1." + to_string(decimals) + "f " + string(quantityUnit) + ")";
			}

			histo = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		}

		//Edit histogram axis
		histo->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), histo->GetBinWidth(0)));
		histo->GetXaxis()->SetTitle(xAxisTitle.data());
	}

	//For consistencyDebugCout()
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

public:
	const char* passingOrFailing = NULL;

	TH1D* hSigBack  = NULL;
	TH1D* hSig 		= NULL;
	TH1D* hBack 	= NULL;

	MassValues* PassFailObj()
	{
		if (strcmp(passingOrFailing, "Passing") == 0)
			return &this->ObjMass.Pass;
			//return &this->ObjMass.All;

		if (strcmp(passingOrFailing, "All") == 0)
			return &this->ObjMass.Pass;
			//return &this->ObjMass.All;

		cerr << "Could not find PassFailObj in PassingFailing class: " << particleType << " " << tagOrProbe << " " << quantityName <<  " " << passingOrFailing << " ERROR" << endl;
		return NULL;
	}

	void subtractSigHistogram()
	{
		this->hSig->Add(this->hSigBack, 1.);
		this->hSig->Add(this->hBack, -PassFailObj()->subtractionFactor());

		//Now normalize yields (to adapt variable binning)
		for (int i = 1; i <= hSig->GetXaxis()->GetNbins(); i++)
		{
			hSig->SetBinContent(i, hSig->GetBinContent(i)/hSig->GetBinWidth(i));
			hSigBack->SetBinContent(i, hSigBack->GetBinContent(i)/hSigBack->GetBinWidth(i));
			hBack->SetBinContent(i, hBack->GetBinContent(i)/hBack->GetBinWidth(i));
		}
	}

	//Fill histogram
	void fillQuantitiesHistograms(double& quantity, double& InvariantMass, bool storeInSignalHistogram = false)
	{
		if (!storeInSignalHistogram)
		{
			if (PassFailObj()->isInSignalRegion(InvariantMass))
				this->hSigBack->Fill(quantity);

			if (PassFailObj()->isInSidebandRegion(InvariantMass))
				this->hBack->Fill(quantity);
		}
		else
		{
			this->hSig->Fill(quantity);
		}

	}

	TCanvas* createDividedCanvas(bool shouldWrite = false, bool shouldSavePNG = true)
	{
		string canvasName 	= string(particleName) + " " + string(passingOrFailing) + " in " + string(particleType) + " " + string(tagOrProbe) + " " + string(quantityName);
		string titleLeft 	= string(extendedQuantityName) + " (" + string(passingOrFailing) + " in " + string(particleType) + " " + string(tagOrProbe) + ")";
		string titleRight 	= string(extendedQuantityName) + " of Signal (" + string(passingOrFailing) + " in " + string(particleType) + " " + string(tagOrProbe) + ")";
		string saveAs 		= string(directoryToSave) + string(particleType) + "_" + string(tagOrProbe) + "_" + string(quantityName) + "_" + string(passingOrFailing) + ".png";

		if (strcmp(passingOrFailing, "Passing") == 0)
			titleRight = string(extendedQuantityName) + " of Signal (" + string(particleType) + " " + string(tagOrProbe) + ")";

		//Create canvas and divide it
		TCanvas* c1 = new TCanvas(canvasName.data(), titleLeft.data(), 1200, 600);
		c1->Divide(2,1);

		//Select canvas part and set margin
		c1->cd(1);
		c1->cd(1)->SetMargin(0.14, 0.03, 0.11, 0.07);

		if (strcmp(quantityName, "Pt") == 0)
		{
			c1->cd(1)->SetLogy();
			hSigBack->SetMaximum(10*hSigBack->GetMaximum());
		}
		else
		{
			hSigBack->SetMinimum(0);
			hSigBack->SetMaximum(1.2*hSigBack->GetMaximum());
		}

		//Draws Main histogram
		hSigBack->SetLineWidth(2);		//Line Width
		hSigBack->Draw();

		//Draws Background histogram
		hBack->SetLineColor(kBlue); 	//Line Color
		hBack->SetLineStyle(kDashed);	//Line Style
		hBack->SetLineWidth(2);			//Line Width
		hBack->Draw("same");

		//Get Y range of draw
		gPad->Update();
		double Ymax = gPad->GetFrame()->GetY2();

		//Not show frame with mean, std dev
		gStyle->SetOptStat(0);

		//Add legend
		TLegend* l1_1 = new TLegend(0.65,0.80,0.92,0.90);
		l1_1->SetTextSize(0.04);
		l1_1->AddEntry(hSigBack,	"Total",		"l");
		l1_1->AddEntry(hBack,		"Background",	"l");
		l1_1->Draw();
		
		//Draws text information
		TLatex* tx1_1 = new TLatex();
		tx1_1->SetTextSize(0.04);
		tx1_1->SetTextFont(42);
		tx1_1->SetNDC(kTRUE);
		if (strcmp(quantityName, "Pt") == 0)
		{
			tx1_1->SetTextAlign(12);	//Align left, center
			tx1_1->DrawLatex(0.48,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
			tx1_1->DrawLatex(0.48,0.45,Form("%g entries (background)",	hBack->GetEntries()));
			tx1_1->DrawLatex(0.25,0.87,Form(canvasWatermark, ""));
		}
		else
		{
			tx1_1->SetTextAlign(22);	//Align center, center
			tx1_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
			tx1_1->DrawLatex(0.55,0.45,Form("%g entries (background)",	hBack->GetEntries()));
			tx1_1->DrawLatex(0.32,0.87,Form(canvasWatermark, ""));
		}

		//Select canvas part and set margin
		c1->cd(2);
		c1->cd(2)->SetMargin(0.14, 0.03, 0.11, 0.07);

		if (strcmp(quantityName, "Pt") == 0)
		{
			c1->cd(2)->SetLogy();
			hSig->SetMaximum(10*hSig->GetMaximum());
		}
		else
		{
			hSig->SetMinimum(0);
			hSig->SetMaximum(1.2*hSig->GetMaximum());
		}

		//Same range as comparision and Draws
		hSig->SetTitle(titleRight.data());
		hSig->SetLineColor(kMagenta); 	//Line Color
		hSig->SetLineWidth(2);			//Line Width
		hSig->Draw("same");

		//Add legend
		TLegend* l1_2 = new TLegend(0.65,0.85,0.92,0.90);
		l1_2->SetTextSize(0.04);
		l1_2->AddEntry(hSig, "Signal","l");
		l1_2->Draw();

		//Draws text information
		TLatex* tx1_2 = new TLatex();
		tx1_2->SetTextSize(0.04);
		tx1_2->SetTextFont(42);
		tx1_2->SetNDC(kTRUE);
		if (strcmp(quantityName, "Pt") == 0)
		{
			tx1_2->SetTextAlign(12);	//Align left, center
			tx1_2->DrawLatex(0.48,0.50,Form("%g entries (signal)", hSig->GetEntries()));
			tx1_2->DrawLatex(0.25,0.87,Form(canvasWatermark, ""));
		}
		else
		{
			tx1_2->SetTextAlign(22);	//Align center, center
			tx1_2->DrawLatex(0.55,0.5,Form("%g entries (signal)", hSig->GetEntries()));
			tx1_2->DrawLatex(0.32,0.87,Form(canvasWatermark, ""));
		}

		//Not show frame with mean, std dev
		gStyle->SetOptStat(0);

		//Writes in file
		if (shouldWrite == true)
		{
			c1->Write();
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
		const int minLegendSpace = 21;

		//Set information what are shown
		string legend = "- #";
		legend += fillAfter(string(particleType),     ' ', 10);
		legend += " ";
		legend += fillAfter(string(tagOrProbe),       ' ', 5);
		legend += " ";
		legend += fillAfter(string(quantityName), 	  ' ', 3);
		legend += " ";
		legend += fillAfter(string(passingOrFailing), ' ', 7);
		legend += " ";
		legend += "= ";

		//Diference calculus
		double diff = hSigBack->GetEntries() - (hSig->GetEntries() + (*PassFailObj()).subtractionFactor()*hBack->GetEntries());

		const char* addSpace = "";
		if (diff >= 0.)
			addSpace = " ";

		//Show information
		cout << legend << fixed << addSpace << diff;
		cout << " (-factor: " << (*PassFailObj()).subtractionFactor() << ")\n";
	}

	void writeQuantitiesHistogramsOnFile(bool hSigBack, bool hSig, bool hBack)
	{
		if (hSigBack == true)
			this->hSigBack->Write();

		if (hSig == true)
			this->hSig->Write();

		if (hBack == true)
			this->hBack->Write();
	}



	PassingFailing(int& method,
		const char*& particleName,
		const char*& canvasWatermark,
		const char*& directoryToSave,
	 	const char*& particleType,
	 	InvariantMass& ObjMass,
	 	const char*& tagOrProbe,
		const char*  passingOrFailing,
		const char*& quantityName,
		const char*& xAxisName,
		const char*& quantityUnit,
		const char*& extendedQuantityName,
		double& 	 xMin,
		double& 	 xMax,
		int&    	 nBins,
		int&    	 decimals)
		  : method(method),
			particleName(particleName),
		    canvasWatermark(canvasWatermark),
		    directoryToSave(directoryToSave),
			particleType(particleType),
		    ObjMass(ObjMass),
			tagOrProbe(tagOrProbe),
			passingOrFailing(passingOrFailing),
			quantityName(quantityName),
			xAxisName(xAxisName),
			quantityUnit(quantityUnit),
			extendedQuantityName(extendedQuantityName),
			nBins(nBins),
			xMin(xMin),
			xMax(xMax),
			decimals(decimals)
	{
		this->createHistogram(this->hSigBack, "SigBack");
		this->createHistogram(this->hSig, 	  "Sig");
		this->createHistogram(this->hBack, 	  "Back");
	}
};