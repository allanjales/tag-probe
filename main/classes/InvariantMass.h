#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"

#include "RooFitResult.h"

#include <iostream>

using namespace std;

#include "FitFunctions.h"
#include "MassValues.h"

class InvariantMass{
private:
	int& method;
	const char*& ressonance;
	const char*& particleName;
	const char*& canvasWatermark;
	const char*& directoryToSave;
	const char*& particleType;

	void createMassHistogram(TH1D* &hMass, const char* passingOrFailing)
	{
		string hName 			= string(passingOrFailing) + "_" + string(particleType) + "_" + string(particleName) + "_InvariantMass";
		string hTitle 			= "Invariant Mass (" + string(passingOrFailing) + " for " + string(particleType) + ")";
		string yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f GeV/c^{2})";

		//Change hTitle name
		if (strcmp(passingOrFailing, "Passing") == 0)
			hTitle = "Invariant Mass (" + string(particleType) + ")";

		if (strcmp(passingOrFailing, "Failing") == 0)
			hTitle = "Invariant Mass (non-" + string(particleType) + ")";

		if (strcmp(passingOrFailing, "All") == 0)
			hTitle = "Invariant Mass (All)";


		//Create histogram
		hMass = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		hMass->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), hMass->GetBinWidth(0)));
		hMass->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass (GeV/c^{2})");
	}

	void drawCanvasQuarter(TCanvas* &canvas, bool drawRegions, int quarter, MassValues* ObjMassValues, int color = kBlue)
	{
		TH1D* &hMass = ObjMassValues->hMass;
		TF1*  &fFit  = ObjMassValues->fitFunction;

		canvas->cd(quarter);
		canvas->cd(quarter)->SetMargin(0.14, 0.02, 0.09, 0.07);

		hMass->SetMarkerStyle(20);		//Set markers style
		hMass->SetMarkerColor(kBlack);	//Set markers colors
		hMass->SetLineColor(kBlack);	//Set errobars color
		hMass->Draw("ep");

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			hMass->SetMinimum(0.);
			hMass->SetMaximum(1.2*hMass->GetMaximum());
		}

		//Add legend
		TLegend* tl = new TLegend(0.70,0.86,0.96,0.92);
		tl->SetTextSize(0.04);
		tl->AddEntry(hMass, "Data", "lp");

		//Draw fit
		if (fFit != NULL)
		{
			//Change the size of TLegend
			tl->SetY1(tl->GetY1() - tl->GetTextSize()*1);

			fFit->SetNpx(1000);
			fFit->SetLineColor(color);
			fFit->SetLineStyle(kSolid);
			fFit->Draw("same");
			tl->AddEntry(fFit, "Total Fit", "l");

			//If is showing pass and fail fit
			if (strcmp(ressonance, "Jpsi") == 0)
			{
				const char* const fittingParName[] = {
						"Gaus(Sg) Height  ",
						"Gaus(Sg) Position",
						"Gaus(Sg) Sigma   ",

						"CB  (Sg) Alpha   ",
						"CB  (Sg) N       ",
						"CB  (Sg) Mean    ",
						"CB  (Sg) Sigma   ",
						"CB  (Sg) Yield   ",

						"Exp1(Bg) Height  ",
						"Exp1(Bg) Width   ",
						"Exp2(Bg) Height  ",
						"Exp2(Bg) Width   "
					};

				//Change the size of TLegend
				tl->SetY1(tl->GetY1() - tl->GetTextSize()*3);

				//Get parameters of fit
				int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
				double fitParameters[arraySize];
				fFit->GetParameters(fitParameters);

				//Signal Gaus Fitting
				TF1* fitGaus = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
				fitGaus->SetNpx(1000);						//Resolution of signal fit function
				fitGaus->SetParameters(fitParameters);		//Get only signal part
				fitGaus->SetLineColor(kMagenta); 			//Fit Color
				fitGaus->SetLineStyle(kDashed);				//Fit Style
				fitGaus->SetLineWidth(3);					//Fit width
				fitGaus->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitGaus->SetParName(i, fittingParName[i]);
				tl->AddEntry(fitGaus, "Gaussian",    "l");
				
				//Signal CB Fitting
				TF1* fitCB = new TF1("FitFunction_CrystalBall", FitFunctions::Primary::CrystalBall, xMin, xMax, 5);
				fitCB->SetNpx(1000);						//Resolution of signal fit function
				fitCB->SetParameters(&fitParameters[3]);	//Get only signal part
				fitCB->SetLineColor(kOrange); 				//Fit Color
				fitCB->SetLineStyle(kDotted);				//Fit Style
				fitCB->SetLineWidth(3);						//Fit width
				fitCB->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitCB->SetParName(i, fittingParName[i+3]);
				tl->AddEntry(fitCB,	 "Crystal Ball", "l");

				//Background Fitting
				TF1* fitExp = new TF1("FitFunction_Background", FitFunctions::Jpsi::Background_InvariantMass, xMin, xMax, 4);
				fitExp->SetNpx(1000);						//Resolution of background fit function
				fitExp->SetParameters(&fitParameters[8]);	//Get only background part
				fitExp->SetLineColor(kBlue); 				//Fit Color
				fitExp->SetLineStyle(kDashDotted);			//Fit style
				fitExp->SetLineWidth(3);					//Fit width
				fitExp->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitExp->SetParName(i, fittingParName[i+8]);
				tl->AddEntry(fitExp, "Exp + Exp",	 "l");
			}

			//If is showing pass and fail fit
			if (strcmp(ressonance, "Upsilon") == 0)
			{
				const char* const fittingParName[] = {
						"Gaus(1S) Height  ",
						"Gaus(1S) Position",
						"Gaus(1S) Sigma   ",

						"Gaus(2S) Height  ",
						"Gaus(2S) Position",
						"Gaus(2S) Sigma   ",

						"Gaus(3S) Height  ",
						"Gaus(3S) Position",
						"Gaus(3S) Sigma   ",

						"Pol3(Bg) Height  ",
						"Pol3(Bg) Width   ",
						"Pol3(Bg) Height  ",
						"Pol3(Bg) Width   "
					};

				//Change the size of TLegend
				tl->SetY1(tl->GetY1() - tl->GetTextSize()*4);

				//Get parameters of fit
				int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
				double fitParameters[arraySize];
				fFit->GetParameters(fitParameters);

				//Signal Gaus Fitting
				TF1* fitGaus1 = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
				fitGaus1->SetNpx(1000);
				fitGaus1->SetParameters(fitParameters);
				fitGaus1->SetLineColor(kMagenta);
				fitGaus1->SetLineStyle(kDashed);
				fitGaus1->SetLineWidth(3);
				fitGaus1->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitGaus1->SetParName(i, fittingParName[i]);
				tl->AddEntry(fitGaus1, "Gauss (1S)", "l");

				//Signal Gaus Fitting
				TF1* fitGaus2 = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
				fitGaus2->SetNpx(1000);
				fitGaus2->SetParameters(&fitParameters[3]);
				fitGaus2->SetLineColor(kOrange);
				fitGaus2->SetLineStyle(kDashed);
				fitGaus2->SetLineWidth(3);
				fitGaus2->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitGaus2->SetParName(i, fittingParName[i+3]);
				tl->AddEntry(fitGaus2, "Gauss (2S)", "l");

				//Signal Gaus Fitting
				TF1* fitGaus3 = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
				fitGaus3->SetNpx(1000);
				fitGaus3->SetParameters(&fitParameters[6]);
				fitGaus3->SetLineColor(kMagenta - 5);
				fitGaus3->SetLineStyle(kDashed);
				fitGaus3->SetLineWidth(3);
				fitGaus3->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitGaus3->SetParName(i, fittingParName[i+6]);
				tl->AddEntry(fitGaus3, "Gauss (3S)", "l");

				//Background Fitting
				TF1* fitPol = new TF1("FitFunction_Background", FitFunctions::Upsilon::Background_InvariantMass, xMin, xMax, 4);
				fitPol->SetNpx(1000);
				fitPol->SetParameters(&fitParameters[9]);
				fitPol->SetLineColor(kBlue);
				fitPol->SetLineStyle(kDashDotted);
				fitPol->SetLineWidth(3);
				fitPol->Draw("same");
				for (int i = 0; i < arraySize; i++)
					fitPol->SetParName(i, fittingParName[i+9]);
				tl->AddEntry(fitPol, "Pol 3", "l");
			}
		}

		//Draw regions
		if (drawRegions == true)
		{
			//Change the size of TLegend
			tl->SetY1(tl->GetY1() - tl->GetTextSize()*2);

			//Get Y range of draw
			gPad->Update();
			double Ymax = gPad->GetFrame()->GetY2();
			double Ymin = gPad->GetFrame()->GetY1();

			//Draw regions
			TBox* side1  = ObjMassValues->createTBox(Ymax, -1, Ymin);
			side1->SetFillColorAlpha(kRed, 0.35);
			side1->Draw();
			TBox* signal = ObjMassValues->createTBox(Ymax, 0,  Ymin);
			signal->SetFillColorAlpha(kGreen, 0.35);
			signal->Draw();
			TBox* side2  = ObjMassValues->createTBox(Ymax, 1,  Ymin);
			side2->SetFillColorAlpha(kRed, 0.35);
			side2->Draw();

			//Add on TLegend
			tl->AddEntry(signal, "Signal R.",   "f");
			tl->AddEntry(side1,  "Sideband R.", "f");
		}

		//Draws information
		TLatex* tx = new TLatex();
		tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);
		tx->DrawLatex(0.16,0.88,Form(canvasWatermark, ""));

		//Draws TLegend
		tl->Draw();
	}

public:
	MassValues Pass;
	MassValues All;

	double 	xMin;
	double	xMax;
	int 	nBins;
	int 	decimals = 3;

	void defineMassHistogramNumbers(double xMin, double xMax, int nBins, int decimals = 3)
	{
		this->xMin     = xMin;
		this->xMax     = xMax;
		this->nBins    = nBins;
		this->decimals = decimals;

		delete this->Pass.hMass;
		delete this->All. hMass;

		this->createMassHistogram(Pass.hMass, "Passing");
		this->createMassHistogram(All. hMass, "All");
	}

	void fillMassHistograms(double& InvariantMass, int& isPassing)
	{
		if (isPassing)
			this->Pass.hMass->Fill(InvariantMass);
		this->All.hMass->Fill(InvariantMass);
	}

	void doFit()
	{
		if (strcmp(ressonance, "Jpsi") == 0)
		{
			cout << endl;
			cout << "Fitting Passing in " << particleType << " " << particleName << "...\n";
			Pass.doFitJpsi();

			cout << endl;
			cout << "Fitting All in " << particleType << " " << particleName << "...\n";
			All.doFitJpsi();
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			cout << endl;
			cout << "Fitting Passing in " << particleType << " " << particleName << "...\n";
			Pass.doFitUpsilon();

			cout << endl;
			cout << "Fitting All in " << particleType << " " << particleName << "...\n";
			All.doFitUpsilon();
		}
	}

	void updateMassValuesFor(MassValues* ObjMassValues, bool isAll = false)
	{
		double value = 0.;
		double fwhm  = 0.;

		//Default: method == 1
		{
			//Get value and uncertain of signal by histogram
			TH1D* &hMass = ObjMassValues->hMass;
			int bin0 = hMass->GetMaximumBin();
			value    = hMass->GetBinCenter(bin0);
			int bin1 = hMass->FindFirstBinAbove(hMass->GetMaximum()/2);
			int bin2 = hMass->FindLastBinAbove(hMass->GetMaximum()/2);
			fwhm     = hMass->GetBinCenter(bin2) - hMass->GetBinCenter(bin1);
		}

		if (this->method == 2)
		{
			//Get value and uncertain of signal by fitting
			TF1* &signalFit = ObjMassValues->fitSignal;
			if (signalFit != NULL)
			{
				value     = signalFit->GetMaximumX();
				double x1 = signalFit->GetX(signalFit->GetMaximum()/2);
				double x2 = signalFit->GetX(signalFit->GetMaximum()/2, x1+0.0001, value + x1*3);
				fwhm      = x2 - x1;
			}
			else
			{
				cerr << "Could not find signal fit function. Detemining singal region by histogram\n";
			}
		}

		double sigma = fwhm/2.355;

		if (strcmp(ressonance, "Jpsi") == 0)
		{
			//Signal region = mass +- 3*sigma
			ObjMassValues->signalRegion_x1 = value - 3*sigma;
			ObjMassValues->signalRegion_x2 = value + 3*sigma;

			//Sideband regions
			ObjMassValues->sidebandRegion1_x1 = value - 7*sigma;
			ObjMassValues->sidebandRegion1_x2 = value - 4*sigma;
			ObjMassValues->sidebandRegion2_x1 = value + 4*sigma;
			ObjMassValues->sidebandRegion2_x2 = value + 7*sigma;
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			//Signal region = mass +- 3*sigma
			ObjMassValues->signalRegion_x1 = value - 3*sigma;
			ObjMassValues->signalRegion_x2 = value + 3*sigma;

			ObjMassValues->sidebandRegion1_x1  = xMin;
			ObjMassValues->sidebandRegion1_x2  = value - 4*sigma;
			ObjMassValues->sidebandRegion2_x1  = 10.6;
			ObjMassValues->sidebandRegion2_x2  = xMax;
		}
	}

	void updateMassValuesAll()
	{
		updateMassValuesFor(&this->Pass);
		updateMassValuesFor(&this->All);
	}

	TCanvas* createCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSavePNG = false)
	{
		string canvasName 	= "InvariantMass_" + string(particleType);
		string canvasTitle	= "Invariant Mass " + string(particleType);
		string saveAs 		= string(directoryToSave) + "InvariantMass_" + string(particleType) + ".png";

		if (drawRegions)
		{
			canvasName 	= "InvariantMass_" + string(particleType) + "_region";
			canvasTitle	= "Invariant Mass " + string(particleType) + " with Regions";
			saveAs 		= string(directoryToSave) + "InvariantMass_" + string(particleType) + "_region" + ".png";
		}

		//Create canvas
		gStyle->SetCanvasPreferGL(kTRUE);
		TCanvas* c1 = new TCanvas(canvasName.data(), canvasTitle.data(), 1200, 600);
		c1->Divide(2,1);

		this->drawCanvasQuarter(c1, drawRegions, 1, &this->Pass, kGreen);
		this->drawCanvasQuarter(c1, drawRegions, 2, &this->All,  kBlue);

		c1->cd(2);

		//Draws information
		TLatex* tx = new TLatex();
		//tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);
		//tx->DrawLatex(0.61,0.60,Form("#chi^{2}/ndf = %.3g",(this->fitResult).Chi2()/(this->fitResult).Ndf()));

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

	void writeMassHistogramsOnFile(bool writehPass, bool writehAll)
	{
		this->Pass.hMass->Write();
		this->All .hMass->Write();
	}
	


	InvariantMass(int& method,
		const char*& ressonance,
		const char*& particleName,
		const char*& canvasWatermark,
		const char*& directoryToSave,
	 	const char*& particleType)
		  : method(method),
		    ressonance(ressonance),
		    particleName(particleName),
		    canvasWatermark(canvasWatermark),
		    directoryToSave(directoryToSave),
		    particleType(particleType)
	{
		if (strcmp(ressonance, "Jpsi") == 0)
		{
			this->xMin  = 2.8;
			this->xMax  = 3.4;
			this->nBins = 240;

			//NEW FOR AVOID "COULD NOT FIT"
			this->xMin  = 2.9;
			this->xMax  = 3.3;
			this->nBins = 160;
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			this->xMin  = 8.7;
			this->xMax  = 11.;
			this->nBins = 60;
		}

		this->createMassHistogram(Pass.hMass, "Passing");
		this->createMassHistogram(All. hMass, "All");
	}
};