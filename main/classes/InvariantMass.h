#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "RooFitResult.h"

#include <iostream>

using namespace std;

#include "FitFunctions.h"
#include "GlobalChi2.h"


struct MassValues
{
	//Histogram and fit function
	TH1D* hMass 	    = NULL;
	TF1*  fitFunction   = NULL;
	TF1*  fitSignal     = NULL;
	TF1*  fitBackground = NULL;

	//Regions
	double sidebandRegion1_x1  = 0.;
	double sidebandRegion1_x2  = 0.;
	double signalRegion_x1     = 0.;
	double signalRegion_x2     = 0.;
	double sidebandRegion2_x1  = 0.;
	double sidebandRegion2_x2  = 0.;

	//For downgrade
	TFitResultPtr fitResult = 0;

	Double_t resultParameters[12];

	const char* const fittingParName[12] = {
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

	//--- For sideband subtraction ---

	bool isInSignalRegion(double InvariantMass)
	{
		if (InvariantMass >= this->signalRegion_x1 && InvariantMass <= this->signalRegion_x2)
			return true;

		return false;
	}

	bool isInSidebandRegion(double InvariantMass)
	{
		if ((InvariantMass >= this->sidebandRegion1_x1 && InvariantMass <= this->sidebandRegion1_x2) ||
			(InvariantMass >= this->sidebandRegion2_x1 && InvariantMass <= this->sidebandRegion2_x2))
			return true;

		return false;
	}

	double subtractionFactor()
	{
		//Simple method (for linear background)
		double signalRegion = abs(signalRegion_x2 - signalRegion_x1);
		double sidebandRegion = abs(sidebandRegion1_x2 - sidebandRegion1_x1) + abs(sidebandRegion2_x2 - sidebandRegion2_x1);

		//Using yield (advanced method)
		if (fitFunction != NULL)
		{		
			signalRegion    = fitBackground->Integral(signalRegion_x1,    signalRegion_x2)   /hMass->GetBinWidth(0);
			sidebandRegion  = fitBackground->Integral(sidebandRegion1_x1, sidebandRegion1_x2)/hMass->GetBinWidth(0);
			sidebandRegion += fitBackground->Integral(sidebandRegion2_x1, sidebandRegion2_x2)/hMass->GetBinWidth(0);
		}
		else
		{
			cerr << "WARNING: not using advanced method for subtraction factor calculation. Using method for linear background." << endl;
		}

		return signalRegion/sidebandRegion;
	}

	void doFit()
	{
		TF1* &f  	= fitFunction;
		TF1* &fs 	= fitSignal;
		TF1* &fb 	= fitBackground;

		double xMin = hMass->GetXaxis()->GetXmin();
		double xMax = hMass->GetXaxis()->GetXmax();

		//Temporary for test
		xMin = 2.9;
		xMax = 3.3;

		//Create fit Function
		f = new TF1("FitFunction", FitFunctions::Merged::InvariantMass, xMin, xMax, 12);
		f->SetNpx(1000);
		f->SetLineStyle(kSolid);
		f->SetLineColor(kBlue);
		f->SetLineWidth(3);

		//Rename parameters
		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
		for (int i = 0; i < arraySize; i++)
		{
			f->SetParName(i, fittingParName[i]);
		}
		
		//Values Signal GS
		f->SetParameter(0,	4098.2);
		f->SetParameter(1,	3.09);
		f->SetParameter(2,	0.020);

		//Values Signal CB
		f->SetParameter(3,	1.58);
		f->SetParameter(4,	1.54);
		f->SetParameter(5,	3.093);
		f->SetParameter(6,	0.032);
		f->SetParameter(7,	42022.27);

		//Values Background
		f->SetParameter(8,	-0.217);
		f->SetParameter(9,	1.915);
		f->SetParameter(10, 263.185);
		f->SetParameter(11,	0.061);

		//Fit function
		fitResult = hMass->Fit(f, "RNS", "", xMin, xMax);
		f->GetParameters(resultParameters);
		
		//Signal Fitting
		fs = new TF1("FitFunction_Signal", FitFunctions::Merged::Signal_InvariantMass, xMin, xMax, 8);
		fs->SetNpx(1000);							//Resolution of signal fit function
		fs->SetParameters(resultParameters);		//Get only signal part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style
		fs->SetLineWidth(3);						//Fit width

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Merged::Background_InvariantMass, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[8]);	//Get only background part
		fb->SetLineColor(kBlue); 					//Fit Color
		fb->SetLineStyle(kDashDotted);				//Fit style
		fb->SetLineWidth(3);						//Fit width
		
		/*
		//TEST
		cout << "Entries  (TH1): " << hMass->GetEntries() << endl;
		cout << "Integral (TH1): " << hMass->Integral(0, hMass->GetNbinsX()+1) << endl;
		cout << "Integral (TF1): " << f->Integral(xMin, xMax)/hMass->GetBinWidth(0) << endl;
		*/
		cout << "chi2/ndf = " << (this->fitResult)->Chi2()/(this->fitResult)->Ndf() << "\n";
	}

	TBox* createTBox(double Ymax, int index = 0, double Ymin = 0.)
	{
		//index = -1 -> left region
		//index = 0 -> signal region
		//index = 1 -> right region

		double x1, x2 = 0;

		switch(index)
		{
			case -1:
				x1 = sidebandRegion1_x1;
				x2 = sidebandRegion1_x2;
				break;
			case 0:
				x1 = signalRegion_x1;
				x2 = signalRegion_x2;
				break;
			case 1:
				x1 = sidebandRegion2_x1;
				x2 = sidebandRegion2_x2;
				break;
		}

		return new TBox(x1, Ymin, x2, Ymax);
	}
};

class InvariantMass{
private:
	int& method;
	const char*& ressonance;
	const char*& particleName;
	const char*& canvasWatermark;
	const char*& directoryToSave;
	const char*& particleType;

	const char* const fittingParName[12] = {
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

		bool shouldDrawAllFitFunctions = true;

		canvas->cd(quarter);
		canvas->cd(quarter)->SetMargin(0.14, 0.02, 0.09, 0.07);

		hMass->SetMarkerStyle(20);		//Set markers style
		hMass->SetMarkerColor(kBlack);	//Set markers colors
		hMass->SetLineColor(kBlack);	//Set errobars color
		hMass->Draw("ep");

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
			if (shouldDrawAllFitFunctions == true)
			{
				//Change the size of TLegend
				tl->SetY1(tl->GetY1() - tl->GetTextSize()*3);

				//Get parameters of fit
				double fitParameters[12];
				fFit->GetParameters(fitParameters);

				//Signal Gaus Fitting
				TF1* fitGaus = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
				fitGaus->SetNpx(1000);						//Resolution of signal fit function
				fitGaus->SetParameters(fitParameters);		//Get only signal part
				fitGaus->SetLineColor(kMagenta); 			//Fit Color
				fitGaus->SetLineStyle(kDashed);				//Fit Style
				fitGaus->SetLineWidth(3);					//Fit width
				fitGaus->Draw("same");
				for (int i = 0; i < 12; i++)
					fitGaus->SetParName(i, this->fittingParName[i]);
				tl->AddEntry(fitGaus, "Gaussian",    "l");
				
				//Signal CB Fitting
				TF1* fitCB = new TF1("FitFunction_CrystalBall", FitFunctions::Primary::CrystalBall, xMin, xMax, 5);
				fitCB->SetNpx(1000);						//Resolution of signal fit function
				fitCB->SetParameters(&fitParameters[3]);	//Get only signal part
				fitCB->SetLineColor(kOrange); 				//Fit Color
				fitCB->SetLineStyle(kDotted);				//Fit Style
				fitCB->SetLineWidth(3);						//Fit width
				fitCB->Draw("same");
				for (int i = 0; i < 12; i++)
					fitCB->SetParName(i, this->fittingParName[i+3]);
				tl->AddEntry(fitCB,	 "Crystal Ball", "l");

				//Background Fitting
				TF1* fitExp = new TF1("FitFunction_Background", FitFunctions::Merged::Background_InvariantMass, xMin, xMax, 4);
				fitExp->SetNpx(1000);						//Resolution of background fit function
				fitExp->SetParameters(&fitParameters[8]);	//Get only background part
				fitExp->SetLineColor(kBlue); 				//Fit Color
				fitExp->SetLineStyle(kDashDotted);			//Fit style
				fitExp->SetLineWidth(3);					//Fit width
				fitExp->Draw("same");
				for (int i = 0; i < 12; i++)
					fitExp->SetParName(i, this->fittingParName[i+8]);
				tl->AddEntry(fitExp, "Exp + Exp",	 "l");
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
			Pass.doFit();

			cout << endl;
			cout << "Fitting All in " << particleType << " " << particleName << "...\n";
			All.doFit();
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			//Upsilon fit should go here
		}
	}

	void updateMassValuesFor(MassValues* ObjMassValues, bool isAll = false)
	{
		double value = 0.;
		double fwhm  = 0.;

		if (this->method == 1)
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
			value     = signalFit->GetMaximumX();
			double x1 = signalFit->GetX(signalFit->GetMaximum()/2);
			double x2 = signalFit->GetX(signalFit->GetMaximum()/2, x1+0.0001, value + x1*3);
			fwhm      = x2 - x1;

			delete signalFit;
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
			ObjMassValues->sidebandRegion1_x1  = 8.50;
			ObjMassValues->sidebandRegion1_x2  = 9.00;
			ObjMassValues->signalRegion_x1     = 9.19;
			ObjMassValues->signalRegion_x2     = 9.70;
			ObjMassValues->sidebandRegion2_x1  = 10.6;		//Run2011
			//ObjMassValues->sidebandRegion2_x1  = 9.70;	//MC
			ObjMassValues->sidebandRegion2_x2  = 11.2;
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

			//NEW FOR AVOID ???
			this->xMin  = 2.9;
			this->xMax  = 3.3;
			this->nBins = 160;
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			this->xMin  = 8.5;
			this->xMax  = 11.4;
			this->nBins = 60;
		}

		this->createMassHistogram(Pass.hMass, "Passing");
		this->createMassHistogram(All. hMass, "All");
	}
};