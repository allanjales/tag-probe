#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"

#include <iostream>

using namespace std;

#include "FitFunctions.h"

//Store invariant mass class
class InvariantMass{
private:
	int *method;
	double *subtractionFactor;
	const char **particleName;
	const char **PassingOrFailing;

	//Bins per each x axis unit (for integrations)
	double scale = 0;

	TF1* fitFunction 		= NULL;
	TF1* fitFunctionSig 	= NULL;
	TF1* fitFunctionGaus 	= NULL;
	TF1* fitFunctionCB 		= NULL;
	TF1* fitFunctionBack 	= NULL;
	TFitResultPtr fitResult = 0;
	
	//References
	TF1* &f  	= fitFunction;
	TF1* &fs 	= fitFunctionSig;
	TF1* &fsG 	= fitFunctionGaus;
	TF1* &fsCB 	= fitFunctionCB;
	TF1* &fb 	= fitFunctionBack;

	//extracted from method
	double M_JPSI = 3.097;
	double W_JPSI = 0.010;
	double signalRegion 	= 3.0;
	double sidebandRegion	= 6.0;

public:
	int			nBins;
	int			decimals = 4;
	double 		xMin;
	double		xMax;
	int color = kBlue;

	TH1D* hMass = NULL;

	//Fitting
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

	Double_t resultParameters[12];

	void defineNumbers(int nBins, double xMin, double xMax, int decimals = 4)
	{
		this->nBins 	= nBins;
		this->xMin 		= xMin;
		this->xMax 		= xMax;
		this->decimals 	= decimals;
	}

	void createMassHistogram()
	{
		string hName 			= string(*PassingOrFailing) + string(*particleName) + "InvariantMass";
		string hTitle 			= "Invariant Mass (" + string(*PassingOrFailing) + ")";
		string yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f GeV/c^{2})";

		//Create histogram
		hMass = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		hMass->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), hMass->GetBinWidth(0)));
		hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");

		//Get scale from histogram
		this->scale = 1/hMass->GetBinWidth(0);
	}

	void fill(double value)
	{
		this->hMass->Fill(value);
	}

	bool isInSignalRegion(Double_t InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) < W_JPSI * signalRegion)
			return true;

		return false;
	}

	bool isInSidebandRegion(Double_t InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) > W_JPSI * signalRegion && fabs(InvariantMass - M_JPSI) < W_JPSI * sidebandRegion)
			return true;

		return false;
	}

	void doFit()
	{
		//Create fit Function
		f = new TF1("FitFunction", FitFunctions::Merged::FFit_InvariantMassAll, xMin, xMax, 12);
		f->SetNpx(1000);
		f->SetLineStyle(kSolid);
		f->SetLineColor(color);
		f->SetLineWidth(3);

		//Rename parameters
		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
		for (int i = 0; i < arraySize; i++)
		{
			f->SetParName(i, fittingParName[i]);
		}
		
		//Values Signal
		f->SetParameter(0,	340.2);
		f->SetParameter(1,	3.09);
		f->SetParameter(2,	0.037);

		f->SetParameter(3,	1.824);
		f->SetParameter(4,	1.034);
		f->SetParameter(5,	3.093);
		f->SetParameter(6,	0.022);
		f->SetParameter(7,	8322.27);

		//Values Background
		f->SetParameter(8,	-0.217);
		f->SetParameter(9,	1.915);
		f->SetParameter(10, 263.185);
		f->SetParameter(11,	0.061); 

		//Fit function
		cout << endl;
		cout << "Fitting " << *PassingOrFailing << " " << *particleName << "..." << endl;
		fitResult = hMass->Fit(f, "RNS", "", xMin, xMax);
		f->GetParameters(resultParameters);
		
		//Signal Fitting
		fs = new TF1("FitFunction_Signal", FitFunctions::Merged::Signal_InvariantMassAll, xMin, xMax, 8);
		fs->SetNpx(1000);							//Resolution of signal fit function
		fs->SetParameters(resultParameters);		//Get only signal part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style
		fs->SetLineWidth(3);						//Fit width
		
		//Signal Gaus Fitting
		fsG = new TF1("FitFunction_Gaussian", FitFunctions::Primary::Gaus, xMin, xMax, 3);
		fsG->SetNpx(1000);							//Resolution of signal fit function
		fsG->SetParameters(resultParameters);		//Get only signal part
		fsG->SetLineColor(kGreen);	 				//Fit Color
		fsG->SetLineStyle(kDashed);					//Fit Style
		fsG->SetLineWidth(3);						//Fit width
		
		//Signal CB Fitting
		fsCB = new TF1("FitFunction_CrystalBall", FitFunctions::Primary::CrystalBall, xMin, xMax, 5);
		fsCB->SetNpx(1000);							//Resolution of signal fit function
		fsCB->SetParameters(&resultParameters[3]);	//Get only signal part
		fsCB->SetLineColor(kRed); 					//Fit Color
		fsCB->SetLineStyle(kDotted);				//Fit Style
		fsCB->SetLineWidth(3);						//Fit width

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Merged::Background_InvariantMassAll, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[8]);	//Get only background part
		fb->SetLineColor(color); 					//Fit Color
		fb->SetLineStyle(kDashDotted);				//Fit style
		fb->SetLineWidth(3);						//Fit width
	}

	void updateMassParameters()
	{
		//Get value and uncertain of signal
		double center = 0;
		double fwhm = 0;

		if (*this->method == 1)
		{
			//Get value and uncertain of signal
			int bin0 = this->hMass->GetMaximumBin();
			center   = this->hMass->GetBinCenter(bin0);
			int bin1 = this->hMass->FindFirstBinAbove(this->hMass->GetMaximum()/2);
			int bin2 = this->hMass->FindLastBinAbove(this->hMass->GetMaximum()/2);
			fwhm     = this->hMass->GetBinCenter(bin2) - this->hMass->GetBinCenter(bin1);
		}

		if (*this->method == 2)
		{
			//Get value and uncertain of signal
			center    = this->fs->GetMaximumX();
			double x1 = this->fs->GetX(this->fs->GetMaximum()/2);
			double x2 = this->fs->GetX(this->fs->GetMaximum()/2, x1+0.0001, center + x1*3);
			fwhm      = x2 - x1;
		}

		double sigma = fwhm/2.355;

		this->M_JPSI = center;
		this->W_JPSI = sigma;

		*this->subtractionFactor = signalRegion/abs(sidebandRegion - signalRegion);
	}

	TBox *createTBox(double Ymax, int index = 0)
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

		TBox *region = new TBox(x1, 0., x2, Ymax);

		return region;
	}

	TCanvas *createCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSave = false)
	{
		string canvasName 	= "InvariantMass_" + string(*PassingOrFailing);
		string canvasTitle	= string(*PassingOrFailing) + " Invariant Mass";
		string saveAs 		= "../result/InvariantMass" + string(*PassingOrFailing) + ".png";

		if (drawRegions)
		{
			canvasName 	= "InvariantMass_" + string(*PassingOrFailing) + "_region";
			canvasTitle	= string(*PassingOrFailing) + " Invariant Mass with Regions";
			saveAs 		= "../result/InvariantMass" + string(*PassingOrFailing) + "_region" + ".png";
		}

		//Create canvas
		gStyle->SetCanvasPreferGL(kTRUE);
		TCanvas *c1 = new TCanvas(canvasName.data(), canvasTitle.data(), 600, 600);

		//Set margin for canvas
		c1->cd(1)->SetMargin(0.15, 0.05, 0.11, 0.07);

		//Set title size for axis and other stuffs for histogram style
		hMass->GetYaxis()->SetTitleSize(0.04);
		hMass->GetXaxis()->SetTitleSize(0.05);
		hMass->GetXaxis()->CenterTitle(true);
		//hMassAll->SetBit(TH1::kNoTitle);		//Not show histogram title
		hMass->SetMarkerStyle(20);				//Set markers style
		hMass->SetMarkerColor(kBlack);			//Set markers colors
		hMass->SetLineColor(kBlack);			//Set lines colors (for errorbars)

		//Draws histogram & fit function
		hMass->Draw("ep");

		if (fitFunction != NULL)
		{
			fb->Draw("same");
			fsG->Draw("same");
			fsCB->Draw("same");
			f->Draw("same");
		}

		//Draws information
		TLatex *tx = new TLatex();
		tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);


		if (fitFunction != NULL)
		{
			//Show chi-squared test
			tx->DrawLatex(0.61,0.60,Form("#chi^{2}/ndf = %.3g",fitResult->Chi2()/fitResult->Ndf()));
			tx->DrawLatex(0.19,0.88,Form("#bf{CMS Open Data}"));
			/*
			//Show number of particles
			tx->DrawLatex(0.61,0.48,Form("%.0f #pm %.0f J/#psi", S, dS));
			tx->DrawLatex(0.64,0.43,"(candidates)");
			*/
		}

		//Add legend
		TLegend *l = new TLegend(0.65,0.68,0.92,0.90);
		l->SetTextSize(0.04);
		l->AddEntry(hMass,	"Data"	,"lp");
		if (fitFunction != NULL)
		{
			l->AddEntry(f,		"Total Fit",	"l");
			l->AddEntry(fsG,	"Gaussian",		"l");
			l->AddEntry(fsCB,	"Crystal Ball",	"l");
			l->AddEntry(fb,		"Exp + Exp",	"l");
		}
		l->Draw();

		//Not show frame with mean, std dev
		gStyle->SetOptStat(0);

		//Get Y range of draw
		gPad->Update();
		Double_t Ymax = gPad->GetFrame()->GetY2();
		
		//Draw regions
		if (drawRegions == true)
		{
			TBox *side1 = this->createTBox(Ymax, -1);
			side1->SetFillColorAlpha(kRed, 0.35);
			side1->Draw();
			TBox *signal = this->createTBox(Ymax, 0);
			signal->SetFillColorAlpha(kGreen, 0.35);
			signal->Draw();
			TBox *side2 = this->createTBox(Ymax, 1);
			side2->SetFillColorAlpha(kRed, 0.35);
			side2->Draw();
		}

		//Writes in file
		if (shouldWrite == true)
		{
			//Here code crashes when it runs for the 2nd time
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

	void debugCout()
	{
		//Show chi-squared test (on prompt)
		cout << endl;
		cout << "Fitting overview for " << *PassingOrFailing << " " << *particleName << endl;
		cout << "- Chi2/ndf          = " << f->GetChisquare()/f->GetNDF() << endl;
		cout << "- Gaus(Sg) Position = " << resultParameters[1] << " +- " << resultParameters[2] << endl;
		cout << "- CB  (Sg) Position = " << resultParameters[5] << " +- " << resultParameters[6] << endl;
		cout << "- Position result   = " << M_JPSI << " +- " << W_JPSI << " [" << M_JPSI - W_JPSI*3 << ", " << M_JPSI + W_JPSI*3  << "]" << endl;
		cout << "- #HistIntegral     = " << hMass->Integral(0, hMass->GetNbinsX()) << endl;
		cout << "- #Total      (fit) = " << f ->Integral(xMin, xMax) * scale << endl;
		cout << "- #Background (fit) = " << fb->Integral(xMin, xMax) * scale << endl;
		cout << "- #Signal     (fit) = " << fs->Integral(xMin, xMax) * scale << endl;
		cout << endl;
	}

	InvariantMass(int *method, double *subtractionFactor, const char **particleName, const char **PassingOrFailing)
		: method(method), subtractionFactor(subtractionFactor), particleName(particleName), PassingOrFailing(PassingOrFailing)
	{}
};