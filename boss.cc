 /*
!--------------------------------
!proposito: ajustar funções sobre histogramas no auxilio do tag and probe
!--------------------------------	
!autor: Allan Jales
!--------------------------------
*/

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"

#include <TF1NormSum.h>

#include <iostream>

using namespace std;



//-------------------------------------
// Fit functions for Invariant Mass
//-------------------------------------

class FitFunctions
{
public:
	class Primary
	{
	public:
		//Gaussian function
		static Double_t Gaus(Double_t *x, Double_t *par)
		{
			//par[0] = height
			//par[1] = position
			//par[2] = sigma
			Double_t gaus = par[0]*TMath::Gaus(x[0],par[1],par[2],1);
			return gaus;
		}

		//Polynomial function
		static Double_t Pol1(Double_t *x, Double_t *par)
		{
			//par[0] = b
			//par[1] = a
			Double_t pol = par[0] + par[1]*x[0];
			return pol;
		}

		//Exponential function
		static Double_t Exp(Double_t *x, Double_t *par)
		{
			//par[0] = height
			//par[1] = width
			Double_t exp = par[0] * TMath::Exp(par[1]*x[0]);
			return exp;
		}

		//crystall ball function
		static Double_t CrystalBall(Double_t *x,Double_t *par)
		{
			//par[0] = alpha
			//par[1] = n
			//par[2] = mean
			//par[3] = sigma
			//par[4] = Yield
			Double_t t = (x[0]-par[2])/par[3];
			if (par[0] < 0) t = -t;
			Double_t absAlpha = fabs((Double_t)par[0]);
			if (t >= -absAlpha)
			{
				return par[4]*exp(-0.5*t*t);
			}
			else
			{
				Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
				Double_t b = par[1]/absAlpha - absAlpha;
				return par[4]*(a/TMath::Power(b - t, par[1]));
			}
		}
	};
	class Merged
	{
	public:
		//Fit function for signal for Invariant Mass Probe
		static Double_t Signal_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Primary::Gaus(x,par) + FitFunctions::Primary::CrystalBall(x, &par[3]);
		}

		//Fit function for background for Invariant Mass Probe
		static Double_t Background_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Primary::Exp(x,par) + FitFunctions::Primary::Exp(x, &par[2]);
		}

		//Fit function for signal & background for Invariant Mass Probe
		static Double_t FFit_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Merged::Signal_InvariantMassAll(x,par) + FitFunctions::Merged::Background_InvariantMassAll(x, &par[8]);
		}
	};
};

//-------------------------------------
// Main functions
//-------------------------------------

class ParticleSelector{
public:
	const char *PassingOrFailing = "Unknown";

	//Jpsi constants (PDG values)
	//const double M_JPSI = 3.097;
	//const double W_JPSI = 0.010;

	//Jpsi constants (PDG values)
	double M_JPSI = 3.097;
	double W_JPSI = 0.010;

	double signalRegionEnd 		= 10.0;
	double sidebandRegionEnd	= 20.0;

	int 	count_signalRegion = 0;
	int 	count_sidebandRegion = 0;

	double 	Signal  = 0;
	double 	dSignal = 0;

	void calculus()
	{
		Signal  = count_signalRegion - count_sidebandRegion;
		dSignal = sqrt(count_signalRegion + count_sidebandRegion);
	}

	bool isInSignalRegion(Double_t InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) < W_JPSI * signalRegionEnd)
			return true;

		return false;
	}

	bool isInSidebandRegion(Double_t InvariantMass)
	{
		if (fabs(InvariantMass - M_JPSI) > W_JPSI * signalRegionEnd && fabs(InvariantMass - M_JPSI) < W_JPSI * sidebandRegionEnd)
			return true;

		return false;
	}

	void debugCout()
	{
		cout << endl;
		cout << "Candidates by sideband subtraction in " << PassingOrFailing 			<< endl;
		cout << "#Total           = "	<< count_signalRegion + count_sidebandRegion	<< endl;
		cout << "#Signal   region = "	<< count_signalRegion							<< endl;
		cout << "#Sideband region = "	<< count_sidebandRegion							<< endl;
		printf( "#Signal          = %.0f +- %.0f\n", Signal, dSignal);
		cout << endl;
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
				dx1 = -sidebandRegionEnd;
				dx2 = -signalRegionEnd;
				break;
			case 0:
				dx1 = -signalRegionEnd;
				dx2 = +signalRegionEnd;
				break;
			case 1:
				dx1 = +sidebandRegionEnd;
				dx2 = +signalRegionEnd;
				break;
		}

		double x1 = M_JPSI + W_JPSI * dx1;
		double x2 = M_JPSI + W_JPSI * dx2;

		TBox *region = new TBox(x1, 0., x2, Ymax);

		return region;
	}
};

class InvariantMassClass{
private:
	TF1* fitFunction 		= NULL;
	TF1* fitFunctionSig 	= NULL;
	TF1* fitFunctionBack 	= NULL;
	TFitResultPtr fitResult = 0;
	
	//References
	TF1* &f  = fitFunction;
	TF1* &fs = fitFunctionSig;
	TF1* &fb = fitFunctionBack;

	//Bins per each x axis unit 
	double scale = 0;

public:
	const char *PassingOrFailing = NULL;
	const char *particle = "Muon";

	int			nBins;
	int			decimals = 4;
	double 		xMin;
	double		xMax;

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

			"Exp1(Bg) Lambda  ",
			"Exp1(Bg) a       ",
			"Exp2(Bg) Lambda  ",
			"Exp2(Bg) a       "
		};

	Double_t resultParameters[12];

	ParticleSelector Sideband;

	void define(const char *PassingOrFailing)
	{
		this->PassingOrFailing 			= PassingOrFailing;
		this->Sideband.PassingOrFailing = PassingOrFailing;
	}

	void defineNumbers(Int_t nBins, Double_t xMin, Double_t xMax, int decimals = 4)
	{
		this->nBins 	= nBins;
		this->xMin 		= xMin;
		this->xMax 		= xMax;
		this->decimals 	= decimals;
	}

	void createMassHistogram()
	{
		string hName 			= string(PassingOrFailing) + string(particle) + "InvariantMass";
		string hTitle 			= "Invariant Mass (" + string(PassingOrFailing) + ")";
		string yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f GeV/c^{2})";

		//Create histogram
		hMass = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		hMass->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), hMass->GetBinWidth(0)));
		hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");

		//Get scale from histogram
		scale = 1/hMass->GetBinWidth(0);
	}

	void fill(double value)
	{
		this->hMass->Fill(value);
	}

	void fit()
	{
		//Create fit Function
		fitFunction = new TF1("FitFunction", FitFunctions::Merged::FFit_InvariantMassAll, xMin, xMax, 12);

		//Get size of array
		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);

		//Rename parameters and set value
		for (int i = 0; i < arraySize; i++)
		{
			f->SetParName(i, fittingParName[i]);
		}

		//Resolution of fit function
		f->SetNpx(1000);
		
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

		//Fit Color
		f->SetLineColor(kBlue);

		//Fit the function
		fitResult = hMass->Fit(fitFunction, "RNS", "", xMin, xMax);

		//Get parameters from fit function and put it in par variable
		fitFunction->GetParameters(resultParameters);
		
		//Signal Fitting
		fs = new TF1("FitFunction_Signal", FitFunctions::Merged::Signal_InvariantMassAll, xMin, xMax, 8);
		fs->SetNpx(1000);							//Resolution of background fit function
		fs->SetParameters(resultParameters);		//Get only background part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Merged::Background_InvariantMassAll, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[8]);	//Get only background part
		fb->SetLineColor(kBlue); 					//Fit Color
		fb->SetLineStyle(kDashed);					//Fit Style

		double position = f->GetMaximumX();
		double x1 = f->GetX(f->GetMaximum()/2);
		double x2 = f->GetX(f->GetMaximum()/2, x1+0.0001, position + x1*3);
		double fwhm = x2 - x1;
		double sigma = fwhm/2;

		this->Sideband.M_JPSI = position;
		this->Sideband.W_JPSI = sigma;
		this->Sideband.signalRegionEnd = 3.;
		this->Sideband.sidebandRegionEnd = 6.;
	}

	TCanvas *createCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSave = false)
	{
		string canvasName 	= "InvariantMass_" + string(PassingOrFailing);
		string canvasTitle	= string(PassingOrFailing) + " Invariant Mass";
		string saveAs 		= "../result/InvariantMass" + string(PassingOrFailing) + ".png";

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
		//fs->Draw("same");
		fb->Draw("same");
		f->Draw("same");

		//Draws information
		TLatex *tx = new TLatex();
		tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);

		//Show chi-squared test
		tx->DrawLatex(0.61,0.60,Form("#chi^{2}/ndf = %.3g",fitResult->Chi2()/fitResult->Ndf()));
		/*
		//Show number of particles
		tx->DrawLatex(0.61,0.48,Form("%.0f #pm %.0f J/#psi", S, dS));
		tx->DrawLatex(0.64,0.43,"(candidates)");
		*/

		//Add legend
		TLegend *l = new TLegend(0.65,0.77,0.92,0.90);
		l->SetTextSize(0.04);
		l->AddEntry(hMass,	"Data"	,"lp");
		l->AddEntry(f,		"Fitting"	,"l");
		//l->AddEntry(fs,	"Signal"	,"l");
		l->AddEntry(fb,		"Background","l");
		l->Draw();

		//Not show frame with mean, std dev
		gStyle->SetOptStat(0);

		//Get Y range of draw
		gPad->Update();
		Double_t Ymax = gPad->GetFrame()->GetY2();
		
		//Draw regions
		if (drawRegions)
		{
			/*
			if (true)
			{
				if (true)
				{
					double position = f->GetMaximumX();
					double x1 = f->GetX(f->GetMaximum()/2);
					double x2 = f->GetX(f->GetMaximum()/2, x1, position + x1*3);
					cout << "Limits: " << x1 << " ~ " << x2 << endl;
					double fwhm = x2 - x1;
					double sigma = fwhm/2;

					TBox *signal = new TBox(position - sigma * 3, 0., position + sigma * 3, Ymax);
					signal->SetFillColorAlpha(kGreen, 0.2);
					signal->Draw();
					TBox *signal1 = new TBox(position - sigma * 2, 0., position + sigma * 2, Ymax);
					signal1->SetFillColorAlpha(kGreen, 0.2);
					signal1->Draw();
					TBox *signal2 = new TBox(position - sigma * 1, 0., position + sigma * 1, Ymax);
					signal2->SetFillColorAlpha(kGreen, 0.2);
					signal2->Draw();
				}
				else
				{
					double position = f->GetMaximumX();
					int bin1 = hMass->FindFirstBinAbove(hMass->GetMaximum()/2);
					int bin2 = hMass->FindLastBinAbove(hMass->GetMaximum()/2);
					double fwhm = hMass->GetBinCenter(bin2) - hMass->GetBinCenter(bin1);
					double sigma = fwhm/2;

					TBox *signal = new TBox(position - sigma * 3, 0., position + sigma * 3, Ymax);
					signal->SetFillColorAlpha(kGreen, 0.2);
					signal->Draw();
					TBox *signal1 = new TBox(position - sigma * 2, 0., position + sigma * 2, Ymax);
					signal1->SetFillColorAlpha(kGreen, 0.2);
					signal1->Draw();
					TBox *signal2 = new TBox(position - sigma * 1, 0., position + sigma * 1, Ymax);
					signal2->SetFillColorAlpha(kGreen, 0.2);
					signal2->Draw();
				}
			}
			else
			*/
			{
				TBox *side1 = Sideband.createTBox(Ymax, -1);
				side1->SetFillColorAlpha(kRed, 0.35);
				side1->Draw();
				TBox *signal = Sideband.createTBox(Ymax, 0);
				signal->SetFillColorAlpha(kGreen, 0.35);
				signal->Draw();
				TBox *side2 = Sideband.createTBox(Ymax, 1);
				side2->SetFillColorAlpha(kRed, 0.35);
				side2->Draw();
			}
		}

		//debugCout();

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

	bool isInSignalRegion(Double_t InvariantMass)
	{
		return this->Sideband.isInSignalRegion(InvariantMass);
	}

	bool isInSidebandRegion(Double_t InvariantMass)
	{
		return this->Sideband.isInSidebandRegion(InvariantMass);
	}

	void debugCout()
	{
		double position = f->GetMaximumX();
		int bin1 = hMass->FindFirstBinAbove(hMass->GetMaximum()/2);
		int bin2 = hMass->FindLastBinAbove(hMass->GetMaximum()/2);
		double fwhm = hMass->GetBinCenter(bin2) - hMass->GetBinCenter(bin1);
		double sigma = fwhm/2;

		//Show chi-squared test (on prompt)
		cout << endl;
		cout << "Fitting overview for " << PassingOrFailing << endl;
		cout << "Chi2/ndf          = " << f->GetChisquare()/f->GetNDF() << endl;
		cout << "Gaus(Sg) Position = " << resultParameters[1] << endl;
		cout << "Gaus(Sg) Sigma    = " << resultParameters[2] << endl;
		cout << "CB  (Sg) Mean     = " << resultParameters[5] << endl;
		cout << "CB  (Sg) Sigma    = " << resultParameters[6] << endl;
		cout << "Position result   = " << position << "+-" << sigma << endl;
		cout << "Limits            = " << position - sigma*3 << " ~ " << position + sigma*3 << endl;
		cout << endl;

		//Show integrals
		cout << "Candidates by integration for " << PassingOrFailing << endl;
		cout << "#HistIntegral = " << hMass->Integral(0, hMass->GetNbinsX()) << endl;
		cout << "#Total        = " << f ->Integral(xMin, xMax) * scale << endl;
		cout << "#Background   = " << fb->Integral(xMin, xMax) * scale << endl;
		cout << "#Signal       = " << fs->Integral(xMin, xMax) * scale << endl;
		cout << endl;
	}
};

class Histograms{
public:
	const char *quantityName 			= NULL;
	const char *extendedQuantityName	= NULL;
	const char *quantityUnit			= NULL;
	const char *xAxisName				= NULL;
	const char *tagOrProbe				= NULL;
	const char *PassingOrFailing		= NULL;
	const char *particle 				= "Muon";

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

	void defineTexts(const char *quantityName, const char *xAxisName, const char *quantityUnit,const char *extendedQuantityName, const char *tagOrProbe, const char *PassingOrFailing, const char *particle = "Muon")
	{
		this->quantityName 			= quantityName;
		this->extendedQuantityName 	= extendedQuantityName;
		this->xAxisName 			= xAxisName;
		this->quantityUnit			= quantityUnit;
		this->tagOrProbe 			= tagOrProbe;
		this->PassingOrFailing 		= PassingOrFailing;
		this->particle 				= particle;
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
		//Define parameters
		string hName 		= PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + "SigBack";
		string hTitle 		= string(extendedQuantityName) + " (" + string(PassingOrFailing) + " " + string(tagOrProbe) + ")";
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
		hSigBack = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		hSigBack->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), hSigBack->GetBinWidth(0)));
		hSigBack->GetXaxis()->SetTitle(xAxisTitle.data());
	}

	void createBackHistogram()
	{
		//Define parameters
		string hName = PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + "Back";
	
		//Create histogram
		hBack = (TH1D*) hSigBack->Clone("hBack");
		hBack->SetName(hName.data());
	}

	void subtractSigHistogram()
	{
		double factor = signalRegionEnd/abs(sidebandRegionEnd - signalRegionEnd);

		//Define parameters
		string hName = PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + "Sig";
	
		//Create histogram
		hSig = (TH1D*) hSigBack->Clone("hSig");
		hSig->SetName(hName.data());
		hSig->Add(hBack,-factor);
	}

	TCanvas *createDividedCanvas(bool shouldWrite = false, bool shouldSave = true)
	{
		string canvasName 	= PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName);
		string titleLeft 	= string(extendedQuantityName) + " (" + string(PassingOrFailing) + " " + string(tagOrProbe) + ")";
		string titleRight 	= string(extendedQuantityName) + " of Signal (" + string(PassingOrFailing) + " " + string(tagOrProbe) + ")";
		string saveAs 		= "../result/" + string(quantityName) + string(PassingOrFailing) + string(tagOrProbe) + ".png";

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
		}
		else
		{
			tx1_1->SetTextAlign(22);	//Align center, center
			tx1_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
			tx1_1->DrawLatex(0.55,0.45,Form("%g entries (signal)",		hSig->GetEntries()));
			tx1_1->DrawLatex(0.55,0.40,Form("%g entries (background)",	hBack->GetEntries()));
		}
		

		//Select canvas part and set margin
		c1->cd(2);
		c1->cd(2)->SetMargin(0.14, 0.03, 0.11, 0.07);

		//Same range as comparision and draws
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
		}
		else
		{
			tx1_2->SetTextAlign(22);	//Align center, center
			tx1_2->DrawLatex(0.55,0.5,Form("%g entries (signal)", hSig->GetEntries()));
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

		string pName 	= PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + "Efficiency";
		string pTitle 	= string(extendedQuantityName) + " Efficiency (" + string(PassingOrFailing) + " " + string(tagOrProbe) + ")";

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
		string canvasName 	= PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + "Efficiency";
		string canvasTitle 	= string(extendedQuantityName) + " Efficiency (" + string(PassingOrFailing) + " " + string(tagOrProbe) + ")";
		string saveAs 		= "../result/" + string(quantityName) + string(PassingOrFailing) + string(tagOrProbe) + "_Efficiency.png";

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

	void debugCout()
	{
		const int minLegendSpace = 21;

		//Set information what are shown
		string legend = "#" + string(PassingOrFailing) + " " + string(tagOrProbe) + " " + string(quantityName);
		if(strlen(legend.data()) < minLegendSpace-3)
		{
			legend.append(minLegendSpace - 3 - strlen(legend.data()), ' ');
		}
		legend += " = ";

		//Show information
		cout << legend << hSigBack->GetEntries() - hSig->GetEntries() - hBack->GetEntries() << endl;
	}
};

//Bunch of 3 histogram class for Tag and Probe quantities
class TagProbe{
public:
	const char *tagOrProbe			= NULL;
	const char *PassingOrFailing	= NULL;

	Histograms Pt;
	Histograms Eta;
	Histograms Phi;

	void define(const char *tagOrProbe, const char *PassingOrFailing)
	{
		this->tagOrProbe 		= tagOrProbe;
		this->PassingOrFailing 	= PassingOrFailing;
	}

	void defineHistogramsTexts()
	{
		this->Pt .defineTexts("Pt",  "P_{t}",	"GeV/c", 	"Transversal Momentum", tagOrProbe, PassingOrFailing);
		this->Eta.defineTexts("Eta", "#eta", 	"", 		"Pseudorapidity", 		tagOrProbe, PassingOrFailing);
		this->Phi.defineTexts("Phi", "#phi", 	"", 		"Azimuthal Angle", 		tagOrProbe, PassingOrFailing);
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

	void subtractSigHistogram()
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

	void setRange(double signalRegionEnd, double sidebandRegionEnd)
	{
		this->Pt .signalRegionEnd 	= signalRegionEnd;
		this->Pt .sidebandRegionEnd = sidebandRegionEnd;
		this->Eta.signalRegionEnd 	= signalRegionEnd;
		this->Eta.sidebandRegionEnd = sidebandRegionEnd;
		this->Phi.signalRegionEnd 	= signalRegionEnd;
		this->Phi.sidebandRegionEnd = sidebandRegionEnd;
	}

	TagProbe(const char *tagOrProbe)
	{
		this->tagOrProbe 		= tagOrProbe;
	}
};

//Holder for tag probe class
class Particle{
public:
	const char *PassingOrFailing = NULL;

	TagProbe Tag{"Tag"};
	TagProbe Probe{"Probe"};
	InvariantMassClass Mass;

	void define(const char *PassingOrFailing)
	{
		this->PassingOrFailing 	= PassingOrFailing;
		this->Tag  .PassingOrFailing 	= PassingOrFailing;
		this->Probe.PassingOrFailing 	= PassingOrFailing;
		this->Mass.define(PassingOrFailing);
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

	void subtractSigHistogram()
	{
		this->Tag  .subtractSigHistogram();
		this->Probe.subtractSigHistogram();
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
		/*
		double &InvariantMass 	= quantities[0];
		double &TagMuon_Pt		= quantities[1];
		double &TagMuon_Eta		= quantities[2];
		double &TagMuon_Phi		= quantities[3];
		double &ProbeMuon_Pt	= quantities[4];
		double &ProbeMuon_Eta	= quantities[5];
		double &ProbeMuon_Phi	= quantities[6];
		*/

		//If is inside signal region
		if (this->Mass.isInSignalRegion(InvariantMass))
		{
			this->fillSigBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			this->Mass.Sideband.count_signalRegion++;
		}

		//If is inside sideband region
		if (this->Mass.isInSidebandRegion(InvariantMass))
		{
			this->fillBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			this->Mass.Sideband.count_sidebandRegion++;
		}
	}

	void setRange()
	{
		this->Tag  .setRange(Mass.Sideband.signalRegionEnd, Mass.Sideband.sidebandRegionEnd);
		this->Probe.setRange(Mass.Sideband.signalRegionEnd, Mass.Sideband.sidebandRegionEnd);
	}

	void debugCout()
	{
		this->Tag  .debugCout();
		this->Probe.debugCout();
	}

	Particle(const char *PassingOrFailing)
	{
		this->define(PassingOrFailing);
	}
};

//Select particles, draws and save histograms
void generateHistograms(bool shouldDrawInvariantMassCanvas = true, bool shouldDrawQuantitiesCanvas = true, bool shouldDrawEfficiencyCanvas = true)
{
	TFile *file0 = TFile::Open("../data_histoall.root");		//Opens the file
	TTree *TreePC = (TTree*)file0->Get("demo/PlotControl");		//Opens TTree of file
	TTree *TreeAT = (TTree*)file0->Get("demo/AnalysisTree");	//Opens TTree of file
	
	//Create variables for PlotControl
	double ProbeMuon_Pt;
	double ProbeMuon_Eta;
	double ProbeMuon_Phi;
	double TagMuon_Pt;
	double TagMuon_Eta;
	double TagMuon_Phi;
	double InvariantMass;

	//Create variables for AnalysisTree
	int PassingProbeTrackingMuon;
	int PassingProbeStandAloneMuon;
	int PassingProbeGlobalMuon;

	//Assign variables
	TreePC->SetBranchAddress("ProbeMuon_Pt",				&ProbeMuon_Pt);
	TreePC->SetBranchAddress("ProbeMuon_Eta",				&ProbeMuon_Eta);
	TreePC->SetBranchAddress("ProbeMuon_Phi",				&ProbeMuon_Phi);
	TreePC->SetBranchAddress("TagMuon_Pt",					&TagMuon_Pt);
	TreePC->SetBranchAddress("TagMuon_Eta",					&TagMuon_Eta);
	TreePC->SetBranchAddress("TagMuon_Phi",					&TagMuon_Phi);
	TreePC->SetBranchAddress("InvariantMass",				&InvariantMass);
	TreeAT->SetBranchAddress("PassingProbeTrackingMuon",	&PassingProbeTrackingMuon);
	TreeAT->SetBranchAddress("PassingProbeStandAloneMuon",	&PassingProbeStandAloneMuon);
	TreeAT->SetBranchAddress("PassingProbeGlobalMuon",		&PassingProbeGlobalMuon);

	//Create a object
	Particle PassingParticles("Passing");
	PassingParticles.defineHistogramsTexts();
	PassingParticles.defineHistogramsNumbers();
	PassingParticles.createSigBackHistograms();
	PassingParticles.createBackHistograms();
	PassingParticles.Mass.defineNumbers(240, 2.8, 3.4);
	PassingParticles.Mass.createMassHistogram();

	Particle FailingParticles("Failing");
	FailingParticles.defineHistogramsTexts();
	FailingParticles.defineHistogramsNumbers();
	FailingParticles.createSigBackHistograms();
	FailingParticles.createBackHistograms();
	FailingParticles.Mass.defineNumbers(240, 2.8, 3.4);
	FailingParticles.Mass.createMassHistogram();

	//Loop between the components
	for (int i = 0; i < TreePC->GetEntries(); i++)
	{
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4 && !PassingProbeStandAloneMuon && !PassingProbeGlobalMuon)
		{
			if (PassingProbeTrackingMuon)
				PassingParticles.Mass.fill(InvariantMass);
			else
				FailingParticles.Mass.fill(InvariantMass);
		}
	}

	/*
	//Assign pointers
	double* quantities[7];
	quantities[0] = &InvariantMass;
	quantities[1] = &TagMuon_Pt;
	quantities[2] = &TagMuon_Eta;
	quantities[3] = &TagMuon_Phi;
	quantities[4] = &ProbeMuon_Pt;
	quantities[5] = &ProbeMuon_Eta;
	quantities[6] = &ProbeMuon_Phi;
	*/

	//Fit
	PassingParticles.Mass.fit();
	FailingParticles.Mass.fit();

	PassingParticles.setRange();
	FailingParticles.setRange();

	//Loop between the components again
	for (int i = 0; i < TreePC->GetEntries(); i++)
	{
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4 && !PassingProbeStandAloneMuon && !PassingProbeGlobalMuon)
		{
			if (PassingProbeTrackingMuon)
				PassingParticles.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			else
				FailingParticles.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}

	//Number of particles in signal and uncertain
	PassingParticles.Mass.Sideband.calculus();
	FailingParticles.Mass.Sideband.calculus();

	//Debug
	PassingParticles.Mass.debugCout();
	FailingParticles.Mass.debugCout();

	//Create signal histograms
	PassingParticles.subtractSigHistogram();
	FailingParticles.subtractSigHistogram();

	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Create file root to store generated files
	TFile *generatedFile = TFile::Open("../result/generated_hist.root","RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->   cd("canvas/");


	if (shouldDrawInvariantMassCanvas)
	{
		PassingParticles.Mass.createCanvas(true, true, true);
		FailingParticles.Mass.createCanvas(true, true, true);
	}

	if (shouldDrawQuantitiesCanvas)
	{
		PassingParticles.createDividedCanvas(true, true);
		FailingParticles.createDividedCanvas(true, true);
	}

	//Debug
	PassingParticles.Mass.Sideband.debugCout();
	FailingParticles.Mass.Sideband.debugCout();

	//Integrate function to get number of particles in it
	cout << endl;
	cout << "Checking histograms number inconsistency (should be 0)" << endl;
	PassingParticles.debugCout();
	FailingParticles.debugCout();


	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	//Write histograms on file
	PassingParticles.write(true, true, true);
	FailingParticles.write(true, true, true);


	//Save plots
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	PassingParticles.createEfficiencyPlot();
	FailingParticles.createEfficiencyPlot();


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");


	if (shouldDrawEfficiencyCanvas)
	{
		PassingParticles.createEfficiencyCanvas();
		FailingParticles.createEfficiencyCanvas();
	}

	//Close files
	generatedFile->Close();
}

//Test
void exibir(double ptr[])
{
	for (int i = 0; i < 4; i++)
	{
		cout << ptr[i] << endl;
	}
}

void test()
{
	cout << "TEST!" << endl;

	double var1 = 1.2;
	double var2 = 1.4;
	double foo  = 1.3;
	double bar  = 1.5;

	double ptr[4];

	ptr[0] = var1;
	ptr[1] = var2;
	ptr[2] = foo;
	ptr[3] = bar;

	var2 = 12.0;

	exibir(ptr);
}

//Call functions
void boss()
{
	generateHistograms(true, true, true);
}