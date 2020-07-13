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

#include <iostream>

using namespace std;

#include "FitFunctions.h"

//GLoal Chi2 operator for simutaneous fit
struct GlobalChi2 {
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;


	int iparPass[12] = {0,
						1,
						2,
						3,
						4,
						5,
						6,
						7,
						8,
						9,
						10,
						11,
					};

	int iparPassFail[24] = {0,
						1,
						2,
						3,
						4,
						5,
						6,
						7,	
						8,
						9,
						10,
						11,
						12,
						13,
						14,
						15,
						16,
						17,
						18,
						19,
						20,
						21,
						22,
						23,
					};

   double operator() (const double *par) const {
      double p1[12];
      for (int i = 0; i < 12; ++i) p1[i] = par[iparPass[i]];

      double p2[24];
      for (int i = 0; i < 24; ++i) p2[i] = par[iparPassFail[i]];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }

   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}
};

//Store invariant mass class
class InvariantMass{
private:
	int *method;
	double *subtractionFactor;
	const char **particleName;

	//Pointer to objects
	PassingFailing* ObjPass = NULL;
	PassingFailing* ObjFail = NULL;
	PassingFailing* ObjAll = NULL;

	//Self fit functions
	TF1* fitFunctionPass 	= NULL;
	TF1* fitFunctionFail 	= NULL;
	TF1* fitFunctionAll		= NULL;

	ROOT::Fit::FitResult fitResult;

public:

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

	int			nBins;
	int			decimals = 4;
	double 		xMin;
	double		xMax;

	void defineNumbers(int nBins, double xMin, double xMax, int decimals = 4)
	{
		this->nBins 	= nBins;
		this->xMin 		= xMin;
		this->xMax 		= xMax;
		this->decimals 	= decimals;
	}

	void createMassHistogram(TH1D* &hMass, const char* PassingOrFailing)
	{
		string hName 			= string(PassingOrFailing) + string(*particleName) + "InvariantMass";
		string hTitle 			= "Invariant Mass (" + string(PassingOrFailing) + ")";
		string yAxisTitleForm 	= "Events / (%1." + to_string(decimals) + "f GeV/c^{2})";

		//Create histogram
		hMass = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		hMass->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), hMass->GetBinWidth(0)));
		hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
	}

	void createAllMassHistograms()
	{
		createMassHistogram((*this->ObjPass).hMass, "Passing tracker");
		createMassHistogram((*this->ObjFail).hMass, "Failing tracker");
		createMassHistogram((*this->ObjAll).hMass,  "All");
	}

	ROOT::Fit::FitResult doFit()
	{
		TH1D* &hPass = (*this->ObjPass).hMass;
		TH1D* &hFail = (*this->ObjFail).hMass;
		TH1D* &hPassFail = (*this->ObjAll).hMass;

		//Passing Fitting
		TF1* &fPass = this->fitFunctionPass;
		fPass = new TF1("FitFunction_Pass", FitFunctions::Merged::Pass_InvariantMass, xMin, xMax, 12);
		
		//Rename parameters
		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
		for (int i = 0; i < arraySize; i++)
		{
			fPass->SetParName(i, this->fittingParName[i]);
		}

		//Failing Fitting
		TF1* &fFail = this->fitFunctionFail;
		fFail = new TF1("FitFunction_Fail", FitFunctions::Merged::Fail_InvariantMass, xMin, xMax, 12);
		
		//Rename parameters
		for (int i = 0; i < arraySize; i++)
		{
			fFail->SetParName(i, this->fittingParName[i]);
		}

		//Both Fitting
		TF1* &fPassFail = this->fitFunctionAll;
		fPassFail = new TF1("FitFunction_Both", FitFunctions::Merged::Both_InvariantMass, xMin, xMax, 24);
		
		//Rename parameters
		for (int i = 0; i < arraySize*2; i++)
		{
			fPassFail->SetParName(i, this->fittingParName[i%arraySize]);
		}

		//Simultaneuos fit
		ROOT::Math::WrappedMultiTF1 wfPass(*fPass, 1);
		ROOT::Math::WrappedMultiTF1 wfPassFail(*fPassFail, 1);

		ROOT::Fit::DataOptions opt;

		ROOT::Fit::DataRange rangePass;
		rangePass.SetRange(xMin, xMax);
		ROOT::Fit::BinData dataPass(opt, rangePass);
		ROOT::Fit::FillData(dataPass, hPass);

		ROOT::Fit::DataRange rangePassFail;
		rangePassFail.SetRange(xMin, xMax);
		ROOT::Fit::BinData dataPassFail(opt, rangePassFail);
		ROOT::Fit::FillData(dataPassFail, hPassFail);

		ROOT::Fit::Chi2Function chi2_Pass(dataPass, wfPass);
		ROOT::Fit::Chi2Function chi2_PassFail(dataPassFail, wfPassFail);

		GlobalChi2 globalChi2(chi2_Pass, chi2_PassFail);

		ROOT::Fit::Fitter fitter;

		double par0[24] = {340.2,
							3.09,
							0.037,
							1.824,
							1.034,
							3.093,
							0.022,
							8322.27,
							-0.217,
							1.915,
							263.185,
							0.061,
							340.2,
							3.09,
							0.037,
							1.824,
							1.034,
							3.093,
							0.022,
							8322.27,
							-0.217,
							1.915,
							263.185,
							0.061
						};

		//Create before the parameter settings in order to fix or set range on them
		fitter.Config().SetParamsSettings(24, par0);

		//Fit FCN function directly
		//(specify optionally data size and flag to indicate that is a chi2 fit)
		fitter.FitFCN(24, globalChi2, 0, dataPass.Size() + dataPassFail.Size(), true);
		ROOT::Fit::FitResult result = fitter.Result();
		result.Print(std::cout);

		//Update Failing fit
		double resultParameters[24];
		fPassFail->GetParameters(resultParameters);
		fFail->SetParameters(&resultParameters[12]);	//Get only failing part

		//Update result
		this->fitResult = result;

		return this->fitResult;
	}

	void updateMassFor(PassingFailing* &Obj, TF1* &signalFit)
	{
		double value;
		double fwhm;

		if (*this->method == 1)
		{
			//Get value and uncertain of signal by histogram
			int bin0 = Obj->hMass->GetMaximumBin();
			value    = Obj->hMass->GetBinCenter(bin0);
			int bin1 = Obj->hMass->FindFirstBinAbove(Obj->hMass->GetMaximum()/2);
			int bin2 = Obj->hMass->FindLastBinAbove(Obj->hMass->GetMaximum()/2);
			fwhm     = Obj->hMass->GetBinCenter(bin2) - Obj->hMass->GetBinCenter(bin1);
		}

		if (*this->method == 2)
		{
			//Get value and uncertain of signal by fitting
			value     = signalFit->GetMaximumX();
			double x1 = signalFit->GetX(signalFit->GetMaximum()/2);
			double x2 = signalFit->GetX(signalFit->GetMaximum()/2, x1+0.0001, value + x1*3);
			fwhm      = x2 - x1;
		}

		double sigma = fwhm/2.355;

		//For regions of sidebande subtraction
		Obj->M_JPSI = value;
		Obj->W_JPSI = sigma;

		//For sideband subtraction
		Obj->subtractionFactor = (*this->ObjPass).signalRegion/abs((*this->ObjPass).sidebandRegion - (*this->ObjPass).signalRegion);
	}

	void updateMassAll()
	{
		updateMassFor(ObjPass, fitFunctionPass);
		updateMassFor(ObjFail, fitFunctionFail);
		updateMassFor(ObjAll,  fitFunctionAll);
	}

	void drawCanvasQuarter(TCanvas* &canvas, bool drawRegions, int quarter, TH1D* &histo, PassingFailing* &Obj, TF1* &fit, int color = kBlue)
	{
		float margins[4] = {0.10, 0.02, 0.09, 0.07};

		canvas->cd(quarter);
		canvas->cd(quarter)->SetMargin(margins[0], margins[1], margins[2], margins[3]);

		histo->SetMarkerStyle(20);		//Set markers style
		histo->SetMarkerColor(kBlack);	//Set markers colors
		histo->SetLineColor(kBlack);	//Set errobars color
		histo->Draw("ep");

		//Draws information
		TLatex *tx = new TLatex();
		tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);
		tx->DrawLatex(0.14,0.88,Form("#bf{CMS Open Data}"));

		//Add legend
		TLegend *l = new TLegend(0.75,0.80,0.92,0.90);
		l->SetTextSize(0.04);
		l->AddEntry(histo, "Data", "lp");
		if (fit != NULL)
		{
			l->AddEntry(fit, "Total Fit", "l");
		}
		l->Draw();

		//Draw fit
		if (fit != NULL)
		{
			fit->SetNpx(1000);
			fit->SetLineColor(color);
			fit->SetLineStyle(kSolid);
			fit->Draw("same");
		}

		//Draw regions
		if (drawRegions == true)
		{
			//Get Y range of draw
			gPad->Update();
			Double_t Ymax = gPad->GetFrame()->GetY2();

			//Draw regions
			TBox *side1  = Obj->createTBox(Ymax, -1);
			side1->SetFillColorAlpha(kRed, 0.35);
			side1->Draw();
			TBox *signal = Obj->createTBox(Ymax, 0);
			signal->SetFillColorAlpha(kGreen, 0.35);
			signal->Draw();
			TBox *side2  = Obj->createTBox(Ymax, 1);
			side2->SetFillColorAlpha(kRed, 0.35);
			side2->Draw();
		}
	}

	TCanvas *createCanvas(bool drawRegions = false, bool shouldWrite = false, bool shouldSave = false)
	{
		string canvasName 	= "InvariantMass_" + string("Tracker");
		string canvasTitle	= "Invariant Mass " + string("Tracker");
		string saveAs 		= "../result/InvariantMass_" + string("Tracker") + ".png";

		if (drawRegions)
		{
			canvasName 	= "InvariantMass_" + string("Tracker") + "_region";
			canvasTitle	= "Invariant Mass " + string("Tracker") + " with Regions";
			saveAs 		= "../result/InvariantMass_" + string("Tracker") + "_region" + ".png";
		}

		//Create canvas
		gStyle->SetCanvasPreferGL(kTRUE);
		TCanvas *c1 = new TCanvas(canvasName.data(), canvasTitle.data(), 1200, 800);
		c1->Divide(2,2);

		this->drawCanvasQuarter(c1, drawRegions, 1, (*this->ObjPass).hMass, ObjPass, fitFunctionPass, kGreen);
		this->drawCanvasQuarter(c1, drawRegions, 2, (*this->ObjFail).hMass, ObjFail, fitFunctionFail, kRed);
		this->drawCanvasQuarter(c1, drawRegions, 3, (*this->ObjAll).hMass,  ObjAll,  fitFunctionAll), kBlue;

		c1->cd(4);

		//Draws information
		TLatex *tx = new TLatex();
		//tx->SetTextSize(0.04);
		tx->SetTextAlign(12);
		tx->SetTextFont(42);
		tx->SetNDC(kTRUE);
		tx->DrawLatex(0.61,0.60,Form("#chi^{2}/ndf = %.3g",(this->fitResult).Chi2()/(this->fitResult).Ndf()));

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

	InvariantMass(int *method, const char **particleName, PassingFailing *Passing, PassingFailing *Failing, PassingFailing *All)
		: method(method), particleName(particleName), ObjPass(Passing), ObjFail(Failing), ObjAll(All)
		{}
};