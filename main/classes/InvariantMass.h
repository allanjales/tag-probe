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
	TH1D* hMass 	  = NULL;
	TF1*  fitFunction = NULL;
	TF1*  fitSignal   = NULL;

	//Value and sigma of invariant mass
	double M_JPSI = 0.;
	double W_JPSI = 0.;

	//Times sigma
	double signalRegion 	= 3.0;
	double sidebandRegion	= 6.0;


	//Regions
	double sidebandRegion1_x1  = 8.50;
	double sidebandRegion1_x2  = 9.19;
	double signalRegion_x1     = 9.19;
	double signalRegion_x2     = 9.70;
	double sidebandRegion2_x1  = 10.6;
	//double sidebandRegion2_x1  = 9.70;
	double sidebandRegion2_x2  = 11.2;

	//--- For sideband subtraction ---

	bool isInSignalRegion(double InvariantMass)
	{
		if (InvariantMass > this->signalRegion_x1 && InvariantMass < this->signalRegion_x2)
			return true;

		return false;
	}

	bool isInSidebandRegion(double InvariantMass)
	{
		if ((InvariantMass > this->sidebandRegion1_x1 && InvariantMass < this->sidebandRegion1_x2) ||
			(InvariantMass > this->sidebandRegion2_x1 && InvariantMass < this->sidebandRegion2_x2))
			return true;

		return false;
	}

	int subtractionFactor()
	{
		double signalRegion = abs(signalRegion_x2 - signalRegion_x1);
		double sidebandRegion = abs(sidebandRegion1_x2 - sidebandRegion1_x1) + abs(sidebandRegion2_x2 - sidebandRegion2_x1);

		return signalRegion/abs(sidebandRegion - signalRegion);
	}

	TBox* createTBox(double Ymax, int index = 0)
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

		return new TBox(x1, 0., x2, Ymax);;
	}
};

class InvariantMass{
private:
	int& method;
	const char*& ressonance;
	const char*& particleName;
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
		hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
	}

	void drawCanvasQuarter(TCanvas* &canvas, bool drawRegions, int quarter, MassValues* ObjMassValues, int color = kBlue)
	{
		TH1D* &hMass = ObjMassValues->hMass;
		TF1*  &fFit  = ObjMassValues->fitFunction;

		bool shouldDrawAllFitFunctions = true;

		float margins[4] = {0.13, 0.02, 0.09, 0.07};

		canvas->cd(quarter);
		canvas->cd(quarter)->SetMargin(margins[0], margins[1], margins[2], margins[3]);

		hMass->SetMarkerStyle(20);		//Set markers style
		hMass->SetMarkerColor(kBlack);	//Set markers colors
		hMass->SetLineColor(kBlack);	//Set errobars color
		hMass->Draw("ep");

		//Add legend
		TLegend* tl = new TLegend(0.70,0.80,0.96,0.92);
		tl->SetTextSize(0.04);
		tl->AddEntry(hMass, "Data", "lp");

		//Draw fit
		if (fFit != NULL)
		{
			fFit->SetNpx(1000);
			fFit->SetLineColor(color);
			fFit->SetLineStyle(kSolid);
			fFit->Draw("same");
			tl->AddEntry(fFit, "Total Fit", "l");

			//If is showing pass and fail fit
			if (quarter < 2 && shouldDrawAllFitFunctions == true)
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

			//Draw regions
			TBox* side1  = ObjMassValues->createTBox(Ymax, -1);
			side1->SetFillColorAlpha(kRed, 0.35);
			side1->Draw();
			TBox* signal = ObjMassValues->createTBox(Ymax, 0);
			signal->SetFillColorAlpha(kGreen, 0.35);
			signal->Draw();
			TBox* side2  = ObjMassValues->createTBox(Ymax, 1);
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
		tx->DrawLatex(0.16,0.88,Form("#bf{CMS Open Data}"));

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

	void fillMassHistograms(double* InvariantMass, int* isPassing)
	{
		if (*isPassing)
			this->Pass.hMass->Fill(*InvariantMass);
		this->All.hMass->Fill(*InvariantMass);
	}

	void doFit()
	{
		if (strcmp(ressonance, "Jpsi") == 0)
		{
			TH1D* &hPass 	 = this->Pass.hMass;
			TH1D* &hPassFail = this->All .hMass;

			//Get size of parNames
			int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);

			//Passing Fitting
			TF1* &fPass = this->Pass.fitFunction;
			fPass = new TF1("FitFunction_Pass", FitFunctions::Merged::Pass_InvariantMass, xMin, xMax, 12);
			for (int i = 0; i < arraySize; i++)
			{
				fPass->SetParName(i, this->fittingParName[i]);
			}

			//Both Fitting
			TF1* &fPassFail = this->All.fitFunction;
			fPassFail = new TF1("FitFunction_Both", FitFunctions::Merged::Both_InvariantMass, xMin, xMax, 24);
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

			//Set initial parameters (Monte Carlo)
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

			//Rename global parameters
			for (int i = 0; i < arraySize*2; i++)
			{
				fitter.Config().ParSettings(i).SetName(this->fittingParName[i%arraySize]);
			}

			//Fit FCN function directly
			//(specify optionally data size and flag to indicate that is a chi2 fit)
			fitter.FitFCN(24, globalChi2, 0, dataPass.Size() + dataPassFail.Size(), true);
			ROOT::Fit::FitResult result = fitter.Result();

			//Save signal fit
			this->Pass.fitSignal = new TF1("Signal_InvariantMass", FitFunctions::Merged::Signal_InvariantMass, xMin, xMax, 8);
			this->All .fitSignal = new TF1("Signal_InvariantMass", FitFunctions::Merged::Both_Signal_InvariantMass, xMin, xMax, 16);
			double fitParameters[24];
			fPass->GetParameters(fitParameters);
			Pass.fitSignal->SetParameters(fitParameters);
			All .fitSignal->SetParameters(fitParameters);

			//Show fit result
			cout << "For " << *particleType << "....";
			result.Print(std::cout);
			cout << endl;
		}

		if (strcmp(ressonance, "Upsilon") == 0)
		{
			/*
			double mass_peak1 = 9.46030;
			double mass_peak2 = 10.02326;
			double mass_peak3 = 10.3552;

			TCanvas* c_all  = new TCanvas;
			TCanvas* c_pass = new TCanvas;

			//DECLARE OBSERVABLE X - INVARIANT MASS BETWEEN _mmin AND _mmax
			RooRealVar mass_all("mass_all","mass_all", xMin, xMax);

			// Create a binned dataset that imports contents of TH1 and associates itscontents to observable 'mass'
			RooDataHist dh("dh","dh",mass_all,RooFit::Import(*this->All.hMass));
			RooDataHist dh_pass("dh_pass","dh_pass",mass_all,RooFit::Import(*this->Pass.hMass));


			RooRealVar lambda("lambda","lambda",-1.3,-10.,10.);
			RooExponential background("background", "background", mass_all, lambda);

			RooRealVar sigma("sigma","sigma",0.05*(xMax - xMin),0.,0.5*(xMax - xMin));

			RooRealVar mean1("mean1","mean1",xMin,(mass_peak1+mass_peak2)/2.);
			RooRealVar mean2("mean2","mean2",(mass_peak1+mass_peak2)/2.,(mass_peak3+mass_peak2)/2.);
			RooRealVar mean3("mean3","mean3",(mass_peak3+mass_peak2)/2.,xMax);
			//FIT FUNCTIONS

			// --Gaussian as the signal pdf
			RooGaussian gaussian1("signal1","signal1",mass_all,mean1,sigma);
			RooGaussian gaussian2("signal2","signal2",mass_all,mean2,sigma);
			RooGaussian gaussian3("signal3","signal3",mass_all,mean3,sigma);

			double n_signal_initial1 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak1)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak1,mass_peak1))) / dh.sumEntries();
			double n_signal_initial2 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak2)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak2,mass_peak2))) / dh.sumEntries();
			double n_signal_initial3 =(dh.sumEntries(TString::Format("abs(mass_all-%g)<0.015",mass_peak3)) -dh.sumEntries(TString::Format("abs(mass_all-%g)<0.030&&abs(mass_all-%g)>.015",mass_peak3,mass_peak3))) / dh.sumEntries();

			double n_signal_initial_total = n_signal_initial1 + n_signal_initial2 + n_signal_initial3;

			//signal = gaussian1*frac1 + gaussian2*frac2 + gaussian3*(1-(frac1 + frac2))
			//S(signal)d mass = 1
			RooRealVar frac1("frac1","frac1",0.333,0.,1.);
			RooRealVar frac2("frac2","frac2",0.333,0.,1.);
			//para RooArgList N-1, assume como frações
			RooAddPdf* signal;
			signal = new RooAddPdf("signal", "signal", RooArgList(gaussian1, gaussian2,gaussian3), RooArgList(frac1, frac2));

			double n_back_initial = 1. - n_signal_initial1 - n_signal_initial2 -n_signal_initial3;

			RooRealVar n_signal_total("n_signal_total","n_signal_total",n_signal_initial_total,0.,dh.sumEntries());

			RooRealVar n_back("n_back","n_back",n_back_initial,0.,dh.sumEntries());
			//para RooArgList N, assume como normalizações
			//modelo_total = n_signal_total*signal + n_back*background
			RooAddPdf* model;
			model = new RooAddPdf("model","model", RooArgList(*signal, background),RooArgList(n_signal_total, n_back));

			RooPlot *frame = mass_all.frame(RooFit::Title("Invariant Mass"));
			RooPlot *frame_new = mass_all.frame(RooFit::Title("Invariant Mass"));

			RooFitResult* fitres = new RooFitResult; //saves fit result
			fitres = model->fitTo(dh, RooFit::Save());

			// collecting fit results in order to differentiate between signal and background
			RooRealVar* pass_mean1 = (RooRealVar*) fitres->floatParsFinal().find("mean1");
			RooRealVar* pass_mean2 = (RooRealVar*) fitres->floatParsFinal().find("mean2");
			RooRealVar* pass_mean3 = (RooRealVar*) fitres->floatParsFinal().find("mean3");

			RooRealVar* pass_sigma = (RooRealVar*) fitres->floatParsFinal().find("sigma");

			double mean1_value = pass_mean1->getVal();
			double mean2_value = pass_mean2->getVal();
			double mean3_value = pass_mean3->getVal();

			double sigma_value = pass_sigma->getVal();
			*/
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

		//For regions of sideband subtraction
		ObjMassValues->M_JPSI = value;
		ObjMassValues->W_JPSI = sigma;
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
		const char*& directoryToSave,
	 	const char*& particleType)
		  : method(method),
		    ressonance(ressonance),
		    particleName(particleName),
		    directoryToSave(directoryToSave),
		    particleType(particleType)
	{
		if (strcmp(ressonance, "Jpsi") == 0)
		{
			this->xMin  = 2.8;
			this->xMax  = 3.4;
			this->nBins = 240;
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