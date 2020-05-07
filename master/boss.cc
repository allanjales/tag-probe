 /*
!--------------------------------
!Purpose: Find efficiency using Tag And Probe method
!--------------------------------	
!author: Allan Jales
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

#include <iostream>

#include "FitFunctions.h"

using namespace std;



class ParticleSelector{
public:
	const char *PassingOrFailing = "Unknown";

	double M_JPSI = 3.097;
	double W_JPSI = 0.010;

	double signalRegionEnd 		= 3.0;
	double sidebandRegionEnd	= 6.0;

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

	//Bins per each x axis unit (for integrations)
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

			"Exp1(Bg) Height  ",
			"Exp1(Bg) Width   ",
			"Exp2(Bg) Height  ",
			"Exp2(Bg) Width   "
		};

	Double_t resultParameters[12];

	ParticleSelector Selector;

	int color = kBlue;

	void define(const char *PassingOrFailing)
	{
		this->PassingOrFailing 			= PassingOrFailing;
		this->Selector.PassingOrFailing = PassingOrFailing;
	}

	void defineNumbers(int nBins, double xMin, double xMax, int decimals = 4)
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
		this->scale = 1/hMass->GetBinWidth(0);
	}

	void fill(double value)
	{
		this->hMass->Fill(value);
	}

	void fit()
	{
		//Create fit Function
		f = new TF1("FitFunction", FitFunctions::Merged::FFit_InvariantMassAll, xMin, xMax, 12);
		f->SetNpx(1000);
		f->SetLineStyle(kSolid);
		f->SetLineColor(color);

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

		//Fit the function
		cout << endl;
		cout << PassingOrFailing << " " << particle << endl;
		fitResult = hMass->Fit(f, "RNS", "", xMin, xMax);
		f->GetParameters(resultParameters);
		
		//Signal Fitting
		fs = new TF1("FitFunction_Signal", FitFunctions::Merged::Signal_InvariantMassAll, xMin, xMax, 8);
		fs->SetNpx(1000);							//Resolution of signal fit function
		fs->SetParameters(resultParameters);		//Get only signal part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Merged::Background_InvariantMassAll, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[8]);	//Get only background part
		fb->SetLineColor(color); 					//Fit Color
		fb->SetLineStyle(kDashed);					//Fit style

		//Get value and uncertain of signal
		double position = fs->GetMaximumX();
		double x1 = fs->GetX(fs->GetMaximum()/2);
		double x2 = fs->GetX(fs->GetMaximum()/2, x1+0.0001, position + x1*3);
		double fwhm = x2 - x1;
		double sigma = fwhm/2;

		this->Selector.M_JPSI = position;
		this->Selector.W_JPSI = sigma;
		this->Selector.signalRegionEnd = 3.;
		this->Selector.sidebandRegionEnd = 6.;
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
		if (drawRegions == true)
		{
			TBox *side1 = Selector.createTBox(Ymax, -1);
			side1->SetFillColorAlpha(kRed, 0.35);
			side1->Draw();
			TBox *signal = Selector.createTBox(Ymax, 0);
			signal->SetFillColorAlpha(kGreen, 0.35);
			signal->Draw();
			TBox *side2 = Selector.createTBox(Ymax, 1);
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

	bool isInSignalRegion(Double_t InvariantMass)
	{
		return this->Selector.isInSignalRegion(InvariantMass);
	}

	bool isInSidebandRegion(Double_t InvariantMass)
	{
		return this->Selector.isInSidebandRegion(InvariantMass);
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
private:
	void createHistogram(TH1D* &histo, const char *histoName)
	{
		//Define parameters
		string hName 		= PassingOrFailing + string(tagOrProbe) + string(particle) + "_" + string(quantityName) + string(histoName);
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
		histo = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
		histo->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), histo->GetBinWidth(0)));
		histo->GetXaxis()->SetTitle(xAxisTitle.data());
	}

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
		this->createHistogram(hSigBack, "SigBack");
	}

	void createBackHistogram()
	{
		this->createHistogram(hBack, "Back");
	}

	void createSigHistogram()
	{
		this->createHistogram(hSig, "Sig");
	}

	void subtractSigHistogram()
	{
		if (this->hSig == NULL)
		{
			//Create histogram
			this->createSigHistogram();;
		}
		else
		{
			cout << "WARNING! Sig Histogram already exists!" << endl;
		}

		this->hSig->Add(hSigBack,1);
		double factor = signalRegionEnd/abs(sidebandRegionEnd - signalRegionEnd);
		this->hSig->Add(hBack,-factor);
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

	void fillSigHistograms(Double_t PtValue, Double_t EtaValue, Double_t PhiValue)
	{
		this->Pt .hSig->Fill(PtValue);
		this->Eta.hSig->Fill(EtaValue);
		this->Phi.hSig->Fill(PhiValue);
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

	void createSigHistograms()
	{
		this->Pt .createSigHistogram();
		this->Eta.createSigHistogram();
		this->Phi.createSigHistogram();
	}

	void subtractSigHistograms()
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
		this->tagOrProbe = tagOrProbe;
	}
};

//Holder for TagProbe class
class PassingFailing{
private:
	int method = 0;
public:
	const char *PassingOrFailing = NULL;

	TagProbe Tag{"Tag"};
	TagProbe Probe{"Probe"};
	InvariantMassClass Mass;

	void define(const char *PassingOrFailing)
	{
		this->PassingOrFailing 	= PassingOrFailing;
		this->Mass.define(PassingOrFailing);
		this->Tag  .PassingOrFailing = PassingOrFailing;
		this->Probe.PassingOrFailing = PassingOrFailing;
	}

	void prepareSideband()
	{
		this->defineHistogramsTexts();
		this->defineHistogramsNumbers();
		this->createSigBackHistograms();
		this->createBackHistograms();
		this->Mass.createMassHistogram();
		this->method = 0;
	}

	void prepareFitting()
	{
		this->defineHistogramsTexts();
		this->defineHistogramsNumbers();
		this->createSigBackHistograms();
		this->createSigHistograms();
		this->createBackHistograms();
		this->Mass.createMassHistogram();
		this->method = 1;
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
		this->Mass.defineNumbers(240, 2.8, 3.4);
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

	void fillSigHistograms(Double_t TagPtValue, Double_t TagEtaValue, Double_t TagPhiValue, Double_t ProbePtValue, Double_t ProbeEtaValue, Double_t ProbePhiValue)
	{
		this->Tag  .fillSigHistograms(TagPtValue, 	TagEtaValue, 	TagPhiValue);
		this->Probe.fillSigHistograms(ProbePtValue, ProbeEtaValue, ProbePhiValue);
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

	void createSigHistograms()
	{
		this->Tag  .createSigHistograms();
		this->Probe.createSigHistograms();
	}

	void subtractSigHistograms()
	{
		this->Tag  .subtractSigHistograms();
		this->Probe.subtractSigHistograms();
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
		if (method == 0)
		{
			//If is inside signal region
			if (this->Mass.isInSignalRegion(InvariantMass))
			{
				this->fillSigBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			}

			//If is inside sideband region
			if (this->Mass.isInSidebandRegion(InvariantMass))
			{
				this->fillBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			}
		}
		else
		{
			//If is inside signal region
			if (this->Mass.isInSignalRegion(InvariantMass))
			{
				this->fillSigHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			}
			else
			{
				this->fillBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			}
			this->fillSigBackHistograms(TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}

	void updateSelectorParameters()
	{
		this->Tag  .setRange(Mass.Selector.signalRegionEnd, Mass.Selector.sidebandRegionEnd);
		this->Probe.setRange(Mass.Selector.signalRegionEnd, Mass.Selector.sidebandRegionEnd);
	}

	void debugCout()
	{
		this->Tag  .debugCout();
		this->Probe.debugCout();
	}

	PassingFailing(const char *PassingOrFailing)
	{
		this->define(PassingOrFailing);
	}
};


//Holder for PassingFailing class
class Particle
{
private:
	int method = 0;
	//if 1 -> sideband
	//if 2 -> fitting

	void prepareMethod()
	{
		if (this->method == 1)
		{
			this->PassingParticles.prepareSideband();
			this->FailingParticles.prepareSideband();
			this->AllParticles.prepareSideband();
			return;
		}
		
		if (this->method == 2)
		{
			this->PassingParticles.prepareFitting();
			this->FailingParticles.prepareFitting();
			this->AllParticles.prepareFitting();
			return;
		}

		cout << "WARNING! Invalid method" << endl;
	}

public:
	PassingFailing PassingParticles{"Passing"};
	PassingFailing FailingParticles{"Failing"};
	PassingFailing AllParticles{"All"};

	void setMethod(int method)
	{
		this->method = method;
		this->prepareMethod();
	}

	void doFit()
	{
		this->PassingParticles.Mass.fit();
		this->FailingParticles.Mass.fit();
		this->AllParticles.Mass.fit();
	}

	void updateSelectionParameters()
	{
		this->PassingParticles.updateSelectorParameters();
		this->FailingParticles.updateSelectorParameters();
		this->AllParticles.updateSelectorParameters();
	}

	void subtractSigHistogramsIfNeeded(bool ignoreCaution = false)
	{	
		if (this->method == 1)
		{
			this->PassingParticles.subtractSigHistograms();
			this->FailingParticles.subtractSigHistograms();
			this->AllParticles.subtractSigHistograms();
			return;
		}

		if (!ignoreCaution)
			cout << "CAUTION! Subtract histogram not needed";
	}

	//Debug methods

	void massDebugCout()
	{
		this->PassingParticles.Mass.debugCout();
		this->FailingParticles.Mass.debugCout();
		this->AllParticles.Mass.debugCout();
	}

	void consistencyDebugCout()
	{
		cout << endl;
		cout << "Checking histograms number inconsistency (should be 0)" << endl;
		this->PassingParticles.debugCout();
		this->FailingParticles.debugCout();
		this->AllParticles.debugCout();
	}
};

//Select particles, draws and save histograms
void generateHistograms(bool shouldDrawInvariantMassCanvas = true, bool shouldDrawQuantitiesCanvas = true, bool shouldDrawEfficiencyCanvas = true)
{
	int method = 0;
	//if 0 -> sideband
	//if 1 -> fitting

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
	Particle Muon;
	Muon.setMethod(1);

	//Loop between the components
	for (int i = 0; i < TreePC->GetEntries(); i++)
	{
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4 && !PassingProbeStandAloneMuon && !PassingProbeGlobalMuon)
		{
			if (PassingProbeTrackingMuon)
				Muon.PassingParticles.Mass.fill(InvariantMass);
			else
				Muon.FailingParticles.Mass.fill(InvariantMass);
			Muon.AllParticles.Mass.fill(InvariantMass);
		}
	}

	Muon.doFit();
	Muon.updateSelectionParameters();

	//Loop between the components again
	for (int i = 0; i < TreePC->GetEntries(); i++)
	{
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4 && !PassingProbeStandAloneMuon && !PassingProbeGlobalMuon)
		{
			if (PassingProbeTrackingMuon)
				Muon.PassingParticles.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			else
				Muon.FailingParticles.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			Muon.AllParticles.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}

	//Debug
	Muon.massDebugCout();

	Muon.subtractSigHistogramsIfNeeded(true);

	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Create file root to store generated files
	TFile *generatedFile = TFile::Open("../result/generated_hist.root","RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->   cd("canvas/");


	if (shouldDrawInvariantMassCanvas)
	{
		Muon.PassingParticles.Mass.createCanvas(true, true, true);
		Muon.FailingParticles.Mass.createCanvas(true, true, true);
		Muon.AllParticles	.Mass.createCanvas(true, true, true);
	}

	if (shouldDrawQuantitiesCanvas)
	{
		Muon.PassingParticles.createDividedCanvas(true, true);
		Muon.FailingParticles.createDividedCanvas(true, true);
		Muon.AllParticles	.createDividedCanvas(true, true);
	}

	//Debug
	Muon.consistencyDebugCout();

	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	//Write histograms on file
	Muon.PassingParticles.write(true, true, true);
	Muon.FailingParticles.write(true, true, true);
	Muon.AllParticles	.write(true, true, true);


	//Save plots
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	Muon.PassingParticles.createEfficiencyPlot();
	Muon.FailingParticles.createEfficiencyPlot();
	Muon.AllParticles	.createEfficiencyPlot();


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");

		Muon.PassingParticles.createEfficiencyCanvas();
	if (shouldDrawEfficiencyCanvas)
	{
		Muon.PassingParticles.createEfficiencyCanvas();
		Muon.FailingParticles.createEfficiencyCanvas();
		Muon.AllParticles	.createEfficiencyCanvas();
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
	generateHistograms(true, false, false);
}