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

#include <iostream>

using namespace std;



//-------------------------------------
// General functions
//-------------------------------------

//Gaussian function
Double_t Gaus(Double_t *x, Double_t *par)
{
	//par[0] = height
	//par[1] = position
	//par[2] = sigma
	Double_t gaus = par[0]*TMath::Gaus(x[0],par[1],par[2],1);
	return gaus;
}

//Polynomial function
Double_t Pol1(Double_t *x, Double_t *par)
{
	//par[0] = b
	//par[1] = a
	Double_t pol = par[0] + par[1]*x[0];
	return pol;
}

//Exponential function
Double_t Exp(Double_t *x, Double_t *par)
{
	//par[0] = height
	//par[1] = width
	Double_t exp = par[0] * TMath::Exp(par[1]*x[0]);
	return exp;
}

//crystall ball function
Double_t CrystalBall(Double_t *x,Double_t *par)
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



//-------------------------------------
// Fit functions
//-------------------------------------

//Fit function for signal for Invariant Mass Probe
Double_t Signal_InvariantMassProbe(Double_t *x, Double_t *par) {
	return Gaus(x,par) + CrystalBall(x, &par[3]);
}

//Fit function for background for Invariant Mass Probe
Double_t Background_InvariantMassProbe(Double_t *x, Double_t *par) {
	return Exp(x,par) + Exp(x, &par[2]);
}

//Fit function for signal & background for Invariant Mass Probe
Double_t FFit_InvariantMassProbe(Double_t *x, Double_t *par) {
	return Signal_InvariantMassProbe(x,par) + Background_InvariantMassProbe(x, &par[8]);
}

//Create canvas and fitting for Invariant Mass
TCanvas *invariantMassProbe(TH1D *hMassAll, double S, double dS, bool shouldWrite = false, const char *saveAs = "")
{
	//Create canvas
	TCanvas *c1 = new TCanvas("AllInvariantMass","Invariant Mass", 600, 600);

	//Set margin for canvas
	c1->SetTopMargin(0.07);
	c1->SetRightMargin(0.05);
	c1->SetBottomMargin(0.11);
	c1->SetLeftMargin(0.15);

	//Set title size for axis and other stuffs for histogram style
	hMassAll->GetYaxis()->SetTitleSize(0.04);
	hMassAll->GetXaxis()->SetTitleSize(0.05);
	hMassAll->GetXaxis()->CenterTitle(true);
	//hMassAll->SetBit(TH1::kNoTitle);			//Not show histogram title
	hMassAll->SetMarkerStyle(20);				//Set markers style
	hMassAll->SetMarkerColor(kBlack);			//Set markers colors
	hMassAll->SetLineColor(kBlack);				//Set lines colors (for errorbars)

	//Fit Function
	TF1 *f = new TF1("f",FFit_InvariantMassProbe,hMassAll->GetXaxis()->GetXmin(),hMassAll->GetXaxis()->GetXmax(),12);
	f->SetParName(0,	"Gaus(Sg) Height");
	f->SetParName(1,	"Gaus(Sg) Position");
	f->SetParName(2,	"Gaus(Sg) Sigma");
	f->SetParName(3,	"CB(Sg) Alpha");
	f->SetParName(4,	"CB(Sg) N");
	f->SetParName(5,	"CB(Sg) Mean");
	f->SetParName(6,	"CB(Sg) Sigma");
	f->SetParName(7,	"CB(Sg) Yield");
	f->SetParName(8,	"Exp1(Bg) Lambda");
	f->SetParName(9,	"Exp1(Bg) a");
	f->SetParName(10,	"Exp2(Bg) Lambda");
	f->SetParName(11,	"Exp2(Bg) a");
	f->SetNpx(1000);	//Resolution of fit function
	
	//Values Signal
	//f->SetParameter(0,	2000.0);
	f->SetParameter(1,	3.1);
	f->SetParameter(2,	0.021);
	f->SetParameter(3,	2.1);
	f->SetParameter(4,	0.40);
	f->SetParameter(5,	3.094);
	f->SetParameter(6,	0.02);
	//f->SetParameter(7,	0);

	//Values Background
	f->SetParameter(8,	-1.7);
	f->SetParameter(9,	2.0);
	//f->SetParameter(10, -1.7);
	//f->SetParameter(11,	2.0);

	//Fit Color
	f->SetLineColor(kBlue);

	//Fit the function
	TFitResultPtr fitr = hMassAll->Fit(f,"RNS","",hMassAll->GetXaxis()->GetXmin(),hMassAll->GetXaxis()->GetXmax());
	//Get parameters from fit function and put it in par variable
   	Double_t par[12];
	f->GetParameters(par);

	//Signal Fitting
	TF1 *fs = new TF1("fs",Signal_InvariantMassProbe,hMassAll->GetXaxis()->GetXmin(),hMassAll->GetXaxis()->GetXmax(),8);
	fs->SetNpx(1000);				//Resolution of background fit function
	fs->SetParameters(par);			//Get only background part
	fs->SetLineColor(kMagenta); 	//Fit Color
	fs->SetLineStyle(kSolid);		//Fit Style

	//Background Fitting
	TF1 *fb = new TF1("fb",Background_InvariantMassProbe,hMassAll->GetXaxis()->GetXmin(),hMassAll->GetXaxis()->GetXmax(),4);
	fb->SetNpx(1000);				//Resolution of background fit function
	fb->SetParameters(&par[8]);		//Get only background part
	fb->SetLineColor(kBlue); 		//Fit Color
	fb->SetLineStyle(kDashed);		//Fit Style

	//Draws histogram & fit function
	hMassAll->Draw("ep");
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
	//tx->DrawLatex(0.6,0.60,Form("#chi^{2}/ndf = %g/%d",fitr->Chi2(),fitr->Ndf()));
	tx->DrawLatex(0.61,0.60,Form("#chi^{2}/ndf = %.3g",fitr->Chi2()/fitr->Ndf()));

	//Show number of particles
	tx->DrawLatex(0.61,0.48,Form("%.0f #pm %.0f J/#psi", S, dS));
	tx->DrawLatex(0.64,0.43,"(candidates)");

	//Add legend
	TLegend *l = new TLegend(0.65,0.75,0.92,0.90);
	l->SetTextSize(0.04);
	l->AddEntry(hMassAll,	"J/#psi","lp");
	l->AddEntry(f,			"Fitting","l");
	l->AddEntry(fs,			"Signal","l");
	l->AddEntry(fb,			"Background","l");
	l->Draw();

	//Not show frame with mean, std dev
	gStyle->SetOptStat(0);
	
	//Show chi-squared test (on prompt)
	cout << endl;
	cout << "Fitting overview" << endl;
	cout << "Chi2/ndf = " << f->GetChisquare()/f->GetNDF() << endl;
	printf( "#Signal  = %.0f +- %.0f\n", S, dS);

	//Integrate function to get number of particles in it
	cout << endl;
	cout << "Candidates by integration" << endl;
	cout << "#Total      = " << f ->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;
	cout << "#Background = " << fb->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;
	cout << "#Signal     = " << fs->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;

	//Writes in file
	if (shouldWrite == true)
	{
		//Aqui dá crash quando roda pela segunda vez no meu ROOT
		c1->Write("",TObject::kOverwrite);
	}

	//If should save
	if (strcmp(saveAs, "") != 0)
	{
		//Saves as image
		c1->SaveAs(saveAs);
	} 

	//return
	return c1;
}



TCanvas *createDividedCanvas(TH1D *hSigBack, TH1D *hSig, TH1D *hBack, const char *canvasName, const char *titleLeft, const char *titleRight, bool shouldWrite = false, const char *saveAs = "")
{
	//Create canvas
	TCanvas *c2 = new TCanvas(canvasName, titleLeft, 1200, 600);

	//Divide canvas
	c2->Divide(2,1);

	//Select canvas part
	c2->cd(1);

	//Set margin for canvas part
	c2->cd(1)->SetTopMargin(0.07);
	c2->cd(1)->SetRightMargin(0.02);
	c2->cd(1)->SetBottomMargin(0.11);
	c2->cd(1)->SetLeftMargin(0.14);

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

	/*
	//Set width of line
	hSigBack->SetLineWidth(2);
	hSig->SetLineWidth(2);
	hBack->SetLineWidth(2);
	*/

	/*
	//Set fill color
	hSig->SetFillColor(kMagenta);
	hSigBack->SetFillColor(kBlue);
	hBack->SetFillColor(kYellow);
	*/

	//Add legend
	TLegend *l2_1 = new TLegend(0.65,0.75,0.92,0.90);
	l2_1->SetTextSize(0.04);
	l2_1->AddEntry(hSigBack,	"All",			"lp");
	l2_1->AddEntry(hSig,		"Signal",		"l");
	l2_1->AddEntry(hBack,		"Background",	"l");
	l2_1->Draw();
	
	//Draws information
	TLatex *tx2_1 = new TLatex();
	tx2_1->SetTextSize(0.04);
	tx2_1->SetTextFont(42);
	tx2_1->SetNDC(kTRUE);

	if (strcmp(canvasName, "ProbeSignal_Pt") == 0 || strcmp(canvasName, "TagSignal_Pt") == 0) 
	{
		tx2_1->SetTextAlign(12);	//Align left, center
		tx2_1->DrawLatex(0.48,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
		tx2_1->DrawLatex(0.48,0.45,Form("%g entries (signal)",		hSig->GetEntries()));
		tx2_1->DrawLatex(0.48,0.40,Form("%g entries (background)",	hBack->GetEntries()));
	}
	else
	{
		tx2_1->SetTextAlign(22);	//Align center, center
		tx2_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hSigBack->GetEntries()));
		tx2_1->DrawLatex(0.55,0.45,Form("%g entries (signal)",		hSig->GetEntries()));
		tx2_1->DrawLatex(0.55,0.40,Form("%g entries (background)",	hBack->GetEntries()));
	}
	

	//Select canvas part
	c2->cd(2);

	//Set margin for canvas part
	c2->cd(2)->SetTopMargin(0.07);
	c2->cd(2)->SetRightMargin(0.02);
	c2->cd(2)->SetBottomMargin(0.11);
	c2->cd(2)->SetLeftMargin(0.13);

	//Same range as comparision
	//hPtSig->GetYaxis()->SetRangeUser(0,hPtSigBack->GetYaxis()->GetYmax());
	hSig->SetMinimum(0);
	hSig->SetTitle(titleRight);

	//Signal
	hSig->Draw("same");

	//Add legend
	TLegend *l2_2 = new TLegend(0.65,0.85,0.92,0.90);
	l2_2->SetTextSize(0.04);
	l2_2->AddEntry(hSig, "Signal","l");
	l2_2->Draw();

	//Draws information
	TLatex *tx2_2 = new TLatex();
	tx2_2->SetTextSize(0.04);
	tx2_2->SetTextFont(42);
	tx2_2->SetNDC(kTRUE);

	//Not show frame with mean, std dev
	gStyle->SetOptStat(0);
	
	if (strcmp(canvasName, "ProbeSignal_Pt") == 0 || strcmp(canvasName, "TagSignal_Pt") == 0)
	{
		tx2_2->SetTextAlign(12);	//Align left, center
		tx2_2->DrawLatex(0.48,0.50,Form("%g entries (signal)", hSig->GetEntries()));
	}
	else
	{
		tx2_2->SetTextAlign(22);	//Align center, center
		tx2_2->DrawLatex(0.55,0.5,Form("%g entries (signal)",hSig->GetEntries()));
	}

	//Writes in file
	if (shouldWrite == true)
	{
		c2->Write();
	}

	//If should save
	if (strcmp(saveAs, "") != 0)
	{
		//Saves as image
		c2->SaveAs(saveAs);
	} 

	//return
	return c2;
}

//-------------------------------------
// Main functions
//-------------------------------------

//Draws and save invariant mass histogram for pp
void step1()
{
	TFile *file0 = TFile::Open("../data_histoall.root");		//Opens the file
	TTree *TreePC = (TTree*)file0->Get("demo/PlotControl");		//Opens TTree of file
	TTree *TreeAT = (TTree*)file0->Get("demo/AnalysisTree");	//Opens TTree of file
	
	//Create variables for PlotControl
	Double_t 	ProbeMuon_Pt;
	Double_t 	ProbeMuon_Eta;
	Double_t 	ProbeMuon_Phi;
	Double_t 	TagMuon_Pt;
	Double_t 	TagMuon_Eta;
	Double_t 	TagMuon_Phi;
	Double_t 	InvariantMass;

	//Create variables for AnalysisTree
	int PassingProbeTrackingMuon;
	int PassingProbeStandAloneMuon;
	int PassingProbeGlobalMuon;

	//Jpsi constants (PDG values)
	const double M_JPSI = 3.097;
	const double W_JPSI = 0.010;

	//Variables for side band subtraction
	int count_sigregion = 0;	//Inside sigma
	int count_sideband = 0;		//Outside sigma

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

	//Create histograms for Invariant Mass
	TH1D *hMassAll    = new TH1D("ProbeMuon_InvariantMass", "Invariant Mass (All);Mass (GeV/c^{2});Events",240,2.8,3.4);
	hMassAll->GetYaxis()->SetTitle(Form("Events / (%1.4f GeV/c^{2})", hMassAll->GetBinWidth(0)));

	//Create histograms for Signal + Background & Background only for Probe Pt
	TH1D *hPtSigBack  = new TH1D("ProbeMuon_PtSigBack",	"Transversal Momentum (Probe);P_{t} (GeV/c);Events",100,0.,100.);
	hPtSigBack->GetYaxis()->SetTitle(Form("Events / (%1.1f GeV/c)", hPtSigBack->GetBinWidth(0)));
	TH1D* hPtBack  = (TH1D*) hPtSigBack->Clone("hPtBack");
	hPtBack->SetTitle("Transversal Momentum (Probe)");
	hPtBack->SetName("ProbeMuon_PtBack");

	//Create histograms for Signal + Background & Background only for Probe Eta
	TH1D *hEtaSigBack = new TH1D("ProbeMuon_EtaSigBack", "Pseudorapidity (Probe);#eta;Events",200,-2.5,2.5);
	hEtaSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", hEtaSigBack->GetBinWidth(0)));
	TH1D* hEtaBack = (TH1D*) hEtaSigBack->Clone("hEtaBack");
	hEtaBack->SetTitle("Pseudorapidity (Probe)");
	hEtaBack->SetName("ProbeMuon_EtaBack");

	//Create histograms for Signal + Background & Background only for Probe Phi
	TH1D *hPhiSigBack = new TH1D("ProbeMuon_PhiSigBack", "Phi (Probe);#phi;Events",79,-3.15,3.15);
	hPhiSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", hPhiSigBack->GetBinWidth(0)));
	TH1D* hPhiBack = (TH1D*) hPhiSigBack->Clone("hPhiBack");
	hPhiBack->SetTitle("Phi (Probe)");
	hPhiBack->SetName("ProbeMuon_PhiBack");

	//Create histograms for Signal + Background & Background only for Tag Pt
	TH1D *hTagPtSigBack  = new TH1D("TagMuon_PtSigBack",	"Transversal Momentum (Tag);P_{t} (GeV/c);Events",100,0.,100.);
	hTagPtSigBack->GetYaxis()->SetTitle(Form("Events / (%1.1f GeV/c)", hTagPtSigBack->GetBinWidth(0)));
	TH1D* hTagPtBack  = (TH1D*) hTagPtSigBack->Clone("hTagPtBack");
	hTagPtBack->SetTitle("Transversal Momentum (Tag)");
	hTagPtBack->SetName("TagMuon_PtBack");

	//Create histograms for Signal + Background & Background only for Tag Eta
	TH1D *hTagEtaSigBack = new TH1D("TagMuon_EtaSigBack", "Pseudorapidity (Tag);#eta;Events",200,-2.5,2.5);
	hEtaSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", hTagEtaSigBack->GetBinWidth(0)));
	TH1D* hTagEtaBack = (TH1D*) hTagEtaSigBack->Clone("hTagEtaBack");
	hTagEtaBack->SetTitle("Pseudorapidity (Tag)");
	hTagEtaBack->SetName("TagMuon_EtaBack");

	//Create histograms for Signal + Background & Background only for Tag Phi
	TH1D *hTagPhiSigBack = new TH1D("TagMuon_PhiSigBack", "Phi (Tag);#phi;Events",79,-3.15,3.15);
	hTagPhiSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", hTagPhiSigBack->GetBinWidth(0)));
	TH1D* hTagPhiBack = (TH1D*) hTagPhiSigBack->Clone("hTagPhiBack");
	hTagPhiBack->SetTitle("Phi (Tag)");
	hTagPhiBack->SetName("TagMuon_PhiBack");

	//Loop between the components
	for (int i = 0; i < TreePC->GetEntries(); i++)
	{
		//Gets the entry from TTree
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Fill invariant mass histogram
		hMassAll->Fill(InvariantMass);

		//if is inside signal region
		if (fabs(InvariantMass - M_JPSI) < W_JPSI*3.0)
		{
			if (PassingProbeTrackingMuon && !PassingProbeStandAloneMuon && !PassingProbeGlobalMuon && abs(TagMuon_Pt) > 7.0 && abs(TagMuon_Eta) < 2.4)
			{
				//Add to Probe histogram
				hPtSigBack->Fill(ProbeMuon_Pt);		//Add to Pt  histogram
				hEtaSigBack->Fill(ProbeMuon_Eta);	//Add to Eta histogram
				hPhiSigBack->Fill(ProbeMuon_Phi);	//Add to Phi histogram

				//Add to Tag histogram
				hTagPtSigBack->Fill(TagMuon_Pt);	//Add to Pt  histogram
				hTagEtaSigBack->Fill(TagMuon_Eta);	//Add to Eta histogram
				hTagPhiSigBack->Fill(TagMuon_Phi);	//Add to Phi histogram
				
				//Count events
				count_sigregion++;
			}
		}

		//If is inside sideband region
		if (fabs(InvariantMass - M_JPSI) > W_JPSI*3.5 && fabs(InvariantMass - M_JPSI) < W_JPSI*6.5)
		{
			//Add to Probe histogram
			hPtBack->Fill(ProbeMuon_Pt);	//Add to Pt  histogram
			hEtaBack->Fill(ProbeMuon_Eta);	//Add to Eta histogram
			hPhiBack->Fill(ProbeMuon_Phi);	//Add to Phi histogram

			//Add to Tag histogram
			hTagPtBack->Fill(TagMuon_Pt);	//Add to Pt  histogram
			hTagEtaBack->Fill(TagMuon_Eta);	//Add to Eta histogram
			hTagPhiBack->Fill(TagMuon_Phi);	//Add to Phi histogram

			//Count events
			count_sideband++;
		}
	}

	//Number of particles in signal and uncertain
	double S = count_sigregion - count_sideband;
	double dS = sqrt(count_sigregion+count_sideband);

	//Create Signal histogram for Probe Pt
	TH1D* hPtSig  = (TH1D*) hPtSigBack->Clone("hPtSig");
	hPtSig->SetName("ProbeMuon_PtSig");
	hPtSig->Add(hPtBack,-1);

	//Create Signal histogram for Probe Eta
	TH1D* hEtaSig = (TH1D*) hEtaSigBack->Clone("hEtaSig");
	hEtaSig->SetName("ProbeMuon_EtaSig");
	hEtaSig->Add(hEtaBack,-1);

	//Create Signal histogram for Probe Phi
	TH1D* hPhiSig = (TH1D*) hPhiSigBack->Clone("hPhiSig");
	hPhiSig->SetName("ProbeMuon_PhiSig");
	hPhiSig->Add(hPhiBack,-1);

	//Create Signal histogram for Tag Pt
	TH1D* hTagPtSig  = (TH1D*) hTagPtSigBack->Clone("hTagPtSig");
	hTagPtSig->SetName("TagMuon_PtSig");
	hTagPtSig->Add(hTagPtBack,-1);

	//Create Signal histogram for Tag Eta
	TH1D* hTagEtaSig = (TH1D*) hTagEtaSigBack->Clone("hTagEtaSig");
	hTagEtaSig->SetName("TagMuon_EtaSig");
	hTagEtaSig->Add(hTagEtaBack,-1);

	//Create Signal histogram for Tag Phi
	TH1D* hTagPhiSig = (TH1D*) hTagPhiSigBack->Clone("hTagPhiSig");
	hTagPhiSig->SetName("TagMuon_PhiSig");
	hTagPhiSig->Add(hTagPhiBack,-1);

	//-------------------------------------
	// Canvas
	//-------------------------------------

	//Create file root to store generated files
	TFile *generatedFile = TFile::Open("../generated_hist.root","RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->cd("canvas/");

	//Create canvas and fitting
	TCanvas *c1 = invariantMassProbe(hMassAll, S, dS, true, "../InvariantMassAll.png");

	//Debug
	cout << endl;
	cout << "Candidates by sideband subtraction" << endl;
	cout << "#Tree Entries    = " 	<< TreePC->GetEntries() 			<< endl;
	cout << "#Signal   region = "	<< count_sigregion 					<< endl;
	cout << "#Sideband region = "	<< count_sideband 					<< endl;
	cout << "#Signal          = " 	<< count_sigregion - count_sideband	<< endl;
	cout << endl;

	//Create canvas for others
	createDividedCanvas(hPtSigBack,  	hPtSig,  	hPtBack,  	 "ProbeSignal_Pt",  "Probe Transversal Momentum", "Transversal Momentum of Signal (Probe)", true, "../PtProbe.png");
	createDividedCanvas(hEtaSigBack, 	hEtaSig, 	hEtaBack, 	 "ProbeSignal_Eta", "Probe Pseudorapidity", 		 "Pseudorapidity of Signal (Probe)", 	true, "../EtaProbe.png");
	createDividedCanvas(hPhiSigBack, 	hPhiSig, 	hPhiBack, 	 "ProbeSignal_Phi", "Probe Angle", 				 "Angle of Signal (Probe)", 				true, "../PhiProbe.png");
	createDividedCanvas(hTagPtSigBack,  hTagPtSig,  hTagPtBack,  "TagSignal_Pt",    "Tag Transversal Momentum",	 "Transversal Momentum of Signal (Tag)",    true, "../PtTag.png");
	createDividedCanvas(hTagEtaSigBack, hTagEtaSig, hTagEtaBack, "TagSignal_Eta",   "Tag Pseudorapidity", 		 "Pseudorapidity of Signal (Tag)", 		    true, "../EtaTag.png");
	createDividedCanvas(hTagPhiSigBack, hTagPhiSig, hTagPhiBack, "TagSignal_Phi",   "Tag Angle", 				 "Angle of Signal (Tag)", 				    true, "../PhiTag.png");

	//Integrate function to get number of particles in it
	cout << endl;
	cout << "Checking histograms number inconsistency (should be 0)" << endl;
	cout << "#Probe Pt  = " << hPtSigBack->GetEntries()  	- hPtSig->GetEntries()     - hPtBack->GetEntries()		<< endl;
	cout << "#Probe Eta = " << hEtaSigBack->GetEntries() 	- hEtaSig->GetEntries()    - hEtaBack->GetEntries()		<< endl;
	cout << "#Probe Phi = " << hPhiSigBack->GetEntries() 	- hPhiSig->GetEntries()    - hPhiBack->GetEntries()		<< endl;
	cout << "#Tag   Pt  = " << hTagPtSigBack->GetEntries()  - hTagPtSig->GetEntries()  - hTagPtBack->GetEntries()  	<< endl;
	cout << "#Tag   Eta = " << hTagEtaSigBack->GetEntries() - hTagEtaSig->GetEntries() - hTagEtaBack->GetEntries() 	<< endl;
	cout << "#Tag   Phi = " << hTagPhiSigBack->GetEntries() - hTagPhiSig->GetEntries() - hTagPhiBack->GetEntries() 	<< endl;

	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->cd("histograms/");
	//Save Probe histograms
	hPtSigBack->Write();
	hPtBack->Write();
	hPtSig->Write();
	hEtaSigBack->Write();
	hEtaBack->Write();
	hEtaSig->Write();
	hPhiSigBack->Write();
	hPhiBack->Write();
	hPhiSig->Write();
	//Save Tag histograms
	hTagPtSigBack->Write();
	hTagPtBack->Write();
	hTagPtSig->Write();
	hTagEtaSigBack->Write();
	hTagEtaBack->Write();
	hTagEtaSig->Write();
	hTagPhiSigBack->Write();
	hTagPhiBack->Write();
	hTagPhiSig->Write();
	//Close files
	generatedFile->Close();
}

//Estimates efficiency
void efficiencyOldMethod()
{
	//Opens the file
	TFile *generatedFile = TFile::Open("../generated_hist.root","UPDATE");

	//Import histograms
	TH1D *hPtSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PtSigBack");
	TH1D *hPtSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PtSig");
	TH1D *hEtaSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_EtaSigBack");
	TH1D *hEtaSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_EtaSig");
	TH1D *hPhiSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PhiSigBack");
	TH1D *hPhiSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PhiSig");

	//Efficiency calculation for Pt
	TH1D* hPtEff  = (TH1D*) hPtSig->Clone("ProbePt_Efficiency");
	hPtEff->SetTitle("Transversal Momentum Efficiency for Probe");
	hPtEff->GetYaxis()->SetTitle("Efficiency");
	hPtEff->Divide(hPtSig, hPtSigBack, 1.0, 1.0, "B");
	hPtEff->SetLineWidth(2);
	hPtEff->SetLineColor(kRed);
	hPtEff->SetMarkerStyle(21);
	hPtEff->SetMarkerSize(0.5);
	hPtEff->SetMarkerColor(kRed);
	hPtEff->SetMaximum(1);
	hPtEff->SetMinimum(0);
	TCanvas *c1 = new TCanvas("ProbePt_Efficiency","Probe Pt Efficiency", 800, 600);
	c1->SetTopMargin(0.07);
	c1->SetLeftMargin(0.12);
	c1->SetTicky(2);
	hPtEff->Draw();

	//Efficiency calculation for Eta
	TH1D* hEtaEff  = (TH1D*) hEtaSig->Clone("ProbeEta_Efficiency");
	hEtaEff->SetTitle("Pseudorapidity Efficiency for Probe");
	hEtaEff->GetYaxis()->SetTitle("Efficiency");
	hEtaEff->Divide(hEtaSig, hEtaSigBack, 1.0, 1.0, "B");
	hEtaEff->SetLineWidth(2);
	hEtaEff->SetLineColor(kRed);
	hEtaEff->SetMarkerStyle(21);
	hEtaEff->SetMarkerSize(0.5);
	hEtaEff->SetMarkerColor(kRed);
	hEtaEff->SetMaximum(1);
	hEtaEff->SetMinimum(0);
	TCanvas *c2 = new TCanvas("ProbeEta_Efficiency","Probe Eta Efficiency", 800, 600);
	c2->SetTopMargin(0.07);
	c2->SetLeftMargin(0.12);
	c2->SetTicky(2);
	hEtaEff->Draw();

	//Efficiency calculation for Phi
	TH1D* hPhiEff  = (TH1D*) hPhiSig->Clone("ProbePhi_Efficiency");
	hPhiEff->SetTitle("Phi Efficiency for Probe");
	hPhiEff->GetYaxis()->SetTitle("Efficiency");
	hPhiEff->Divide(hPhiSig, hPhiSigBack, 1.0, 1.0, "B");
	hPhiEff->SetLineWidth(2);
	hPhiEff->SetLineColor(kRed);
	hPhiEff->SetMarkerStyle(21);
	hPhiEff->SetMarkerSize(0.5);
	hPhiEff->SetMarkerColor(kRed);
	hPhiEff->SetMaximum(1);
	hPhiEff->SetMinimum(0);
	TCanvas *c3 = new TCanvas("ProbePhi_Efficiency","Probe Phi Efficiency", 800, 600);
	c3->SetTopMargin(0.07);
	c3->SetLeftMargin(0.12);
	c3->SetTicky(2);
	hPhiEff->Draw();



	//Not show frame
	gStyle->SetOptStat(0);

	//Saves as image
	c1->SaveAs("../PtProbe_Efficiency.png");
	c2->SaveAs("../EtaProbe_Efficiency.png");
	c3->SaveAs("../PhiProbe_Efficiency.png");

	//Saves new histograms and canvas in file
	generatedFile->Delete("efficiency/");
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");
	c1->Write("",TObject::kOverwrite);
	c2->Write("",TObject::kOverwrite);
	c3->Write("",TObject::kOverwrite);
	generatedFile->mkdir("efficiency/histograms/");
	generatedFile->cd("efficiency/histograms/");
	hPtEff->Write("",TObject::kOverwrite);
	hEtaEff->Write("",TObject::kOverwrite);
	hPhiEff->Write("",TObject::kOverwrite);
}

//Creates a efficiency plot with histograms
TEfficiency *efficiencyPlot(TH1D *hPass, TH1D *hTotal, const char *name, const char *title, bool shouldWrite = false)
{
	//Creates TEfficiency object
	TEfficiency* pEff = 0;

	//Set Y axis title for efficiency plot
	hTotal->GetYaxis()->SetTitle("Efficiency");

	//Check if are valid and consistent histograms
	if(TEfficiency::CheckConsistency(*hPass, *hTotal))
	{
		//Fills histogram
		pEff = new TEfficiency(*hPass, *hTotal);
	}

	//Set plot config
	pEff->SetTitle(title);
	pEff->SetLineWidth(2);
	pEff->SetLineColor(kRed);
	pEff->SetMarkerStyle(21);
	pEff->SetMarkerSize(0.5);
	pEff->SetMarkerColor(kRed);

	//Writes in file
	if (shouldWrite == true)
	{
		pEff->Write("",TObject::kOverwrite);
	}

	//return
	return pEff;
}

//Creates canvas for efficiency plots
TCanvas *createEfficiencyCanvas(TEfficiency* pEff, const char *canvasName, const char *title, bool shouldWrite = false, const char *saveAs = "")
{
	//Draw on canvas
	TCanvas *c1 = new TCanvas(canvasName, title, 800, 600);
	c1->SetRightMargin(0.05);
	pEff->Draw();

	//Writes in file
	if (shouldWrite == true)
	{
		c1->Write("",TObject::kOverwrite);
	}

	//If should save
	if (strcmp(saveAs, "") != 0)
	{
		//Saves as image
		c1->SaveAs(saveAs);
	} 

	//return
	return c1;
}

//Estimates efficiency
void efficiency()
{
	//Opens the file
	TFile *generatedFile = TFile::Open("../generated_hist.root","UPDATE");

	//Import Probe histograms
	TH1D *hProbePtSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PtSigBack");
	TH1D *hProbePtSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PtSig");
	TH1D *hProbeEtaSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_EtaSigBack");
	TH1D *hProbeEtaSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_EtaSig");
	TH1D *hProbePhiSigBack 	= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PhiSigBack");
	TH1D *hProbePhiSig 		= (TH1D*)generatedFile->Get("histograms/ProbeMuon_PhiSig");
	//Import Tag histograms
	TH1D *hTagPtSigBack 	= (TH1D*)generatedFile->Get("histograms/TagMuon_PtSigBack");
	TH1D *hTagPtSig 		= (TH1D*)generatedFile->Get("histograms/TagMuon_PtSig");
	TH1D *hTagEtaSigBack 	= (TH1D*)generatedFile->Get("histograms/TagMuon_EtaSigBack");
	TH1D *hTagEtaSig 		= (TH1D*)generatedFile->Get("histograms/TagMuon_EtaSig");
	TH1D *hTagPhiSigBack 	= (TH1D*)generatedFile->Get("histograms/TagMuon_PhiSigBack");
	TH1D *hTagPhiSig 		= (TH1D*)generatedFile->Get("histograms/TagMuon_PhiSig");

	//Deletes old dir and creates another
	generatedFile->Delete("efficiency/");
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	TEfficiency* pProbePtEff  = efficiencyPlot(hProbePtSig,  hProbePtSigBack,  	"ProbeMuon_PtEfficiency",  "Transversal Momentum Efficiency for Probe", true);
	TEfficiency* pProbeEtaEff = efficiencyPlot(hProbeEtaSig, hProbeEtaSigBack, 	"ProbeMuon_EtaEfficiency", "Pseudorapidity Efficiency for Probe", 		true);
	TEfficiency* pProbePhiEff = efficiencyPlot(hProbePhiSig, hProbePhiSigBack, 	"ProbeMuon_PhiEfficiency", "Angle Efficiency for Probe", 				true);
	TEfficiency* pTagPtEff    = efficiencyPlot(hTagPtSig,  	 hTagPtSigBack,  	"TagMuon_PtEfficiency",    "Transversal Momentum Efficiency for Tag", 	true);
	TEfficiency* pTagEtaEff   = efficiencyPlot(hTagEtaSig, 	 hTagEtaSigBack, 	"TagMuon_EtaEfficiency",   "Pseudorapidity Efficiency for Tag", 		true);
	TEfficiency* pTagPhiEff   = efficiencyPlot(hTagPhiSig, 	 hTagPhiSigBack, 	"TagMuon_PhiEfficiency",   "Angle Efficiency for Tag",					true);

	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/histograms/");
	generatedFile->cd("efficiency/histograms/");

	//Create canvas for others
	createEfficiencyCanvas(pProbePtEff,   "ProbeMuon_PtEfficiency",  "Transversal Momentum Efficiency for Probe", 	true,	"../PtProbe_Efficiency.png");
	createEfficiencyCanvas(pProbeEtaEff,  "ProbeMuon_EtaEfficiency", "Pseudorapidity Efficiency for Probe", 		true,	"../EtaProbe_Efficiency.png");
	createEfficiencyCanvas(pProbePhiEff,  "ProbeMuon_PhiEfficiency", "Angle Efficiency for Probe", 					true,	"../PhiProbe_Efficiency.png");
	createEfficiencyCanvas(pTagPtEff,     "TagMuon_PtEfficiency", 	 "Transversal Momentum Efficiency for Tag", 	true,	"../PtTag_Efficiency.png");
	createEfficiencyCanvas(pTagEtaEff,    "TagMuon_EtaEfficiency",	 "Pseudorapidity Efficiency for Tag", 			true,	"../EtaTag_Efficiency.png");
	createEfficiencyCanvas(pTagPhiEff,    "TagMuon_PhiEfficiency",	 "Angle Efficiency for Tag", 					true,	"../PhiTag_Efficiency.png");
}

//Call functions
void boss() {
	step1();
	//efficiencyOldMethod();
	efficiency();
}