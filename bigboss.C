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



//-------------------------------------
// Main functions
//-------------------------------------

//Draws and save invariant mass histogram for pp
void invariantMassProbe()
{
	TFile *file0 = TFile::Open("data_histoall.root");				//Opens the file
	TTree *Tree = (TTree*)file0->Get("demo/AnalysisTree");			//Opens TTree from file
	
	//Creates variables
	Double_t 	ProbeMuon_Pt;
	Double_t 	ProbeMuon_Eta;
	Double_t 	InvariantMass;

	int i;	//Itineration variable

	//Assign variables
	Tree->SetBranchAddress("ProbeMuon_Pt",		&ProbeMuon_Pt);
	Tree->SetBranchAddress("ProbeMuon_Eta",		&ProbeMuon_Eta);
	Tree->SetBranchAddress("InvariantMass",		&InvariantMass);

	//Create histogram
	TH1F *h1 = new TH1F("InvariantMassProbe","Invariant Mass (Probe);Mass (GeV/c^{2});Events / (0.1 GeV/c^{2})",240,2.8,3.4);	//Creates histogram

	//Reset title for Y axis wth the right bin width
	h1->GetYaxis()->SetTitle(Form("Events / (%1.4f GeV/c^{2})", h1->GetBinWidth(0)));

	//Loop between the components
	for (i = 0; i < Tree->GetEntries(); i++)
	{
		//Gets the entry from TTree
		Tree->GetEntry(i);

		//Fill the histogram then
		h1->Fill(InvariantMass);
	}



	//Create canvas
	TCanvas *c1 = new TCanvas("c1","Invariant Mass for pp", 600, 600);	//Creates Canvas for drawing

	//Set margin for canvas
	c1->SetTopMargin(0.07);
	c1->SetRightMargin(0.05);
	c1->SetBottomMargin(0.11);
	c1->SetLeftMargin(0.15);

	//Set title size for axis and other stuffs for histogram style
	h1->GetYaxis()->SetTitleSize(0.04);
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->CenterTitle(true);
	//h1->SetBit(TH1::kNoTitle);		//Not show histogram title
	h1->SetMarkerStyle(20);				//Set markers style
	h1->SetMarkerColor(kBlack);			//Set markers colors
	h1->SetLineColor(kBlack);			//Set lines colors (for errorbars)



	//Fit Function
	TF1 *f = new TF1("f",FFit_InvariantMassProbe,h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax(),12);
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
	TFitResultPtr fitr = h1->Fit(f,"RNS","",h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());


	

	//Get parameters from fit function and put it in par variable
   	Double_t par[12];
	f->GetParameters(par);

	//Signal Fitting
	TF1 *fs = new TF1("fs",Signal_InvariantMassProbe,h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax(),8);
	fs->SetNpx(1000);				//Resolution of background fit function
	fs->SetParameters(par);			//Get only background part
	fs->SetLineColor(kMagenta); 	//Fit Color
	fs->SetLineStyle(kSolid);		//Fit Style

	//Background Fitting
	TF1 *fb = new TF1("fb",Background_InvariantMassProbe,h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax(),4);
	fb->SetNpx(1000);				//Resolution of background fit function
	fb->SetParameters(&par[8]);		//Get only background part
	fb->SetLineColor(kBlue); 		//Fit Color
	fb->SetLineStyle(kDashed);		//Fit Style

	//Draws histogram & fit function
	h1->Draw("ep");
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
	tx->DrawLatex(0.6,0.60,Form("#chi^{2}/ndf = %g/%d",fitr->Chi2(),fitr->Ndf()));

	//Show number of particles
	tx->DrawLatex(0.6,0.48,Form("%3.0f J/#psi",fs->Integral(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax()) * (1/h1->GetBinWidth(0))));


	//Add legend
	TLegend *l = new TLegend(0.65,0.75,0.92,0.90);
	l->SetTextSize(0.04);
	l->AddEntry(h1,	"J/#psi","lp");
	l->AddEntry(f,	"Fitting","l");
	//l->AddEntry(fs,	"Signal","l");
	l->AddEntry(fb,	"Background","l");
	l->Draw();


	//Not show frame with mean, std dev
	gStyle->SetOptStat(0);

	//Saves as image
	c1->SaveAs("InvariantMassProbe.png");



	//SHow chi-squared test (on prompt)
	cout << endl;
	cout << "Chi2/ndf = " << f->GetChisquare()/f->GetNDF() << endl;
	cout << endl;
}

//Calls functions to draw and save invariant mass histogram
void bigboss() {
	invariantMassProbe();	//Draws and save invariant mass histogram for pp collision
}
