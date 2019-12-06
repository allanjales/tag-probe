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
	TTree *Tree = (TTree*)file0->Get("demo/PlotControl");	//Opens TTree from file
	
	//Creates variables
	Double_t 	ProbeMuon_Pt;
	Double_t 	ProbeMuon_Eta;
	Double_t 	ProbeMuon_Phi;
	Double_t 	InvariantMass;

	int i;	//Itineration variable

	//Variables for side band subtraction
	int count_sigregion = 0;	//Inside sigma
	int count_sideband = 0;		//Outside sigma

	//Jpsi constants (PDG values)
	const double M_JPSI = 3.097;
	const double W_JPSI = 0.010;

	//Assign variables
	Tree->SetBranchAddress("ProbeMuon_Pt",		&ProbeMuon_Pt);
	Tree->SetBranchAddress("ProbeMuon_Eta",		&ProbeMuon_Eta);
	Tree->SetBranchAddress("ProbeMuon_Phi",		&ProbeMuon_Phi);
	Tree->SetBranchAddress("InvariantMass",		&InvariantMass);

	//Create histograms for Singla + Background
	TH1D *hMassAll = new TH1D("InvariantMassProbe","All Invariant Mass (Probe);Mass (GeV/c^{2});Events",240,2.8,3.4);
	TH1D *hPtSigBack = new TH1D("TransversalMomentumProbeAll","Transversal Momentum (Probe);Transversal Momentum (GeV/c);Events",100,0.,100.);
	TH1D *hEtaSigBack = new TH1D("PseudorapidityProbeAll","Pseudorapidity (Probe);#eta;Events",200,-2.5,2.5);
	TH1D *hPhiSigBack = new TH1D("PhiProbeAll","Phi (Probe);#phi;Events",100,-3.15,3.15);

	//Clone histogram for Background
	TH1D* hPtBack  = (TH1D*) hPtSigBack->Clone("hPtBack");
	TH1D* hEtaBack = (TH1D*) hEtaSigBack->Clone("hEtaBack");
	TH1D* hPhiBack = (TH1D*) hPhiSigBack->Clone("hPhiBack");

	//Rename histogram
	hPtBack->SetTitle("Transversal Momentum (Probe)");
	hPtBack->SetTitle("Pseudorapidity (Probe)");
	hPtBack->SetTitle("Phi (Probe)");

	//Reset title for Y axis with the right bin width
	hMassAll->GetYaxis()->SetTitle(Form("Events / (%1.4f GeV/c^{2})", 	hMassAll->GetBinWidth(0)));
	hPtSigBack->GetYaxis()->SetTitle(Form("Events / (%1.1f GeV/c)", 	hPtSigBack->GetBinWidth(0)));
	hEtaSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", 			hEtaSigBack->GetBinWidth(0)));
	hPhiSigBack->GetYaxis()->SetTitle(Form("Events / (%1.3f)", 			hPhiSigBack->GetBinWidth(0)));


	//Loop between the components
	for (i = 0; i < Tree->GetEntries(); i++)
	{
		//Gets the entry from TTree
		Tree->GetEntry(i);

		//Fill the histogram then
		hMassAll->Fill(InvariantMass);

		//if is inside signal region
		if (fabs(InvariantMass - M_JPSI) < W_JPSI*3.0)
		{
			hPtSigBack->Fill(ProbeMuon_Pt);		//Add to Pt histogram
			hEtaSigBack->Fill(ProbeMuon_Eta);	//Add to Eta histogram
			hPhiSigBack->Fill(ProbeMuon_Phi);	//Add to Phi histogram
			count_sigregion++;
		}

		//If is inside sideband region
		if (fabs(InvariantMass - M_JPSI) > W_JPSI*3.5 && fabs(InvariantMass - M_JPSI) < W_JPSI*6.5)
		{
			hPtBack->Fill(ProbeMuon_Pt);	//Add to Pt histogram
			hEtaBack->Fill(ProbeMuon_Eta);	//Add to Eta hitogram
			hPhiBack->Fill(ProbeMuon_Phi);	//Add to Phi hitogram
			count_sideband++;
		}
	}

	///Clone histogram for Signal
	TH1D* hPtSig  = (TH1D*) hPtSigBack->Clone("hPtSig");
	TH1D* hEtaSig = (TH1D*) hEtaSigBack->Clone("hEtaSig");
	TH1D* hPhiSig = (TH1D*) hPhiSigBack->Clone("hPhiSig");

	//Subtract histograms
	hPtSig->Add(hPtBack,-1);
	hEtaSig->Add(hEtaBack,-1);
	hPhiSig->Add(hPhiBack,-1);

	//Number of signal and uncertain
	double S = count_sigregion - count_sideband;
	double dS = sqrt(count_sigregion+count_sideband);

	//-------------------------------------
	// Canvas Invariant Mass
	//-------------------------------------

	//Create canvas
	TCanvas *c1 = new TCanvas("c1","Invariant Mass", 600, 600);			//Creates Canvas for drawing

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



	//-------------------------------------
	// Canvas Transversal Momentum
	//-------------------------------------

	//Create canvas
	TCanvas *c2 = new TCanvas("c2","Transversal Momentum", 1200, 600);

	//Divide canvas
	c2->Divide(2,1);

	//Select canvas part
	c2->cd(1);

	//Set margin for canvas part
	c2->cd(1)->SetTopMargin(0.07);
	c2->cd(1)->SetRightMargin(0.02);
	c2->cd(1)->SetBottomMargin(0.11);
	c2->cd(1)->SetLeftMargin(0.13);

	//Draws Main histogram
	hPtSigBack->Draw();

	//Draws Background histogram
	hPtBack->SetLineColor(kBlue); 		//Fit Color
	hPtBack->SetLineStyle(kDashed);		//Fit Style
	hPtBack->Draw("same");

	//Draws Signal histogram
	hPtSig->SetLineColor(kMagenta); 	//Fit Color
	hPtSig->Draw("same");

	//Add legend
	TLegend *l2_1 = new TLegend(0.65,0.75,0.92,0.90);
	l2_1->SetTextSize(0.04);
	l2_1->AddEntry(hEtaSigBack,	"All","lp");
	l2_1->AddEntry(hEtaSig,		"Signal","l");
	l2_1->AddEntry(hEtaBack,	"Background","l");
	l2_1->Draw();
	
	//Draws information
	TLatex *tx2_1 = new TLatex();
	tx2_1->SetTextSize(0.04);
	tx2_1->SetTextAlign(12);	//Align left, center
	tx2_1->SetTextFont(42);
	tx2_1->SetNDC(kTRUE);
	tx2_1->DrawLatex(0.48,0.50,Form("%g entries (total)",		hPtSigBack->GetEntries()));
	tx2_1->DrawLatex(0.48,0.45,Form("%g entries (signal)",		hPtSig->GetEntries()));
	tx2_1->DrawLatex(0.48,0.40,Form("%g entries (background)",	hPtBack->GetEntries()));
	//tx2_1->DrawLatex(0.48,0.35,Form("%g entries (inverse)",		hPtInverse->GetEntries()));

	//Select canvas part
	c2->cd(2);

	//Set margin for canvas part
	c2->cd(2)->SetTopMargin(0.07);
	c2->cd(2)->SetRightMargin(0.02);
	c2->cd(2)->SetBottomMargin(0.11);
	c2->cd(2)->SetLeftMargin(0.13);

	//Same range as comparision
	//hPtSig->GetYaxis()->SetRangeUser(0,hPtSigBack->GetYaxis()->GetYmax());
	hPtSig->SetMinimum(0);
	hPtSig->SetTitle("Transversal Momentum of Signal (Probe)");

	//Signal
	hPtSig->SetLineColor(kMagenta); 		//Fit Color
	hPtSig->Draw("same");

	//Add legend
	TLegend *l2_2 = new TLegend(0.65,0.85,0.92,0.90);
	l2_2->SetTextSize(0.04);
	l2_2->AddEntry(hEtaSig,		"Signal","l");
	l2_2->Draw();

	//Draws information
	TLatex *tx2_2 = new TLatex();
	tx2_2->SetTextSize(0.04);
	tx2_2->SetTextAlign(12);	//Align left, center
	tx2_2->SetTextFont(42);
	tx2_2->SetNDC(kTRUE);
	tx2_2->DrawLatex(0.48,0.50,Form("%g entries (signal)", hPtSig->GetEntries()));



	//-------------------------------------
	// Canvas Pseudorapidity
	//-------------------------------------

	//Create canvas
	TCanvas *c3 = new TCanvas("c3","Pseudorapidity", 1200, 600);

	//Same range as comparision
	//hPtSig->GetYaxis()->SetRangeUser(0,hPtSigBack->GetYaxis()->GetYmax());
	hEtaSig->SetMinimum(0);
	hEtaSig->SetTitle("Pseudorapidity of Signal (Probe)");

	//Divide canvas
	c3->Divide(2,1);

	//Select canvas part
	c3->cd(1);

	//Set margin for canvas part
	c3->cd(1)->SetTopMargin(0.07);
	c3->cd(1)->SetRightMargin(0.02);
	c3->cd(1)->SetBottomMargin(0.11);
	c3->cd(1)->SetLeftMargin(0.13);
	
	//Draws Main histogram
	hEtaSigBack->Draw();

	//Draws Signal histogram
	hEtaSig->SetLineColor(kMagenta); 	//Fit Color
	hEtaSig->Draw("same");

	//Draws Background histogram
	hEtaBack->SetLineColor(kBlue); 		//Fit Color
	hEtaBack->SetLineStyle(kDashed);	//Fit Style
	hEtaBack->Draw("same");

	//Add legend
	TLegend *l3_1 = new TLegend(0.65,0.75,0.92,0.90);
	l3_1->SetTextSize(0.04);
	l3_1->AddEntry(hEtaSigBack,	"All","lp");
	l3_1->AddEntry(hEtaSig,		"Signal","l");
	l3_1->AddEntry(hEtaBack,	"Background","l");
	l3_1->Draw();
	
	//Draws information
	TLatex *tx3_1 = new TLatex();
	tx3_1->SetTextSize(0.04);
	tx3_1->SetTextAlign(22);	//Align center, center
	tx3_1->SetTextFont(42);
	tx3_1->SetNDC(kTRUE);
	tx3_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hEtaSigBack->GetEntries()));
	tx3_1->DrawLatex(0.55,0.45,Form("%g entries (signal)",		hEtaSig->GetEntries()));
	tx3_1->DrawLatex(0.55,0.40,Form("%g entries (background)",	hEtaBack->GetEntries()));

	//Select canvas part
	c3->cd(2);

	//Set margin for canvas part
	c3->cd(2)->SetTopMargin(0.07);
	c3->cd(2)->SetRightMargin(0.02);
	c3->cd(2)->SetBottomMargin(0.11); 
	c3->cd(2)->SetLeftMargin(0.13);

	//Signal
	hEtaSig->SetLineColor(kMagenta);	//Fit Color
	hEtaSig->Draw("same");

	/*
	//Set width of line
	hEtaSigBack->SetLineWidth(2);
	hEtaSig->SetLineWidth(2);
	hEtaBack->SetLineWidth(2);
	*/
	
	/*
	//Set fill color
	hEtaSig->SetFillColor(kMagenta);
	hEtaSigBack->SetFillColor(kBlue);
	hEtaBack->SetFillColor(kYellow);
	*/

	//Add legend
	TLegend *l3_2 = new TLegend(0.65,0.85,0.92,0.90);
	l3_2->SetTextSize(0.04);
	l3_2->AddEntry(hEtaSig,		"Signal","l");
	l3_2->Draw();

	//Draws information
	TLatex *tx3_2 = new TLatex();
	tx3_2->SetTextSize(0.04);
	tx3_2->SetTextAlign(22);	//Align center, center
	tx3_2->SetTextFont(42);
	tx3_2->SetNDC(kTRUE);
	tx3_2->DrawLatex(0.55,0.5,Form("%g entries (signal)",hEtaSig->GetEntries()));



	//-------------------------------------
	// Canvas Pseudorapidity
	//-------------------------------------

	//Create canvas
	TCanvas *c4 = new TCanvas("c4","Phi", 1200, 600);

	//Same range as comparision
	//hPtSig->GetYaxis()->SetRangeUser(0,hPtSigBack->GetYaxis()->GetYmax());
	hPhiSig->SetMinimum(0);
	hPhiSig->SetTitle("Phi of Signal (Probe)");

	//Divide canvas
	c4->Divide(2,1);

	//Select canvas part
	c4->cd(1);

	//Set margin for canvas part
	c4->cd(1)->SetTopMargin(0.07);
	c4->cd(1)->SetRightMargin(0.02);
	c4->cd(1)->SetBottomMargin(0.11);
	c4->cd(1)->SetLeftMargin(0.13);
	
	//Draws Main histogram
	hPhiSigBack->SetMinimum(0);
	hPhiSigBack->Draw();

	//Draws Signal histogram
	hPhiSig->SetLineColor(kMagenta); 	//Fit Color
	hPhiSig->Draw("same");

	//Draws Background histogram
	hPhiBack->SetLineColor(kBlue); 		//Fit Color
	hPhiBack->SetLineStyle(kDashed);	//Fit Style
	hPhiBack->Draw("same");

	//Add legend
	TLegend *l4_1 = new TLegend(0.65,0.75,0.92,0.90);
	l4_1->SetTextSize(0.04);
	l4_1->AddEntry(hPhiSigBack,	"All","lp");
	l4_1->AddEntry(hPhiSig,		"Signal","l");
	l4_1->AddEntry(hPhiBack,	"Background","l");
	l4_1->Draw();
	
	//Draws information
	TLatex *tx4_1 = new TLatex();
	tx4_1->SetTextSize(0.04);
	tx4_1->SetTextAlign(22);	//Align center, center
	tx4_1->SetTextFont(42);
	tx4_1->SetNDC(kTRUE);
	tx4_1->DrawLatex(0.55,0.50,Form("%g entries (total)",		hPhiSigBack->GetEntries()));
	tx4_1->DrawLatex(0.55,0.45,Form("%g entries (signal)",		hPhiSig->GetEntries()));
	tx4_1->DrawLatex(0.55,0.40,Form("%g entries (background)",	hPhiBack->GetEntries()));

	//Select canvas part
	c4->cd(2);

	//Set margin for canvas part
	c4->cd(2)->SetTopMargin(0.07);
	c4->cd(2)->SetRightMargin(0.02);
	c4->cd(2)->SetBottomMargin(0.11); 
	c4->cd(2)->SetLeftMargin(0.13);

	//Signal
	hPhiSig->SetLineColor(kMagenta);	//Fit Color
	hPhiSig->Draw("same");

	/*
	//Set width of line
	hPhiSigBack->SetLineWidth(2);
	hPhiSig->SetLineWidth(2);
	hPhiBack->SetLineWidth(2);
	*/

	/*
	//Set fill color
	hPhiSig->SetFillColor(kMagenta);
	hPhiSigBack->SetFillColor(kBlue);
	hPhiBack->SetFillColor(kYellow);
	*/

	//Add legend
	TLegend *l4_2 = new TLegend(0.65,0.85,0.92,0.90);
	l4_2->SetTextSize(0.04);
	l4_2->AddEntry(hPhiSig,		"Signal","l");
	l4_2->Draw();

	//Draws information
	TLatex *tx4_2 = new TLatex();
	tx4_2->SetTextSize(0.04);
	tx4_2->SetTextAlign(22);	//Align center, center
	tx4_2->SetTextFont(42);
	tx4_2->SetNDC(kTRUE);
	tx4_2->DrawLatex(0.55,0.5,Form("%g entries (signal)",hPhiSig->GetEntries()));



	//Saves as image
	c1->SaveAs("InvariantMassProbe.png");
	c2->SaveAs("PtProbe.png");
	c3->SaveAs("EtaProbe.png");
	c4->SaveAs("PhiProbe.png");



	//Integrate function to get number of particles in it
	cout << endl;
	cout << "Candidates by integrating" << endl;
	cout << "#Signal     = " << fs->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;
	cout << "#Background = " << fb->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;
	cout << "#Total      = " << f ->Integral(hMassAll->GetXaxis()->GetXmin(), hMassAll->GetXaxis()->GetXmax()) * hMassAll->GetNbinsX() << endl;

	//Debug
	cout << endl;
	cout << "Candidates by sideband" << endl;
	cout << "#Tree Entries     = " << Tree->GetEntries() << endl;
	cout << "#Total            = "	<< count_sigregion 					<< endl;
	cout << "#Sideband         = "	<< count_sideband 					<< endl;
	cout << "#Total - Sideband = " 	<< count_sigregion - count_sideband << endl;


	//Show chi-squared test (on prompt)
	cout << endl;
	cout << "Overview"<< endl;
	cout << "Chi2/ndf = " << f->GetChisquare()/f->GetNDF() << endl;
	printf( "#Signal  = %.0f +- %.0f\n",S,dS);
	cout << endl;
}

//Calls functions to draw and save invariant mass histogram
void bigboss() {
	invariantMassProbe();	//Draws and save invariant mass histogram for pp collision
}