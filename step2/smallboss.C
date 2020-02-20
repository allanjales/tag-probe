 /*
!--------------------------------
!proposito: criar gráficos de eficiência dado histogramas
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
// Main functions
//-------------------------------------

//Estimates efficiency
void efficiency()
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
	hPtEff->Divide(hPtSig, hPtSigBack, 1.0, 1.0, "B");
	hPtEff->SetLineWidth(2);
	hPtEff->SetLineColor(kRed);
	hPtEff->SetMarkerStyle(21);
	hPtEff->SetMarkerSize(0.5);
	hPtEff->SetMarkerColor(kRed);
	TCanvas *c1 = new TCanvas("ProbePt_Efficiency","Probe Pt Efficiency", 800, 600);
	c1->SetTopMargin(0.07);
	c1->SetLeftMargin(0.12);
	c1->SetTicky(2);
	hPtEff->Draw();

	//Efficiency calculation for Eta
	TH1D* hEtaEff  = (TH1D*) hEtaSig->Clone("ProbeEta_Efficiency");
	hEtaEff->SetTitle("Pseudorapidity Efficiency for Probe");
	hEtaEff->Divide(hEtaSig, hEtaSigBack, 1.0, 1.0, "B");
	hEtaEff->SetLineWidth(2);
	hEtaEff->SetLineColor(kRed);
	hEtaEff->SetMarkerStyle(21);
	hEtaEff->SetMarkerSize(0.5);
	hEtaEff->SetMarkerColor(kRed);
	TCanvas *c2 = new TCanvas("ProbeEta_Efficiency","Probe Eta Efficiency", 800, 600);
	c2->SetTopMargin(0.07);
	c2->SetLeftMargin(0.12);
	c2->SetTicky(2);
	hEtaEff->Draw();

	//Efficiency calculation for Phi
	TH1D* hPhiEff  = (TH1D*) hPhiSig->Clone("ProbePhi_Efficiency");
	hPhiEff->SetTitle("Phi Efficiency for Probe");
	hPhiEff->Divide(hPhiSig, hPhiSigBack, 1.0, 1.0, "B");
	hPhiEff->SetLineWidth(2);
	hPhiEff->SetLineColor(kRed);
	hPhiEff->SetMarkerStyle(21);
	hPhiEff->SetMarkerSize(0.5);
	hPhiEff->SetMarkerColor(kRed);
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

//Calls functions
void smallboss() {
	efficiency();
}