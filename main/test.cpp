 /*
!--------------------------------
!Purpose: Test root files created by this program
!--------------------------------	
!author: Allan Jales
!--------------------------------
*/

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <string>

#include "FitFunctions.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

using namespace std;
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

///Objeto que vai dizer o chi global
struct GlobalChi2 {

   ///Bagulho de duas funcoes
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)

   ///funcao que soma os parametros
   double operator() (const double *par) const {
      double p1[12];
      for (int i = 0; i < 12; ++i) p1[i] = par[iparPass[i]];

      double p2[24];
      for (int i = 0; i < 24; ++i) p2[i] = par[iparPassFail[i]];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }

	///Construtor do objeto! Define o que sao as funoes
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}
      ///Saporra de cima é so pra iniciar os objetos antes de fazer o construtor
      ///Vai la e associa as funcoes com esses objetos
};

void openFile()
{
	const char *files[3] = {"raphael.root",
							"run2011.root",
							"mc2020.root"};

	int useFile = 1;

	//Open and read files
	TFile *file0 = TFile::Open((string("resultsForTesting/") + string(files[useFile])).data());

	//Create histograms
	TH1D* hPass 	= new TH1D("PassingMuonInvariantMass", "InvariantMass (Passing);Mass (GeV/c^{2}); Events", 240, 2.8, 3.4);
	TH1D* hFail 	= new TH1D("FailingMuonInvariantMass", "InvariantMass (Failing);Mass (GeV/c^{2}); Events", 240, 2.8, 3.4);
	TH1D* hPassFail = new TH1D("BothMuonInvariantMass", "InvariantMass (Both);Mass (GeV/c^{2}); Events", 240, 2.8, 3.4);

	//Seta parametros
	//Passing Fitting
	TF1* fPass = new TF1("FitFunction_Pass", FitFunctions::Merged::Pass_InvariantMass, 2.8, 3.4, 12);
	fPass->SetNpx(1000);
	fPass->SetLineColor(kMagenta);
	fPass->SetLineStyle(kSolid);
	fPass->SetParameter(0,		340.2);
	fPass->SetParameter(1,		3.09);
	fPass->SetParameter(2,		0.037);
	fPass->SetParameter(3,		1.824);
	fPass->SetParameter(4,		1.034);
	fPass->SetParameter(5,		3.093);
	fPass->SetParameter(6,		0.022);
	fPass->SetParameter(7,		8322.27);
	fPass->SetParameter(8,		-0.217);
	fPass->SetParameter(9,		1.915);
	fPass->SetParameter(10, 	263.185);
	fPass->SetParameter(11,		0.061);
	
	//Rename parameters
	int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);
	for (int i = 0; i < arraySize; i++)
	{
		fPass->SetParName(i, fittingParName[i]);
	}

	//Failing Fitting
	TF1* fFail = new TF1("FitFunction_Fail", FitFunctions::Merged::Fail_InvariantMass, 2.8, 3.4, 12);
	fFail->SetNpx(1000);
	fFail->SetLineColor(kBlue);
	fFail->SetLineStyle(kDashed);
	fFail->SetParameter(0,		340.2);
	fFail->SetParameter(1,		3.09);
	fFail->SetParameter(2,		0.037);
	fFail->SetParameter(3,		1.824);
	fFail->SetParameter(4,		1.034);
	fFail->SetParameter(5,		3.093);
	fFail->SetParameter(6,		0.022);
	fFail->SetParameter(7,		8322.27);
	fFail->SetParameter(8,		-0.217);
	fFail->SetParameter(9,		1.915);
	fFail->SetParameter(10, 	263.185);
	fFail->SetParameter(11,		0.061);
	
	//Rename parameters
	for (int i = 0; i < arraySize; i++)
	{
		fFail->SetParName(i, fittingParName[i]);
	}

	hPass = (TH1D*)file0->Get("histograms/PassingMuonInvariantMass");
	hFail = (TH1D*)file0->Get("histograms/FailingMuonInvariantMass");
	hPassFail = (TH1D*)file0->Get("histograms/AllMuonInvariantMass");

	//Both Fitting
	TF1* fPassFail = new TF1("FitFunction_Both", FitFunctions::Merged::Both_InvariantMass, 2.8, 3.4, 24);
	fPassFail->SetNpx(1000);
	fPassFail->SetLineColor(kBlue);
	fPassFail->SetLineStyle(kSolid);
	
	//Rename parameters
	for (int i = 0; i < arraySize*2; i++)
	{
		fPassFail->SetParName(i, fittingParName[i%arraySize]);
	}

	//Importantes não entendidos

	ROOT::Math::WrappedMultiTF1 wfPass(*fPass, 1);
	ROOT::Math::WrappedMultiTF1 wfPassFail(*fPassFail, 1);

	ROOT::Fit::DataOptions opt;

	ROOT::Fit::DataRange rangePass;
	rangePass.SetRange(2.8, 3.4);
	ROOT::Fit::BinData dataPass(opt, rangePass);
	ROOT::Fit::FillData(dataPass, hPass);

	ROOT::Fit::DataRange rangePassFail;
	rangePassFail.SetRange(2.8, 3.4);
	ROOT::Fit::BinData dataPassFail(opt, rangePassFail);
	ROOT::Fit::FillData(dataPassFail, hPassFail);

	ROOT::Fit::Chi2Function chi2_Pass(dataPass, wfPass);
	ROOT::Fit::Chi2Function chi2_PassFail(dataPassFail, wfPassFail);

	GlobalChi2 globalChi2(chi2_Pass, chi2_PassFail);

	ROOT::Fit::Fitter fitter;

	//?
	double par0[24] = {340.2,
						3.09,
						0.037,
						1.824,	//Daqui pra baixo começa a repetir
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


   TCanvas * c1 = new TCanvas("Simfit","Simultaneous fit of two histograms", 10,10,700,700);
   c1->Divide(1,2);
   c1->cd(1);
   gStyle->SetOptFit(1111);

   fPass->SetFitResult(result, iparPass);
   fPass->SetRange(rangePass().first, rangePass().second);
   fPass->SetLineColor(kBlue);
   hPass->GetListOfFunctions()->Add(fPass);
   hPass->Draw();

   c1->cd(2);
   fPassFail->SetFitResult(result, iparPassFail);
   fPassFail->SetRange(rangePassFail().first, rangePassFail().second);
   fPassFail->SetLineColor(kRed);
   hPassFail->GetListOfFunctions()->Add(fPassFail);
   hPassFail->Draw();
}

//Call functions
void test()
{
	openFile();
}