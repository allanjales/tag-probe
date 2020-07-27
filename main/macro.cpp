 /*
!--------------------------------
!Purpose: Find efficiency using Tag And Probe method
!--------------------------------	
!author: Allan Jales
!--------------------------------
*/

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

#include "classes/Particle.h"

using namespace std;

//Select particles, draws and save histograms
void generateHistograms(bool shouldDrawInvariantMassCanvas = true, bool shouldDrawInvariantMassCanvasRegion = true, bool shouldDrawQuantitiesCanvas = true, bool shouldDrawEfficiencyCanvas = true)
{
	const char *files[3] = {"data_histoall.root",
							"Run2011AMuOnia_mergeNtuple.root",
							"JPsiToMuMu_mergeMCNtuple.root"};

	int useFile = 0;

	//Compatibility adjusts on file read
	string folderName = "tagandprobe/";
	if (useFile == 0)
		folderName = "demo/";

	//Open and read files
	TFile *file0  = TFile::Open(("../" + string(files[useFile])).data());
	TTree *TreePC = (TTree*)file0->Get((folderName + "PlotControl").data());
	TTree *TreeAT = (TTree*)file0->Get((folderName + "AnalysisTree").data());

	int limitData = 0; //0 -> do not limit
	
	//Create variables
	double ProbeMuon_Pt;
	double ProbeMuon_Eta;
	double ProbeMuon_Phi;
	double TagMuon_Pt;
	double TagMuon_Eta;
	double TagMuon_Phi;
	double InvariantMass;
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
	if (useFile == 0)
	TreePC->SetBranchAddress("InvariantMass",				&InvariantMass);
	else
	TreeAT->SetBranchAddress("InvariantMass",				&InvariantMass);
	TreeAT->SetBranchAddress("PassingProbeTrackingMuon",	&PassingProbeTrackingMuon);
	TreeAT->SetBranchAddress("PassingProbeStandAloneMuon",	&PassingProbeStandAloneMuon);
	TreeAT->SetBranchAddress("PassingProbeGlobalMuon",		&PassingProbeGlobalMuon);

	//Create a object
	Particle Muon;
	Muon.setMethod(1);

	int numberEntries = TreePC->GetEntries();

	if (limitData > 0)
		numberEntries = limitData;

	string progressFormat = "%d/2 progress: %05.2f%% %0"+to_string(strlen(to_string(numberEntries).data()))+"d/%d\r";

	//Loop between the components
	for (int i = 0; i < numberEntries; i++)
	{
		//Show progress
		printf((progressFormat).data(), 1, (float)i/(float)numberEntries*100, i, numberEntries);

		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4)
		{
			if (PassingProbeTrackingMuon)
				Muon.Pass.Mass.fill(InvariantMass);
			else
				Muon.Fail.Mass.fill(InvariantMass);
			Muon.Both.Mass.fill(InvariantMass);
		}
	}

	Muon.doFit();
	Muon.updateSelectionParameters();

	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Create file root to store generated files
	TFile *generatedFile = TFile::Open("../result/generated_hist.root","RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->   cd("canvas/");

	if (shouldDrawInvariantMassCanvas)
	{
		bool drawRegions 	= false;
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.Pass.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
		Muon.Fail.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
		Muon.Both.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
	}

	if (shouldDrawInvariantMassCanvasRegion)
	{
		bool drawRegions 	= true;
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.Pass.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
		Muon.Fail.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
		Muon.Both.Mass.createCanvas(drawRegions, shouldWrite, shouldSave);
	}

	//Loop between the components again
	for (int i = 0; i < numberEntries; i++)
	{
		//Show progress
		printf((progressFormat).data(), 2, (float)i/(float)numberEntries*100, i, numberEntries);

		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4)
		{
			if (PassingProbeTrackingMuon)
				Muon.Pass.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			else
				Muon.Fail.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			Muon.Both.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}
	cout << endl;

	//Debug
	Muon.massDebugCout();
	
	Muon.subtractSigHistograms();

	if (shouldDrawQuantitiesCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.Pass.createDividedCanvas(shouldWrite, shouldSave);
		Muon.Fail.createDividedCanvas(shouldWrite, shouldSave);
		Muon.Both.createDividedCanvas(shouldWrite, shouldSave);
	}

	//Debug
	Muon.consistencyDebugCout();

	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	//Write histograms on file
	bool writehSigBack 	= true;
	bool writehSig 		= true;
	bool writehBack 	= true;

	Muon.Pass.write(writehSigBack, writehSig, writehBack);
	Muon.Fail.write(writehSigBack, writehSig, writehBack);
	Muon.Both.write(writehSigBack, writehSig, writehBack);

	//Write mass histograms
	Muon.writehMass();

	//Save plots
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	{
		bool shouldWrite 	= true;

		Muon.Pass.createEfficiencyPlot(shouldWrite);
		Muon.Fail.createEfficiencyPlot(shouldWrite);
		Muon.Both.createEfficiencyPlot(shouldWrite);
	}


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");

	if (shouldDrawEfficiencyCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.Pass.createEfficiencyCanvas(shouldWrite, shouldSave);
		Muon.Fail.createEfficiencyCanvas(shouldWrite, shouldSave);
		Muon.Both.createEfficiencyCanvas(shouldWrite, shouldSave);
	}

	//Close files
	generatedFile->Close();
}

//Call functions
void macro()
{
	bool shouldDrawInvariantMassCanvas 			= false;
	bool shouldDrawInvariantMassCanvasRegion 	= false;
	bool shouldDrawQuantitiesCanvas 			= false;
	bool shouldDrawEfficiencyCanvas 			= false;

	generateHistograms(shouldDrawInvariantMassCanvas, shouldDrawInvariantMassCanvasRegion, shouldDrawQuantitiesCanvas, shouldDrawEfficiencyCanvas);
}