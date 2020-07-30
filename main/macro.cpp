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
#include <TSystem.h>

#include "classes/Particle.h"

using namespace std;

//Select particles, draws and save histograms
void generateHistograms()
{
	//List of files
	const char *files[3] = {"data_histoall.root",
							"Run2011AMuOnia_mergeNtuple.root",
							"JPsiToMuMu_mergeMCNtuple.root"};



	//Options to change

	//Which file of files (variable above) should use
	int useFile = 0;

	//Path where is going to save results 
	const char* directoryToSave = "../result/";

	//Should limit data?
	Long64_t limitData = 0; //0 -> do not limit

	//Canvas drawing
	bool shouldDrawInvariantMassCanvas 			= true;
	bool shouldDrawInvariantMassCanvasRegion 	= false;
	bool shouldDrawQuantitiesCanvas 			= false;
	bool shouldDrawEfficiencyCanvas 			= false;


	
	//Check if dir exists and create
 	if (gSystem->AccessPathName(directoryToSave))
	{
		if (gSystem->mkdir(directoryToSave))
		{
			cout << "\"" << directoryToSave << "\" directory not found could not be created" << endl;
		}
		else
		{
			cout << "\"" << directoryToSave << "\" directory created OK" << endl;
		}
	}
	else
	{
		cout << "\"" << directoryToSave << "\" directory OK" << endl;
	}

	//Compatibility adjusts on file read (for data_histoall ntupples)
	string folderName = "tagandprobe/";
	if (useFile == 0)
		folderName = "demo/";

	//Open and read files
	TFile *file0  = TFile::Open(("../" + string(files[useFile])).data());
	TTree *TreePC = (TTree*)file0->Get((folderName + "PlotControl").data());
	TTree *TreeAT = (TTree*)file0->Get((folderName + "AnalysisTree").data());
	cout << "Using \"" << files[useFile] << "\" ntupple" << endl;
	
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

	//Get data size and set data limit if has
	Long64_t numberEntries = TreePC->GetEntries();
	if (limitData > 0 && limitData < numberEntries)
		numberEntries = limitData;
	printf("Data size analysed = %lld of %lld\n", numberEntries, TreePC->GetEntries());
	cout << endl;

	//Progress cout format
	string progressFormat = "%d/2 progress: %05.2f%% %0"+to_string(strlen(to_string(numberEntries).data()))+"d/%d\r";

	//Loop between the components
	for (int i = 0; i < numberEntries; i++)
	{
		//Show progress
		printf(progressFormat.data(), 1, (float)i/(float)numberEntries*100, i, numberEntries);

		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4)
		{
			if (PassingProbeTrackingMuon)
				Muon.TrackerPass.hMass->Fill(InvariantMass);
			Muon.All        .hMass->Fill(InvariantMass);
		}
	}

	Muon.doFit();

	Muon.updateSelectionParameters();


	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Create file root to store generated files
	TFile *generatedFile = TFile::Open((string(directoryToSave) + "generated_hist.root").data(),"RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->   cd("canvas/");

	if (shouldDrawInvariantMassCanvas)
	{
		bool drawRegions 	= false;
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.MassTracker.createCanvas(drawRegions, shouldWrite, shouldSave);
	}

	if (shouldDrawInvariantMassCanvasRegion)
	{
		bool drawRegions 	= true;
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.MassTracker.createCanvas(drawRegions, shouldWrite, shouldSave);
	}

	//Loop between the components again
	for (int i = 0; i < numberEntries; i++)
	{
		//Show progress
		printf(progressFormat.data(), 2, (float)i/(float)numberEntries*100, i, numberEntries);

		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Accepted particles
		if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4)
		{
			if (PassingProbeTrackingMuon)
				Muon.TrackerPass.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
			Muon.All        .fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
		}
	}
	cout << endl;

	Muon.subtractSigHistograms();

	if (shouldDrawQuantitiesCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.TrackerPass.createDividedCanvas(shouldWrite, shouldSave);
		Muon.All        .createDividedCanvas(shouldWrite, shouldSave);
	}

	//Debug consistency for histograms
	Muon.consistencyDebugCout();

	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	//Write histograms on file
	{
		//Write histograms on file
		bool writehSigBack 	= true;
		bool writehSig 		= true;
		bool writehBack 	= true;
	
		Muon.TrackerPass.write(writehSigBack, writehSig, writehBack);
		Muon.All        .write(writehSigBack, writehSig, writehBack);
	}

	//Write mass histograms on file
	{
		bool writehPass = true;
		bool writehAll 	= true;

		Muon.writeMassHistograms(writehPass, writehAll);
	}

	//Save plots
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	{
		bool shouldWrite 	= true;

		Muon.TrackerPass.createEfficiencyPlot(shouldWrite);
		Muon.All.createEfficiencyPlot(shouldWrite);
	}


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");

	if (shouldDrawEfficiencyCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSave 	= true;

		Muon.TrackerPass.createEfficiencyCanvas(shouldWrite, shouldSave);
		Muon.All.createEfficiencyCanvas(shouldWrite, shouldSave);
	}

	//Close files
	generatedFile->Close();
}

//Call functions
void macro()
{
	generateHistograms();
}