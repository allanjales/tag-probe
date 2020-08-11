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
#include <chrono>
#include <TSystem.h>

#include "classes/Particle.h"
#include "config/cutsAndFill.h"

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

	//Choose method
	//if 1 -> sideband by histogram || if 2 -> sideband by fitting
	int method = 1;

	//Path where is going to save results 
	const char* directoryToSave = "../result/";

	//Should limit data?
	long long limitData = 0; //0 -> do not limit

	//Canvas drawing
	bool shouldDrawInvariantMassCanvas 			= true;
	bool shouldDrawInvariantMassCanvasRegion 	= false;
	bool shouldDrawQuantitiesCanvas 			= false;
	bool shouldDrawEfficiencyCanvas 			= false;





	//Check if the name of dir is ok
	if (string(directoryToSave).back() != string("/"))
	{
		cerr << "To avoid errors, please, end the result directory with a \"/\"" << endl;
		abort();
	}
	
	//Check if dir exists and create
	if (gSystem->AccessPathName(directoryToSave))
	{
		if (gSystem->mkdir(directoryToSave))
		{
			cerr << "\"" << directoryToSave << "\" directory not found and could not be created ERROR" << endl;
			abort();
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

	double* quantities[7] = {&ProbeMuon_Pt,
							&ProbeMuon_Eta,
							&ProbeMuon_Phi,
							&TagMuon_Pt,
							&TagMuon_Eta,
							&TagMuon_Phi,
							&InvariantMass,
		};
	
	int* types[3] = {&PassingProbeTrackingMuon,
					&PassingProbeStandAloneMuon,
					&PassingProbeGlobalMuon
		};

	//Create a object and set method
	Particle Muon;
	Muon.setMethod(method);

	//Get data size and set data limit if has
	long long numberEntries = TreePC->GetEntries();
	if (limitData > 0 && limitData < numberEntries)
		numberEntries = limitData;
	printf("Data size analysed = %lld of %lld\n", numberEntries, TreePC->GetEntries());
	cout << endl;



	//Prepare for showing progress
	string progressFormat = "%d/%d progress: %05.2f%% %0"+to_string(strlen(to_string(numberEntries).data()))+"lld/%lld\r";
	auto lastTime = std::chrono::high_resolution_clock::now();

	//Loop between the components
	for (long long i = 0; i < numberEntries; i++)
	{
		//Select particle pair
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Show progress on screen
		if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - lastTime).count() >= 1000 || i == numberEntries - 1)
		{
			printf(progressFormat.data(), 1, 2, (float)(i+1)/(float)numberEntries*100, i+1, numberEntries);
			lastTime = std::chrono::high_resolution_clock::now();
		}

		//Fill histograms
		if (applyCuts(quantities, types))
		{
			fillMassHistograms(&Muon, quantities, types);
		}
	}
	cout << endl;



	Muon.doFit();

	Muon.updateSelectionParameters();


	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Supress canvas
	//gROOT->SetBatch(1);

	//Create file root to store generated files
	TFile* generatedFile = TFile::Open((string(directoryToSave) + "generated_hist.root").data(),"RECREATE");
	generatedFile->mkdir("canvas/");
	generatedFile->   cd("canvas/");

	if (shouldDrawInvariantMassCanvas)
	{
		bool drawRegions 	= false;
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		Muon.MassTracker   .createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
		Muon.MassStandalone.createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
		Muon.MassGlobal    .createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
	}

	if (shouldDrawInvariantMassCanvasRegion)
	{
		bool drawRegions 	= true;
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		Muon.MassTracker   .createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
		Muon.MassStandalone.createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
		Muon.MassGlobal    .createCanvas(drawRegions, shouldWrite, directoryToSave, shouldSavePNG);
	}



	//Loop between the components again
	for (long long i = 0; i < numberEntries; i++)
	{
		//Select particle pair
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Show progress on screen
		if (std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - lastTime).count() >= 1000 || i == numberEntries - 1)
		{
			printf(progressFormat.data(), 1, 2, (float)(i+1)/(float)numberEntries*100, i+1, numberEntries);
			lastTime = std::chrono::high_resolution_clock::now();
		}

		//Fill histograms
		if (applyCuts(quantities, types))
		{	
			fillQuantitiesHistograms(&Muon, quantities, types);
		}
	}
	cout << endl;



	Muon.subtractSigHistograms();

	if (shouldDrawQuantitiesCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		Muon.Tracker   .createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.Standalone.createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.Global    .createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.All       .createDividedCanvas(shouldWrite, directoryToSave, shouldSavePNG);
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
	
		Muon.Tracker   .write(writehSigBack, writehSig, writehBack);
		Muon.Standalone.write(writehSigBack, writehSig, writehBack);
		Muon.Global    .write(writehSigBack, writehSig, writehBack);
		Muon.All       .write(writehSigBack, writehSig, writehBack);
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

		Muon.Tracker   .createEfficiencyPlot(shouldWrite);
		Muon.Standalone.createEfficiencyPlot(shouldWrite);
		Muon.Global    .createEfficiencyPlot(shouldWrite);
		Muon.All       .createEfficiencyPlot(shouldWrite);
	}


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");

	if (shouldDrawEfficiencyCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		Muon.Tracker   .createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.Standalone.createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.Global    .createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
		Muon.All       .createEfficiencyCanvas(shouldWrite, directoryToSave, shouldSavePNG);
	}

	//Close files
	generatedFile->Close();
}

//Call functions
void macro()
{
	generateHistograms();
}