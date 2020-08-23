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

#include "classes/TagAndProbe.h"
#include "config/cuts.h"

using namespace std;

//Select particles, draws and save histograms
void generateHistograms()
{
	//List of files
	const char *files[] = {"../data_histoall.root",
							"../Run2011AMuOnia_mergeNtuple.root",
							"../JPsiToMuMu_mergeMCNtuple.root",
							"../data_histoall.root",
							"../Upsilon1SToMuMu_MC_full.root"};



	//Options to change

	//Which file of files (variable above) should use
	int useFile = 2;

	//Choose method
	//if 1 -> sideband by histogram || if 2 -> sideband by fitting
	int method = 2;

	//Path where is going to save results 
	const char* directoryToSave = "../result/";

	//Should limit data?
	long long limitData = 100000; //0 -> do not limit

	//Canvas drawing
	bool shouldDrawInvariantMassCanvas 			= false;
	bool shouldDrawInvariantMassCanvasRegion 	= true;
	bool shouldDrawQuantitiesCanvas 			= false;
	bool shouldDrawEfficiencyCanvas 			= false;





	//Check if the name of dir is ok
	if (string(directoryToSave).back() != string("/"))
	{
		cerr << "To avoid errors, please end the result directory with a \"/\"" << endl;
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
	TFile *file0  = TFile::Open(files[useFile]);
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

	//Create a object and set configs
	TagAndProbe TNP;
	TNP.method 			= method;
	TNP.directoryToSave = directoryToSave;

	cout << "Using method " << TNP.method << endl;

	//Get data size and set data limit if has
	long long numberEntries = TreePC->GetEntries();
	if (limitData > 0 && limitData < numberEntries)
		numberEntries = limitData;
	printf("Data analysed = %lld of %lld\n", numberEntries, TreePC->GetEntries());



	//Prepare for showing progress
	string progressFormat = "Progress: %05.2f%% %0"+to_string(strlen(to_string(numberEntries).data()))+"lld/%lld\r";
	auto lastTime = std::chrono::steady_clock::now();

	cout << endl;
	cout << "Filling Invariant Mass Histograms..... (1/2)" << endl;

	//Loop between the components
	for (long long i = 0; i < numberEntries; i++)
	{
		//Select particle pair
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Show progress on screen
		if (chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - lastTime).count() >= 1000 || i == numberEntries - 1)
		{
			printf(progressFormat.data(), (float)(i+1)/(float)numberEntries*100, i+1, numberEntries);
			lastTime = chrono::steady_clock::now();
		}

		//Fill histograms
		if (applyCuts(quantities, types))
		{
			TNP.fillMassHistograms(quantities, types);
		}
	}
	cout << endl << endl;

	//Do function fit ober the histogram
	TNP.doFit();
	
	//Get values for invariant mass and sigma from plot
	TNP.updateMassValuesAll();



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

		TNP.createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}

	if (shouldDrawInvariantMassCanvasRegion)
	{
		bool drawRegions 	= true;
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		TNP.createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}


	
	//Prepare for showing progress
	lastTime = chrono::steady_clock::now();

	cout << endl;
	cout << "Filling Quantities Histograms..... (2/2)" << endl;

	//Loop between the components again
	for (long long i = 0; i < numberEntries; i++)
	{
		//Select particle pair
		TreePC->GetEntry(i);
		TreeAT->GetEntry(i);

		//Show progress on screen
		if (chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - lastTime).count() >= 1000 || i == numberEntries - 1)
		{
			printf(progressFormat.data(), (float)(i+1)/(float)numberEntries*100, i+1, numberEntries);
			lastTime = chrono::steady_clock::now();
		}

		//Fill histograms
		if (applyCuts(quantities, types))
		{	
			TNP.fillQuantitiesHistograms(quantities, types);
		}
	}
	cout << endl << endl;

	//For sideband subtraction
	TNP.subtractSigHistograms();



	if (shouldDrawQuantitiesCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		TNP.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	//Debug consistency for histograms
	TNP.consistencyDebugCout();

	//Save histograms
	generatedFile->mkdir("histograms/");
	generatedFile->   cd("histograms/");

	//Write quantities histograms on file
	{
		bool writehSigBack 	= true;
		bool writehSig 		= true;
		bool writehBack 	= true;
	
		TNP.writeQuantitiesHistogramsOnFile(writehSigBack, writehSig, writehBack);
	}

	//Write mass histograms on file
	{
		bool writehPass = true;
		bool writehAll 	= true;

		TNP.writeMassHistogramsOnFile(writehPass, writehAll);
	}

	//Save plots
	generatedFile->mkdir("efficiency/plots/");
	generatedFile->cd("efficiency/plots/");

	//Creates efficiency plots
	{
		bool shouldWrite 	= true;

		TNP.createEfficiencyPlot(shouldWrite);
	}


	//Saves new histograms and canvas in file
	generatedFile->mkdir("efficiency/canvas/");
	generatedFile->cd("efficiency/canvas/");

	if (shouldDrawEfficiencyCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		TNP.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}

	//Close files
	generatedFile->Close();
}

//Call functions
void macro()
{
	generateHistograms();
}