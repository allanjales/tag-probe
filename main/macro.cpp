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
void macro()
{
	//List of files
	const char *files[] = {"../data_histoall.root",
							"../Run2011AMuOnia_mergeNtuple.root",""
							"../JPsiToMuMu_mergeMCNtuple.root",
							"../Run2011A_MuOnia_Upsilon.root",
							"../Upsilon1SToMuMu_MC_full.root"};

	const char* directoriesToSave[] = {"../results/result/",
										"../results/Jpsi Run 2011/",
										"../results/Jpsi MC 2020/",
										"../results/Upsilon Run 2011/",
										"../results/Upsilon MC 2020/"};

	//Options to change

	//Which file of files (variable above) should use
	int useFile = 3;

	//Choose method
	//if 1 -> sideband by histogram || if 2 -> sideband by fitting
	int method = 2;

	//Set the canvasW wtermark
	const char* canvasWatermark = "#bf{CMS Open Data}";

	//Path where is going to save results 
	const char* directoryToSave = directoriesToSave[useFile];
	directoryToSave = "../result/";

	//Should limit data?
	long long limitData = 0; //0 -> do not limit

	//Canvas drawing
	bool shouldDrawInvariantMassCanvas 			= true;
	bool shouldDrawInvariantMassCanvasRegion 	= true;
	bool shouldDrawQuantitiesCanvas 			= true;
	bool shouldDrawEfficiencyCanvas 			= true;

	//Muon id anlyse
	bool doTracker    = true;
	bool doStandalone = true;
	bool doGlobal     = true;

	//Auto detect resonance due file index
	const char* resonance = "Jpsi";
	if (useFile > 2)
		resonance = "Upsilon";

	//For saving in log file
	//freopen((string(directoryToSave) + "log.txt").data(), "w", stdout);
	//freopen((string(directoryToSave) + "log.txt").data(), "w", stderr);

	//Auto detect if is MC due file index
	bool isMC = false;
	if (useFile == 2 || useFile == 4)
		isMC = true;


	//Check if the name of dir is ok
	if (string(directoryToSave).back() != string("/"))
	{
		cerr << "To avoid errors, please end the result directory with a \"/\"" << endl;
		abort();
	}

	//Check if dir exists and create
	if (gSystem->AccessPathName(directoryToSave))
	{
		if (gSystem->mkdir(directoryToSave, true))
		{
			cerr << "\"" << directoryToSave << "\" path could not be found and could not be created ERROR" << endl;
			cerr << "Try to create manually this folder path" << endl;
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

	double* quantities[] = {&ProbeMuon_Pt,
							&ProbeMuon_Eta,
							&ProbeMuon_Phi,
							&TagMuon_Pt,
							&TagMuon_Eta,
							&TagMuon_Phi,
							&InvariantMass,
		};
	
	int* types[] = {&PassingProbeTrackingMuon,
					&PassingProbeStandAloneMuon,
					&PassingProbeGlobalMuon
		};

	//Create a object and set configs
	TagAndProbe TNP{resonance};
	TNP.method 			= method;
	TNP.canvasWatermark	= canvasWatermark;
	TNP.directoryToSave = directoryToSave;
	TNP.doTracker       = doTracker;
	TNP.doStandalone    = doStandalone;
	TNP.doGlobal        = doGlobal;

	cout << "resonance: " << TNP.resonance << endl;
	cout << "Using method " << TNP.method << endl;

	//Get data size and set data limit if has
	long long numberEntries = TreePC->GetEntries();
	if (limitData > 0 && limitData < numberEntries)
		numberEntries = limitData;
	printf("Data analysed = %lld of %lld\n", numberEntries, TreePC->GetEntries());



	//Prepare for showing progress
	string progressFormat = "Progress: %05.2f%% %0"+to_string(strlen(to_string(numberEntries).data()))+"lld/%lld\r";
	auto lastTime = std::chrono::steady_clock::now();
	auto start    = std::chrono::steady_clock::now();

/*
	//TEMPORARY FOR TEST
	TFile *file1 = TFile::Open("../results/Jpsi Run 2011/generated_hist.root");
	TNP.Tracker.Mass.Pass.hMass = (TH1D*)file1->Get("histograms/Passing_Tracker_Muon_InvariantMass");
	TNP.Tracker.Mass.All.hMass = (TH1D*)file1->Get("histograms/All_Tracker_Muon_InvariantMass");
	TNP.Standalone.Mass.Pass.hMass = (TH1D*)file1->Get("histograms/Passing_Standalone_Muon_InvariantMass");
	TNP.Standalone.Mass.All.hMass = (TH1D*)file1->Get("histograms/All_Standalone_Muon_InvariantMass");
	TNP.Global.Mass.Pass.hMass = (TH1D*)file1->Get("histograms/Passing_Global_Muon_InvariantMass");
	TNP.Global.Mass.All.hMass = (TH1D*)file1->Get("histograms/All_Global_Muon_InvariantMass");
*/

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

	cout << "\nTook " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count() << " ms" << endl;

	//Do function fit ober the histogram
	if (!isMC)
		TNP.doFit();

	//Get values for invariant mass and sigma from plot
	if (!isMC)
		TNP.updateMassValuesAll();

	//-------------------------------------
	// Generate and save files
	//-------------------------------------

	//Supress canvas PROBLEMS BELOW:
	//- Does not store TBox of
	//- Does not save anything in created .root
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

	if (shouldDrawInvariantMassCanvasRegion && !isMC)
	{
		bool drawRegions 	= true;
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		TNP.createMassCanvas(drawRegions, shouldWrite, shouldSavePNG);
	}
	
	//Prepare for showing progress
	lastTime = std::chrono::steady_clock::now();
	start    = std::chrono::steady_clock::now();

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
			TNP.fillQuantitiesHistograms(quantities, types, isMC);
		}
	}
	cout << "\nTook " << chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count() << " ms" << endl;

	//For sideband subtraction
	if (!isMC)
		TNP.subtractSigHistograms();



	if (shouldDrawQuantitiesCanvas)
	{
		bool shouldWrite 	= true;
		bool shouldSavePNG 	= true;

		cout << endl;
		TNP.createQuantitiesCanvas(shouldWrite, shouldSavePNG);
	}

	//Debug consistency for histograms
	if (!isMC)
		TNP.consistencyDebugCout();
	else
		cout << "No consistency check needed. It's MC simulations.\n";

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

		cout << endl;
		TNP.createEfficiencyCanvas(shouldWrite, shouldSavePNG);
	}

	//Close files
	generatedFile->Close();

	cout << "\nDone. All result files can be found at \"" << TNP.directoryToSave << "\"\n" << endl;
}