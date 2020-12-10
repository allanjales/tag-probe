/*
What is this file: This file is called in macro.cpp
What does it do:   It set some main options in this code
*/

//List of files
const char *files[] = {"../data_histoall.root",
						"../Run2011AMuOnia_mergeNtuple.root",
						"../JPsiToMuMu_mergeMCNtuple.root",
						"../Run2011A_MuOnia_Upsilon.root",
						"../Upsilon1SToMuMu_MC_full.root"};

const char* directoriesToSave[] = {"../results/result/",
									"../results/Jpsi_Run_2011/",
									"../results/Jpsi_MC_2020_sbs/",
									"../results/Upsilon_Run_2011/",
									"../results/Upsilon_MC_2020_sbs/"};

//MAIN OPTIONS

//Which file of files (variable above) should use
int useFile = 1;

//Set the canvasW wtermark
const char* canvasWatermark = "#bf{CMS Open Data}";

//Path where is going to save results 
const char* directoryToSave = directoriesToSave[useFile];
//directoryToSave = "../result/";

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

//ENDED MAIN OPTIONS



//Auto detect resonance due file index
const char* resonance = "Jpsi";
if (useFile > 2)
	resonance = "Upsilon";

//Auto detect if is MC due file index
bool isMC = false;
if (useFile == 2 || useFile == 4)
	isMC = true;
isMC = false;

//Auto detect limit of data
if (limitData > 0)
	directoryToSave = "../partial_result/";



//Compatibility adjusts on file read (for data_histoall ntupples)
bool needsRetroCompatibility = false;

if (useFile == 0)
	needsRetroCompatibility = true;