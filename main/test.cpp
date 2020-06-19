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

#include "classes/Particle.h"

using namespace std;

void openFile()
{
	//Create object
	Particle Muon;
	Muon.setMethod(1);

	const char *files[3] = {"raphael.root",
							"run2011.root",
							"mc2020.root"};

	int useFile = 1;

	//Open and read files
	TFile *file0 = TFile::Open((string("resultsForTesting/") + string(files[useFile])).data());
	
	Muon.Pass.Mass.hMass = (TH1D*)file0->Get("histograms/PassingMuonInvariantMass");
	Muon.Fail.Mass.hMass = (TH1D*)file0->Get("histograms/FailingMuonInvariantMass");
	Muon.Both.Mass.hMass = (TH1D*)file0->Get("histograms/AllMuonInvariantMass");

	Muon.doFit();

	Muon.Pass.Mass.createCanvas();
	Muon.Fail.Mass.createCanvas();
	Muon.Both.Mass.createCanvas();
}

//Call functions
void test()
{
	openFile();
}