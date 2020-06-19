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

#include "classes/Particle.h"

using namespace std;

void openFile()
{
	//Create object
	Particle Muon;
	Muon.setMethod(1);

	int useData = 0;

	//Open and read files
	TFile *file0 = TFile::Open("resultsForTesting/run2011.root");
	
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