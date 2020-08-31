 /*
!--------------------------------
!Purpose: Compare efficiency of files
!--------------------------------	
!author: Allan Jales
!--------------------------------
*/

void compare_efficiency()
{
	TFile *file0  = TFile::Open("../result/jpsi_run2011.root");
	TFile *file1  = TFile::Open("../aresult/jpsi_mc2020.root");

	if (file0 == NULL || file1 == NULL)
	{
		std::cerr << "ABORTING...\n";
		abort();
	}

	TEfficiency* pEff = (TEfficiency*)file0->Get("efficiency/plots/Muon Probe Pt Tracker Efficiency");
}