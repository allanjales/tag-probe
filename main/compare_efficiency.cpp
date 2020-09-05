 /*
!--------------------------------
!Purpose: Compare efficiency of files
!--------------------------------	
!author: Allan Jales
!--------------------------------
*/

void compare_plot(TFile *file0, TFile *file1, const char* path)
{
	TEfficiency* pEff0 = (TEfficiency*)file0->Get(path);
	TEfficiency* pEff1 = (TEfficiency*)file1->Get(path);

	if (pEff0 == NULL)
	{
		cerr << "Could not read the path in file0\n";
		abort();
	}
	
	if (pEff1 == NULL)
	{
		cerr << "Could not read the path in file1\n";
		abort();
	}

	//Create canvas
	TCanvas* c1 = new TCanvas();
	gStyle->SetOptTitle(0);
	c1->SetMargin(0.10, 0.03, 0.11, 0.03);

	//Plot
	pEff0->SetMarkerColor(kRed);
	pEff0->SetLineColor(kRed);
	pEff0->Draw();

	pEff1->SetMarkerColor(kBlue);
	pEff1->SetLineColor(kBlue);
	pEff1->Draw("same");
	
	//Set range in y axis
	gPad->Update();
	auto graph = pEff0->GetPaintedGraph(); 
	graph->SetMinimum(0.0);
	graph->SetMaximum(1.2);
	gPad->Update();

	//Set range if is pT
	if (regex_match(path, regex(".*Pt.*")))
	{
		pEff0->GetPaintedGraph()->GetHistogram()->GetXaxis()->SetRangeUser(0.,80.);
		//graph->SetMinimum(0.0);
		graph->SetMaximum(1.4);
	}
	
	//Set range if is pT
	if (regex_match(path, regex(".*Eta.*")))
	{
		pEff0->GetPaintedGraph()->GetHistogram()->GetXaxis()->SetRangeUser(-3.,3.);
		//graph->SetMinimum(0.8);
		graph->SetMaximum(1.3);
	}
	
	//Set range if is pT
	if (regex_match(path, regex(".*Phi.*")))
	{
		graph->SetMinimum(0.40);
		graph->SetMaximum(1.25);
	}

	//Legenda
	TLegend* tl = new TLegend(0.68,0.83,0.94,0.93);
	tl->SetTextSize(0.04);
	tl->AddEntry(pEff0, "Real data",      "lp");
	tl->AddEntry(pEff1, "Simulated data", "lp");
	tl->Draw();

	//CMS Open Data
	TLatex* txCOD = new TLatex();
	txCOD->SetTextSize(0.04);
	txCOD->SetTextAlign(12);
	txCOD->SetTextFont(42);
	txCOD->SetNDC(kTRUE);
	txCOD->DrawLatex(0.14,0.90,Form("#bf{CMS Open Data}"));







	//Saving as png


	//Path where is going to save results 
	const char* directoryToSave = "../result_comparison/";

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

	//Path of file
	string saveAs = string(directoryToSave) + string(pEff0->GetName()) + ".png";

	c1->SaveAs(saveAs.data());
}

//Compare efficiency
void compare_efficiency()
{
	TFile *file0 = TFile::Open("../resultRUN/generated_hist.root");
	TFile *file1 = TFile::Open("../resultMC/generated_hist.root");

	if (file0 == NULL || file1 == NULL)
	{
		std::cerr << "ABORTING...\n";
		abort();
	}

	compare_plot(file0, file1, "efficiency/plots/Muon_Pt_Tracker_Probe_Efficiency");
	compare_plot(file0, file1, "efficiency/plots/Muon_Eta_Tracker_Probe_Efficiency");
	compare_plot(file0, file1, "efficiency/plots/Muon_Phi_Tracker_Probe_Efficiency");;
}