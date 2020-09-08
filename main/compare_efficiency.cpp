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

	int colorScheme[][2] = {
		{kGreen - 2, kBlue},
		{kBlue,      kRed},
		{kGreen - 2, kRed}
	};

	const char* nameScheme[][2] = {
		{"#Upsilon data",      "J/#psi data"},
		{"J/#psi real data",   "Simulated data"},
		{"#Upsilon real data", "Simulated data"}
	};

	const char* directoriesToSave[] = {
		"../Comparison Jpsi vs Upsilon/",
		"../Comparison Jpsi Run vs MC/",
		"../Comparison Upsilon Run vs MC/",
	};

	int useScheme = 0;
	//Upsilon vs Jpsi
	//Jpsi    Run vs MC
	//Upsilon Run vs MC

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
	//gStyle->SetOptTitle(0);
	c1->SetMargin(0.10, 0.03, 0.11, 0.07);

	//Plot
	pEff0->SetMarkerColor(colorScheme[useScheme][0]);
	pEff0->SetLineColor(colorScheme[useScheme][0]);
	pEff0->Draw();

	pEff1->SetMarkerColor(colorScheme[useScheme][1]);
	pEff1->SetLineColor(colorScheme[useScheme][1]);
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
		graph->SetMinimum(0.5);
		graph->SetMaximum(1.15);
	}
	
	//Set range if is pT
	if (regex_match(path, regex(".*Eta.*")))
	{
		pEff0->GetPaintedGraph()->GetHistogram()->GetXaxis()->SetRangeUser(-3.,3.);
		graph->SetMinimum(0.8);
		graph->SetMaximum(1.05);
	}
	
	//Set range if is pT
	if (regex_match(path, regex(".*Phi.*")))
	{
		graph->SetMinimum(0.8);
		graph->SetMaximum(1.05);
	}

	//Legenda
	TLegend* tl = new TLegend(0.68,0.78,0.94,0.88);
	tl->SetTextSize(0.04);
	tl->AddEntry(pEff0, nameScheme[useScheme][0], "lp");
	tl->AddEntry(pEff1, nameScheme[useScheme][1],   "lp");
	tl->Draw();

	//CMS Open Data
	TLatex* txCOD = new TLatex();
	txCOD->SetTextSize(0.04);
	txCOD->SetTextAlign(12);
	txCOD->SetTextFont(42);
	txCOD->SetNDC(kTRUE);
	txCOD->DrawLatex(0.14,0.85,Form("#bf{CMS Open Data}"));







	//Saving as png


	//Path where is going to save results 
	const char* directoryToSave = directoriesToSave[useScheme];

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
	TFile *file0 = TFile::Open("../Upsilon Run 2011/generated_hist.root");
	TFile *file1 = TFile::Open("../Jpsi Run 2011/generated_hist.root");

	if (file0 == NULL || file1 == NULL)
	{
		std::cerr << "ABORTING...\n";
		abort();
	}

	compare_plot(file0, file1, "efficiency/plots/Muon_Pt_Tracker_Probe_Efficiency");
	compare_plot(file0, file1, "efficiency/plots/Muon_Eta_Tracker_Probe_Efficiency");
	compare_plot(file0, file1, "efficiency/plots/Muon_Phi_Tracker_Probe_Efficiency");;
}