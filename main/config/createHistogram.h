/*
What is this file: This file is called in PassingFailing.cpp
What does it do:   It set quantity histograms bins and create the hitogram
*/

void createHistogram(TH1D* &histo, const char* histoName)
{
	//Set parameters
	string hName 		  = string(particleType) + string(passingOrFailing) + string(tagOrProbe) + string(particleName) + "_" + string(quantityName) + string(histoName);
	string hTitle 		  = string(passingOrFailing) + " in " + string(particleType) + " " + string(tagOrProbe);
	string xAxisTitle 	  = string(xAxisName);
	string yAxisTitleForm = "Events";

	//Add unit if has
	if (strcmp(quantityUnit, "") != 0)
		xAxisTitle += " [" + string(quantityUnit) + "]";

	//Change title is passing
	if (strcmp(passingOrFailing, "Passing") == 0)
		hTitle = string(particleType) + " " + string(particleName) + " " + string(tagOrProbe);

	if (strcmp(passingOrFailing, "All") == 0)
		hTitle = "All " + string(particleName) + " " + string(tagOrProbe);


	//Variable bin for pT
	if (strcmp(quantityName, "Pt") == 0)
	{
		double xbins[10000];
		xbins[0] = .0;
		int nbins = 0;
		double binWidth = 1.;
		for (int i = 1; xbins[i-1] < xMax+binWidth; i++)
		{
			xbins[i] = xbins[i-1] < 1. ? 1. : xbins[i-1] *(1+binWidth);
			nbins++;
		}

		histo = new TH1D(hName.data(), hTitle.data(), nbins, xbins);
	}

	//Variable bin for eta
	else if (strcmp(quantityName, "Eta") == 0)
	{
		double xbins[10000];
		xbins[0] = .5;
		int nbins = 0;
		double binWidth = 0.2;

		//For positive
		for (int i = 1; xbins[i-1] < xMax+binWidth; i++)
		{
			xbins[i] = xbins[i-1] < 1. ? 1. : xbins[i-1] *(1+binWidth);
			nbins++;
		}

		//Duplicate array and create another
		double rxbins[nbins*2+1];
		int entry = 0;
		for (int i = nbins; i >= 0; i--)
		{
			rxbins[entry] = -xbins[i];
			entry++;
		}
		rxbins[entry] = 0.;
		entry++;
		for (int i = 0; i <= nbins; i++)
		{
			rxbins[entry] = xbins[i];
			entry++;
		}
		
		histo = new TH1D(hName.data(), hTitle.data(), entry-1, rxbins);
	}

	//Bins for phi 
	else
	{
		if (strcmp(quantityUnit, "") == 0)
		{
			yAxisTitleForm += " / (%1." + to_string(decimals) + "f)";
		}
		else
		{
			yAxisTitleForm += " / (%1." + to_string(decimals) + "f " + string(quantityUnit) + ")";
		}

		histo = new TH1D(hName.data(), hTitle.data(), nBins, xMin, xMax);
	}

	//Edit histogram axis
	histo->GetYaxis()->SetTitle(Form(yAxisTitleForm.data(), histo->GetBinWidth(0)));
	histo->GetXaxis()->SetTitle(xAxisTitle.data());
}