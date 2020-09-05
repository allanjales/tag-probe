//Contains pass or all fit infos
struct MassValues
{
	//Histogram and fit function
	TH1D* hMass 	    = NULL;
	TF1*  fitFunction   = NULL;
	TF1*  fitSignal     = NULL;
	TF1*  fitBackground = NULL;

	//Regions
	double sidebandRegion1_x1  = 0.;
	double sidebandRegion1_x2  = 0.;
	double signalRegion_x1     = 0.;
	double signalRegion_x2     = 0.;
	double sidebandRegion2_x1  = 0.;
	double sidebandRegion2_x2  = 0.;

	//For downgrade
	TFitResultPtr fitResult = 0;

	//--- For sideband subtraction ---

	bool isInSignalRegion(double InvariantMass)
	{
		if (InvariantMass >= this->signalRegion_x1 && InvariantMass <= this->signalRegion_x2)
			return true;

		return false;
	}

	bool isInSidebandRegion(double InvariantMass)
	{
		if ((InvariantMass >= this->sidebandRegion1_x1 && InvariantMass <= this->sidebandRegion1_x2) ||
			(InvariantMass >= this->sidebandRegion2_x1 && InvariantMass <= this->sidebandRegion2_x2))
			return true;

		return false;
	}

	double subtractionFactor()
	{
		//Simple method (for linear background)
		double signalRegion = abs(signalRegion_x2 - signalRegion_x1);
		double sidebandRegion = abs(sidebandRegion1_x2 - sidebandRegion1_x1) + abs(sidebandRegion2_x2 - sidebandRegion2_x1);

		//Using yield (advanced method)
		if (fitFunction != NULL)
		{		
			signalRegion    = fitBackground->Integral(signalRegion_x1,    signalRegion_x2)   /hMass->GetBinWidth(0);
			sidebandRegion  = fitBackground->Integral(sidebandRegion1_x1, sidebandRegion1_x2)/hMass->GetBinWidth(0);
			sidebandRegion += fitBackground->Integral(sidebandRegion2_x1, sidebandRegion2_x2)/hMass->GetBinWidth(0);
		}
		else
		{
			cerr << "WARNING: not using advanced method for subtraction factor calculation. Using method for linear background." << endl;
		}

		return signalRegion/sidebandRegion;
	}

	void doFitJpsi()
	{
		TF1* &f  	= fitFunction;
		TF1* &fs 	= fitSignal;
		TF1* &fb 	= fitBackground;

		double xMin = hMass->GetXaxis()->GetXmin();
		double xMax = hMass->GetXaxis()->GetXmax();

		//Temporary for test
		xMin = 2.9;
		xMax = 3.3;

		const char* const fittingParName[] = {
				"Gaus(Sg) Height  ",
				"Gaus(Sg) Position",
				"Gaus(Sg) Sigma   ",

				"CB  (Sg) Alpha   ",
				"CB  (Sg) N       ",
				"CB  (Sg) Mean    ",
				"CB  (Sg) Sigma   ",
				"CB  (Sg) Yield   ",

				"Exp1(Bg) Height  ",
				"Exp1(Bg) Width   ",
				"Exp2(Bg) Height  ",
				"Exp2(Bg) Width   "
			};

		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);

		//Fit Function
		f = new TF1("FitFunction", FitFunctions::Upsilon::InvariantMass, xMin, xMax, arraySize);
		f->SetNpx(1000);
		f->SetLineStyle(kSolid);
		f->SetLineColor(kBlue);
		f->SetLineWidth(3);

		//Rename parameters
		double resultParameters[arraySize];
		for (int i = 0; i < arraySize; i++)
		{
			f->SetParName(i, fittingParName[i]);
		}
		
		//Values Signal GS
		f->SetParameter(0,	4098.2);
		f->SetParameter(1,	3.09);
		f->SetParameter(2,	0.020);

		//Values Signal CB
		f->SetParameter(3,	1.58);
		f->SetParameter(4,	1.54);
		f->SetParameter(5,	3.093);
		f->SetParameter(6,	0.032);
		f->SetParameter(7,	42022.27);

		//Values Background
		f->SetParameter(8,	-0.217);
		f->SetParameter(9,	1.915);
		f->SetParameter(10, 263.185);
		f->SetParameter(11,	0.061);

		//Fit function
		fitResult = hMass->Fit(f, "RNS", "", xMin, xMax);
		f->GetParameters(resultParameters);
		
		//Signal Fitting
		fs = new TF1("FitFunction_Signal", FitFunctions::Jpsi::Signal_InvariantMass, xMin, xMax, 8);
		fs->SetNpx(1000);							//Resolution of signal fit function
		fs->SetParameters(resultParameters);		//Get only signal part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style
		fs->SetLineWidth(3);						//Fit width

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Jpsi::Background_InvariantMass, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[8]);	//Get only background part
		fb->SetLineColor(kBlue); 					//Fit Color
		fb->SetLineStyle(kDashDotted);				//Fit style
		fb->SetLineWidth(3);						//Fit width
		
		/*
		//TEST
		cout << "Entries  (TH1): " << hMass->GetEntries() << endl;
		cout << "Integral (TH1): " << hMass->Integral(0, hMass->GetNbinsX()+1) << endl;
		cout << "Integral (TF1): " << f->Integral(xMin, xMax)/hMass->GetBinWidth(0) << endl;
		*/
		cout << "chi2/ndf = " << (this->fitResult)->Chi2()/(this->fitResult)->Ndf() << "\n";
	}

	void doFitUpsilon()
	{
		cout << fixed;
		TF1* &f  	= fitFunction;
		TF1* &fs 	= fitSignal;
		TF1* &fb 	= fitBackground;

		double xMin = hMass->GetXaxis()->GetXmin();
		double xMax = hMass->GetXaxis()->GetXmax();

		const char* const fittingParName[] = {
				"CB  (1S) Alpha   ",
				"CB  (1S) N       ",
				"CB  (1S) Mean    ",
				"CB  (1S) Sigma   ",
				"CB  (1S) Yield   ",

				"Gaus(2S) Height  ",
				"Gaus(2S) Position",
				"Gaus(2S) Sigma   ",

				"Gaus(3S) Height  ",
				"Gaus(3S) Position",
				"Gaus(3S) Sigma   ",

				"Pol3(Bg) a       ",
				"Pol3(Bg) b       ",
				"Pol3(Bg) c       ",
				"Pol3(Bg) d       "
			};

		int arraySize = sizeof(fittingParName)/sizeof(*fittingParName);

		//Fit Function
		f = new TF1("FitFunction", FitFunctions::Upsilon::InvariantMass, xMin, xMax, arraySize);
		f->SetNpx(1000);
		f->SetLineStyle(kSolid);
		f->SetLineColor(kBlue);
		f->SetLineWidth(3);

		//Rename parameters
		double resultParameters[arraySize];
		for (int i = 0; i < arraySize; i++)
		{
			f->SetParName(i, fittingParName[i]);
		}

		//Values Y(1S)
		f->SetParameter(0, -10.);
		f->SetParameter(1,  9.4603);
		f->SetParameter(2,  9.4);
		f->SetParameter(3, -0.08);
		f->SetParameter(4,  4673.);

		//Values Y(2S)
		f->SetParameter(5,  374.);
		f->SetParameter(6,  10.02326);
		f->SetParameter(7,  0.09);

		//Values Y(3S)
		f->SetParameter(8,  237.);
		f->SetParameter(9,  10.3552);
		f->SetParameter(10, 0.08);

		//Values Background
		f->SetParameter(11, -193041.);
		f->SetParameter(12,  57792.);
		f->SetParameter(13, -5729.);
		f->SetParameter(14,  190.);

		//Set par limits
		f->SetParLimits(1, xMin, xMax);
		f->SetParLimits(6, xMin, xMax);
		f->SetParLimits(7, 0., xMax-xMin);
		f->SetParLimits(9, xMin, xMax);
		f->SetParLimits(10, 0., xMax-xMin);

		//Fit
		fitResult = hMass->Fit(f, "RNS", "", xMin, xMax);
		f->GetParameters(resultParameters);

		//Signal of analysis Fitting
		fs = new TF1("FitFunction_CrystalBall", FitFunctions::Primary::CrystalBall, xMin, xMax, 5);
		fs->SetNpx(1000);							//Resolution of signal fit function
		fs->SetParameters(resultParameters);		//Get only signal part
		fs->SetLineColor(kMagenta); 				//Fit Color
		fs->SetLineStyle(kSolid);					//Fit Style
		fs->SetLineWidth(3);						//Fit width

		//Background Fitting
		fb = new TF1("FitFunction_Background", FitFunctions::Upsilon::Background_InvariantMass, xMin, xMax, 4);
		fb->SetNpx(1000);							//Resolution of background fit function
		fb->SetParameters(&resultParameters[11]);	//Get only background part
		fb->SetLineColor(kBlue); 					//Fit Color
		fb->SetLineStyle(kDashDotted);				//Fit style
		fb->SetLineWidth(3);						//Fit width

		cout << "chi2/ndf = " << (this->fitResult)->Chi2()/(this->fitResult)->Ndf() << "\n";
	}

	TBox* createTBox(double Ymax, int index = 0, double Ymin = 0.)
	{
		//index = -1 -> left region
		//index = 0 -> signal region
		//index = 1 -> right region

		double x1, x2 = 0;

		switch(index)
		{
			case -1:
				x1 = sidebandRegion1_x1;
				x2 = sidebandRegion1_x2;
				break;
			case 0:
				x1 = signalRegion_x1;
				x2 = signalRegion_x2;
				break;
			case 1:
				x1 = sidebandRegion2_x1;
				x2 = sidebandRegion2_x2;
				break;
		}

		return new TBox(x1, Ymin, x2, Ymax);
	}
};