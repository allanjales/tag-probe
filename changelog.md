# Changelog
> This is for internal porpuse.

* 2020-08-09
	* Initial parameters set to pointers at declaration in all classes;
	* Changed to standard pointer form (`const char *example` -> `const char* example`);
	* lowercase `const char* PassingOrFailing` -> `const char* passingOrFailing`;
	* Moved ROOT types from `Double_t`, `Int_t` and `Long64_t` for C++ types `double`, `int` and `long long`;
	* Changed the name of function `createAllMassHistograms()` in invariantMass.h to `createPassingAndAllMassHistograms()`;
	* Added argument in `createMassHistogram(...)`: `bool alertIfCant = true`;
	* Changed function `createMassHistogram(...)` and `createPassingAndAllMassHistograms()`;
	* Changed inital parameters on `doFit()` in invariantMass.h;
	* Changed macro.cpp to show progress based on time.

* 2020-08-10
	* Changed `Histogram` class to `PtEtaPhi` class. Also Histogram.h -> PtEtaPhi.h;
	* Added `int method` variable on macro.cpp;
	* Changed argument in `createCanvas(...)` function in InvariantMass class and PtEtaPhi class: `shouldSave` -> `shouldSavePNG`;
	* Added argument in `createCanvas(...)` in InvariantMass class and PtEtaPhi class: `const char* directoryToSave`;
	* Added `abort()` function in macro.cpp if could not find or create the directory;
	* Changed text on failed creating directory output;
	* Added a function to check if the directory on macro.cpp was set right;
	* Added list of pointers to variables from readed file on macro.cpp: `double* quantities[7]` and `int* types[3]`;
	* Added file cutsAndFill.h and included in macro.cpp;
	* Added function `applyCuts(...)` in cutsAndFill.h;
	* Added function `fillMassHistograms(..)` in cutsAndFill.h;
	* Added function `fillQuantitiesHistograms(..)` in cutsAndFill.h;
	* Changed both loops on macro.cpp to use new functions;
	* Changed fatal errors in macroc.cpp to output with `cerr`, not `cout`. 
