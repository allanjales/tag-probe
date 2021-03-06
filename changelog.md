# Changelog
> This is for internal porpuse. (In macro.cpp structure)

For v1

* 2020-08-19
	* Added signal and sideband region on invariant mass canvas;
	* Changed legends on invariant mass canvas;

* 2020-08-20
	* Legend now shows above result histogram on InvariantMass class;
	* Now macro.cpp shows method as output;
	* Fixed method 2 bug;
	* Changed to user input full path in const char \*files[] on macro.cpp;

* 2020-08-21
	* Changed progress format on macro.cpp;
	* Changed consitencyDebugCout() format and corrected it;
	* Removed `std::` from chrono uses;
	* Changed from `high_resolution_clock` to `steady_clock` in chrono uses on macro.cpp;
	* Updated README.md;
	* Changed data analysed cout form;

* 2020-08-21
	* Added Upsilon files names on macro.cpp;
	* Now MassValues stores signal function;
	* Changed `updateMassValuesFor(...)`. Now do not need create a function for signal, just read;
	* Changed `doFit()` to store signal fit in MassValue structure;

* 2020-08-21 (upsilon)
	* Added function `defineMassHistogramNumbers(...)` in TagAndProbe, Type and InvariantMass class;
	* Changed macro.cpp to call `defineMassHistogramNumbers(...)` to set upsilon histograms
	* Added varaible `const char* ressonance` in TagAndProbe class
	* Added pointer to `ressonance` in Type and InvariantMass class
	* Changed constructor of InvariantMass class
	* Added time to completion output after loops in macro.cpp

* 2020-08-21 (pointer-to-reference)
	* Changed all `int* method` to `int& method`;
	* Changed all `const char** ressonance` to `const char*& ressonance`;
	* Changed all `const char** particleName` to `const char*& particleName`;
	* Changed all `const char** directoryToSave` to `const char*& directoryToSave`;
	* Changed all `const char** particleType` to `const char*& particleType`;
	* Changed all `const char** tagOrProbe` to `const char*& tagOrProbe`;
	* Changed all `const char** quantityName` to `const char*& quantityName`;
	* Changed all `const char** xAxisName` to `const char*& xAxisName`;
	* Changed all `const char** quantityUnit` to `const char*& quantityUnit`;
	* Changed all `const char** extendedQuantityName` to `const char*& extendedQuantityName`;
	* Changed all `double* xMin` to `double& xMin`;
	* Changed all `double* xMax` to `double& xMax`;
	* Changed all `int* nBins` to `int& nBins`;
	* Changed all `int* decimals` to `int& decimals`;
	* Changed all `double** quantities` to `double*& quantities`;
	* Changed all `int** types` to `int*& types`;
	* Changed all `double* quantity` to `double& quantity`;
	* Changed all `double* InvariantMass` to `double& InvariantMass`;
	* Changed all `int* isPassing` to `int& isPassing`
	* Changed all `double* quantity` to `double& quantity`;
	* Changed all `double* InvariantMass` to `double& InvariantMass`;
	* Changed all `double* isPassing` to `double& isPassing`;
	* Changed `updateMassValuesFor(...)` in InvariantMass class;
	* Removed M_JPSI, W_JPSI, signalRegion, sidebandRegion variables in MassValues struct;
	* Added sidebandRegion1_x1, sidebandRegion1_x2, signalRegion_x1, signalRegion_x2, sidebandRegion2_x1, sidebandRegion2_x2 variables in MassValues struct;

* 2020-08-26
	* Added `doTracker`, `doStandalone` and `doGlobal` variables to TagAndProbe class;
	* Corrected `subtractionFactor(...)`;
	* Added error message in `PassFailObj()` in PassingFailing class;
	* Added `Primary::Merged::Both_Background_InvariantMass(...)` in FitFunctions class;
	* Added `TF1* fitBackground` in MassValues struct;
	* Corrected cout of `*particleType` to `particleType` at `doFit()` cout in InvariatMass class;
	* Added advanced method in `subtractionFactor()` calculus in MassValue class (but not using);
	* Changed `drawCanvasQuarter(...)` in InvariantMass class to show background fit for All fit too (removed);
	* Now it prepares for upsilon depending which file index did you choose in macro.cpp;
	* Changed histgogram names of PassingFailing class;
	* Changed result histogram names of PtEtaPhi and PassingFailing class;
	* Changed TLegend of histogram in PassingFailing class changed. All -> Total;
	* Added Ymin argument in `createTBox(...)` in MassValue struct;
	* Changed `drawCanvasQuarter(...)` in InvariantMass class to pass Ymin parameter at creating TBoxes;
	* Corrected `consistencyDebugCout()` in PassingFailing class;
	* Corrected `consistencyDebugCout()` in TagAndProbe class;
	* Corrected `subtractionFactor()` type in MassValues struct. int -> double;

* 2020-08-27
	* Changed hTitle at `createHistogram(...)` in PassingFailing class if is passing;
	* Changed `createDividedCanvas(...)` to show signal only on the right plot;
	* Changed hTitleRight at `createDividedCanvas(...)` in PassingFailing class if is passing;
	* Added `const char* watermark` in TagnAndProbe class;

* 2020-08-28
	* Removed singal counter in left histogram of canvas in PassingFailing class;
	* Removed `float margin[4]` at `drawCanvasQuarter(...)` in InvariantMass class;
	* Removed background fit from canvas at `drawCanvasQuarter(...)` in InvariantMass class;
	* Changed `const char* watermark` to `const char* canvasWatermark` in every file;
	* Removed `gROOT->SetBatch(0);` at `createEfficiencyCanvas(...)` in PassingFailing class;
	* Added canvas supression in macro.cpp;
	* Fixed error at `createEfficiencyCanvas(...)` in PtEtaPhi class with canvas supression;
	* Changed cout of macro.cpp. Added `cout << endl;` before saving canvas;
	* Added end message in macro.cpp;

* 2020-08-28
	* Changed plot name at `createEfficiencyPlot(...)` in PtEtaPhi class;
	* Corrected name of pT canvas: Transversal -> Transverse;

* 2020-09-01
	* Removed `macro()` function in macro.cpp;
	* Renamed `generateHistograms()` to `macro()`;
	* Removed some functions in FitFunction class;
	* Changed `doFit()` function in InvariantMass class. Added this function to MassValue too;
	* Added some variables in MassValue class; 
	* Changed range of pT histograms to 40 GeV max;
	* Changed margins in invariant mass canvas;
	* Corrected Efficiency plot created by signal histograms;
	* Edited error messages in efficiency consistency