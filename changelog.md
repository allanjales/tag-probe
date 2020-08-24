# Changelog
> This is for internal porpuse.

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