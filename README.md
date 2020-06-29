# Tag & Probe Fitting

> Tag &amp; probe efficiency calculus project

## Necessary files and preparation

The analysed datas are from those files:
* [0] [DoubleMu_data_ntuples.tar]() - obsolete (Old ntupple. Need to be merged. Current name after merge should be `data_histoall.root`)
* [1] [Run2011AMuOnia_mergeNtuple.root]()
* [2] [JPsiToMuMu_mergeMCNtuple.root]()

After download one of those files, you can run the code. Don't forget to use `useNewData` (int) var in `macro.cpp` file to set which ntupple did you choose. Just set the current id above of it.

It is necessary to have a folder named `result` on `main` folder side.

## Preferences

You can change the method to estimate signal region by modifying `Muon.setMethod(1)` line by choosing 1 (estimate by FWHM of histograms) or 2 (estimate by FWHM of fitting):

```cpp
Muon.setMethod(1);
```

Change this line to specify the ntupple you are analysing by choosing 0 (Old ntupple), 1 (Run 2011 ntupple) or 2 (Monte Carlo ntupple):

```cpp
int useNewData = 0;
```

## Development setting

It is necessary to have [ROOT](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html), CERN's software, installed on your machine.

Go on your folder where the file code is downloaded and run:

```sh
$ cd main
$ root -l -n
root[0] .L macro.cpp+
root[1] macro()
```

## Output
Images are created in `result` folder. In addition a .root file is generated named `generated_hist.root` with all canvas and histograms.
