# Tag & Probe Fitting

> Fitting for tag &amp; probe project
> Teste pelo git

## Note

Link to our [Google Drive](https://drive.google.com/drive/folders/1KZ0OyHnHObX_z6l_ZQ3LN4n7lWHzJ9Fy).

## Necessary files and preparation

The analysed datas are from this file:
* [DoubleMu_data_ntuples.tar]() (It is not the one on Google Drive. I need to replace)

After download the file, you will need to merge all `.root` files  in one `data_histoall.root` file and put it in the same diretory of downloaded code. Here follows how to merge the histogram files. Run in the same directory of hisograms:

```sh
$ hadd data_histoall.root *.root
```

It is necessary to have a folder named `result` on `master` folder side.

## Development setting

It is necessary to have [ROOT](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html), CERN's software, installed on your machine.

Go on your folder where the file code is downloaded and run:

```sh
$ cd master
$ root -l -n
root[0] .L boss.cpp
root[1] boss()
```

## Generated files and images

![](result/InvariantMassPassing.png)

![](result/InvariantMassFailing.png)


![](result/PtPassingProbe.png)

![](result/PtPassingTag.png)

![](result/EtaPassingProbe.png)

![](result/EtaPassingTag.png)

![](result/PhiPassingProbe.png)

![](result/PhiPassingTag.png)


![](result/PtPassingProbe_Efficiency.png)

![](result/PtPassingTag_Efficiency.png)

![](result/EtaPassingProbe_Efficiency.png)

![](result/EtaPassingTag_Efficiency.png)

![](result/PhiPassingProbe_Efficiency.png)

![](result/PhiPassingTag_Efficiency.png)


![](result/PtFailingProbe.png)

![](result/PtFailingTag.png)

![](result/EtaFailingProbe.png)

![](result/EtaFailingTag.png)

![](result/PhiFailingProbe.png)

![](result/PhiFailingTag.png)


![](result/PtFailingProbe_Efficiency.png)

![](result/PtFailingTag_Efficiency.png)

![](result/EtaFailingProbe_Efficiency.png)

![](result/EtaFailingTag_Efficiency.png)

![](result/PhiFailingProbe_Efficiency.png)

![](result/PhiFailingTag_Efficiency.png)

## Output
In addition a .root file is generated named `generated_hist.root` with all canvas above and histograms inside.
