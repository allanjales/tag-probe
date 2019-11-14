# Tag & Probe Fitting

> Fitting for tag &amp; probe project

## Note

Link to our [Google Drive](https://drive.google.com/drive/folders/1KZ0OyHnHObX_z6l_ZQ3LN4n7lWHzJ9Fy)

## Arquivos necessários

The analysed datas are from this file:
* [DoubleMu_data_ntuples.tar](https://drive.google.com/file/d/1z4oNmr3Vcv2JOtH-iBxXOFuWCd4llTNe/view?usp=sharing)

After download the file, you will need to merge all `.root` files  in one `.root` file and put it in the same diretory of downloaded files.

## Development setting

It is necessary to have [ROOT](https://root.cern.ch/root/html534/guides/users-guide/InstallandBuild.html), CERN's software, installed on your machine.

Go on your folder where the file is downloaded and run:

```sh
$ root -l -n
root[0] .x bigboss.C
```

## Generated image

![](InvariantMassProbe.png)

## Mudança na versão

* INTEGRACAO
    * Melhoria nos comentários
    * ADD: Um fit separado para cada sinal `f1`,`f2`,`f3`
    * ADD: Área dos ajuste no terminal
