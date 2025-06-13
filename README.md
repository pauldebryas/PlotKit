# PlotKit

## How to install

Clone the repository:
```sh
git clone git@github.com:pauldebryas/PlotKit.git
```

The main conda environment needed to produce the plot:
```sh
conda env create -f common/config/all/env_PlotKit.yml
```

An other conda environment is needed to compute DNN score:
```sh
conda env create -f common/config/all/env_MLenv.yml
```

Don't forget to change the path to your miniconda in env_PlotKit.sh and env_MLenv.sh

This code is using law (luigi analysis workflow) to produce plots (see documentation) so first thing to do is to index the task:
```sh
source env_plot.sh
law index --verbose
```
If you decide to implement new task, you need to redo this step.

## How to run

load environment 
```sh
source env_PlotKit.sh
```

- Parameters:
    Here is the list of some global parameters and a short description:
    - period: period of data taking: can be 2018, 2017, 2016 and 2016_HIPM
    - channel: channel you want to plot: can be ttm, tmm, tte, tee and tem
    - tag = tag used for anatuple production

- Useful informations:
    - All regions used to do background estimation are defined in common/regions/regions.py
    - Methods used to implement background estimation are in hist_makers/hnl_maker_{BCKestimation}.py
    - Do not forget to specify histogram options in common/config/all/histograms/ folder

### Compute MainXsec File
- You can monitor main MC backgrounds composition in signal region (all MC and only true lepton)
```sh
python A_cutflow/produce_MCbackground_PieFig.py
```
- You can monitor anatuple cutflow for signal and background
```sh
python A_cutflow/produce_cutflow_fig.py 
```
- Then you need to produce Xsec file, needed for unc. computation
```sh
python A_cutflow/produce_XsecFile.py
```

### Compute FakeRate
- You can compute MC FFs using MC information only (“baseline method”) optional
```sh
python B_FakeRate/compute_MCFakeFactors.py
```

- You need to compute data derived FFs (used in the analysis):
    - You need to compute proportion of DY/ttbar in MC background in application region where leptons are not prompt
    ```sh
    python B_FakeRate/compute_weightsDYttbar.py
    ```
    - You need to compute correction factor for Light Lepton since we use ptcorr
    ```sh
    python B_FakeRate/compute_LightLeptonsCorrectionFactors.py
    ```
    - Compute Tau FFs in data
    ```sh
    python B_FakeRate/compute_TauFakeFactors.py
    ```
    - Compute Light Lepton FFs in data
    ```sh
    python B_FakeRate/compute_LightLeptonsFakeFactors.py
    ```
    - Then compute data derived FFs systematic uncertainty (evaluated in RegionUnc parameter) 
    ```sh
    python B_FakeRate/compute_FakeFactorsSystLL.py
    ```
    ```sh
    python B_FakeRate/compute_FakeFactorsSystTau.py
    ```

- You can plot the difference with between 2 sets of Fake rates, or look at FR variation 
```sh
python B_FakeRate/plot_diff_FR.py 
```
```sh
python B_FakeRate/plot_FRvar.py 
```

### Compute DNN score
load environment 
```sh
source C_DNNvar/env_MLenv.sh
```

- First you need to import model from GPUpc119 (where training was performed) in saved_models folder
- Then compute and save DNN score in json files
```sh
python C_DNNvar/save_DNN_score.py
```

### Produce TH1F hits (for limit and plot)//Plot CMS like histograms

#### General examples
Here some examples on how to:
- Produce TH1F hits
```sh
law run MakeTH1Hist --period 2018 --channel tem --tag AddJETcorr --BCKestimation FakeRate --PlotRegion SignalRegion 
```

- Produce CMS like plots
```sh
law run RunCMSPlot --period 2018 --channel tem --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```

#### For Fake rate monitoring
- Tau (!!!change tem config!!!)
    - ttbar
    FR perfect closure 
    ```sh
    law run RunCMSPlot --period 2018 --channel tem --tag AddJETcorr --BCKestimation ValidationTauFR --PlotRegion ttbarRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation
    ```sh
    law run RunCMSPlot --period 2018 --channel tll --tag AddJETcorr --BCKestimation ValidationTauFR --PlotRegion ttbarRegionValidation --branch 0 --UnblindData
    ```
    - DY
    FR perfect closure
    ```sh
    law run RunCMSPlot --period 2018 --channel tll --tag AddJETcorr --BCKestimation ValidationTauFR --PlotRegion DYRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation
    ```sh
    law run RunCMSPlot --period 2018 --channel tll --tag AddJETcorr --BCKestimation ValidationTauFR --PlotRegion DYRegionValidation --branch 0 --UnblindData
    ```

- Muon
    - ttbar
    FR perfect closure
    ```sh
    law run RunCMSPlot --period 2018 --channel llmu --tag LightLeptFFV2 --BCKestimation ValidationMuFR --PlotRegion ttbarRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation
    ```sh
    law run RunCMSPlot --period 2018 --channel llmu --tag LightLeptFFV2 --BCKestimation ValidationMuFR --PlotRegion ttbarRegionValidation --branch 0 --UnblindData
    ```
    - DY
    FR perfect closure
    ```sh
    law run RunCMSPlot --period 2018 --channel Zmu --tag LightLeptFFV2 --BCKestimation ValidationMuFR --PlotRegion DYRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation
    ```sh
    law run RunCMSPlot --period 2018 --channel Zmu --tag LightLeptFFV2 --BCKestimation ValidationMuFR --PlotRegion DYRegionValidation --branch 0 --UnblindData
    ```

- Electron
    - ttbar
    FR perfect closure
    ```sh
    law run RunCMSPlot --period 2018 --channel lle --tag LightLeptFFV2 --BCKestimation ValidationEleFR --PlotRegion ttbarRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation 
    ```sh
    law run RunCMSPlot --period 2018 --channel lle --tag LightLeptFFV2 --BCKestimation ValidationEleFR --PlotRegion ttbarRegionValidation --branch 0 --UnblindData
    ```
    - DY
    FR perfect closure
    ```sh
    law run RunCMSPlot --period 2018 --channel Ze --tag LightLeptFFV2 --BCKestimation ValidationEleFR --PlotRegion DYRegionFR --branch 0 --UnblindData
    ```
    syst unc estimation
    ```sh
    law run RunCMSPlot --period 2018 --channel Ze --tag LightLeptFFV2 --BCKestimation ValidationEleFR --PlotRegion DYRegionValidation --branch 0 --UnblindData
    ```

#### Check closure in InvertedBjetsVetoRegion
change channel config
- tem
```sh
law run RunCMSPlot --period 2018 --channel tem --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```
- tmm
```sh
law run RunCMSPlot --period 2018 --channel tmm --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```
- ttm
```sh
law run RunCMSPlot --period 2018 --channel ttm --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```
- tee
```sh
law run RunCMSPlot --period 2018 --channel tee --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```
- tte
```sh
law run RunCMSPlot --period 2018 --channel tte --tag AddJETcorr --BCKestimation FakeRate --PlotRegion InvertedBjetsVetoRegion --branch 8 --UnblindData
```

### Compute ExpLimit
For each NP, check wether we should use lnN or the shape of the NP (we commute Chi2 of the ratio of between nominal value and the NP, and if p-value high —> Shape)
```sh
python E_lnN_or_Shape/lnN_or_Shape.py
```

### go to limit code and compute the limits (link to github repo here)

### Limit plots
check the limits by discriminating to check which one is the best at each mass point
```sh
python F_Exp_limit/run_ExpLimitStudy.py
```
produce Brazil plot by channel
```sh
python Exp_limit/run_brazPlot.py
```
produce Brazil plot channel combined
```sh
python run_BDVbrazPlot_allchannels.py
```

## Additional useful scripts
### cleanup jobs folder
```sh
python other_scripts/cleanup_jobs_folder.py 
```
### combine tmm and tee channels
```sh
python other_scripts/make_tll_cha
```

## Documentation
- law documentation: https://luigi.readthedocs.io/en/stable/#
