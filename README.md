# PlotKit

## How to install

Clone the repository:
```sh
git clone git@github.com:pauldebryas/PlotKit.git
```

The main conda environment needed to produce the plot:
```sh
conda env create -f environment_plot.yml
```

An other conda environment is needed to produce fake factors with BDT (using hep_ml package):
```sh
conda env create -f environment_MLHEP.yml
```

Don't forget to change the path to your miniconda in env_plot.sh and env_MLHEP.sh

This code is using law (luigi analysis workflow) to produce plots (see documentation) so first thing to do is to index the task:
```sh
source env_plot.sh
law index --verbose
```
If you decide to implement new task, you need to redo this step.

## How to run

### Compute FakeRate
- Load plot env:
```sh
source env_plot.sh
```
- Change Fake Rate computation parameters in run_fakeRate.py
- Run the code
    - Computing P_fake(T|L) as function of matching Jet pt/eta of the tau in the MCsample:
    ```sh
    python FakeRate/run_fakeRate.py
    ```
    - Computing P_fake(Sideband) as function of pt/eta of the tau in the MCsample:
    ```sh
    python FakeRate/run_FakePropInSideband.py
    ```
- Fake rates are stored in json files in 'FakeRate/results' folder and control plots in 'FakeRate/figures'

### Compute FakeFactors
- Load MLHEP env:
```sh
source env_MLHEP.sh
```
- Change Fake Factors computation parameters in run_weights_computation.py
- Run the code
    - Computing weights of BDT for background estimation:
    ```sh
    python GBReweighter/run_weights_computation.py
    ```
    - Plot additional figures for monitoring:
    ```sh
    python GBReweighter/run_plot.py
    ```
- Fake Factors are stored in json files in 'GBReweighter/results' folder and control plots in 'GBReweighter/figures'

### Plot CMS like histograms
- Load plot env:
```sh
source env_plot.sh
```

- Run the code on Condor (example):
```sh
law run RunPlot --period 2018 --channel tmm --tag AllHNLSamples --BCKestimation None --OutputTag None
```
    Here is the list of the parameters and a short description:
    - period: period of data taking: can be 2018, 2017, 2016 and 2016_HIPM
    - channel: channel you want to plot: can be ttm, tmm, tte, tee and tem
    - tag = tag used for anatuple production
    - BCKestimation = method used to estimate background. It can be:
        - None -> No data-driven background estimation (only MC)
        - ABCD -> Classical ABCD method
        - FakeRate -> using fake rate
        - GBReweighter -> using BDT
    - OutputTag: tag where to store the histograms (in output/{period}/{tag}_{OutputTag}/{channel})

- Option to add for debugging (run the code locally and only 1 plot)
```sh
--RunPlot-workflow local --branch 0
```

- Other useful informations:
    - All regions used to do background estimation are defined in hist_makers/regions/
    - Methods used to implement background estimation are in hist_makers/hnl_maker_{BCKestimation}.py
    - If you want to add more variable to plot, go in plotting/tasks.py and change 'self.var_to_plot' paramter in the corresponding channel
    - Do not forget to specify histogram options in config/hnl/hnl_{channel}/histograms.yaml
    - You can also tune what to plot in the histogram in config/hnl/hnl_{channel}/inputs_{BCKestimation}.yaml

## Documentation
- hep_ml package: https://arogozhnikov.github.io/hep_ml/index.html
- Akward array: https://awkward-array.readthedocs.io/en/latest/index.html
- law documentation: https://luigi.readthedocs.io/en/stable/#
