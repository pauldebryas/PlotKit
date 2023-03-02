# PlotKit

# How to install

```sh
git clone git@github.com:kandrosov/PlotKit.git
```

Example of a minimal conda environment setup:
```sh
conda create -n plot python=3.9
conda activate plot
conda install boost hist ROOT uproot
```

# How to run

Example how to run from a command line:
```sh
python Plotter.py --page-cfg config/examples/hnl_tau/cms_stacked.yaml --page-cfg-custom config/examples/hnl_tau/2018.yaml --hist-cfg config/examples/hnl_tau/histograms.yaml --inputs-cfg config/examples/hnl_tau/inputs.yaml --hist-name Tau1_pt --custom "cat_text=#tau_{h}#tau_{h}#mu" --output output/Tau1_pt.pdf --hist-maker hist_makers/hnl_tau.py output/ttm_DeepTau2p1
```

Example how to import Plotter in python:
```python
from PlotKit.Plotter import Plotter
plotter = Plotter(page_cfg=page_cfg, page_cfg_custom=page_cfg_custom, hist_cfg=hist_cfg, inputs_cfg=inputs_cfg)
plotter.plot(hist_name1, hists1, output1, custom=custom1)
plotter.plot(hist_name2, hists2, output2, custom=custom2)
...
```
