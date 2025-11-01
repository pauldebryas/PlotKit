import hist
import numpy as np
import correctionlib.convert
import os
import yaml

from common.helpers import load_ntuples
from D_RootHist.hist_makers.helpers import ToRootHist_val, compute_var_to_plot, load_inputs, load_xbins
from common.regions.regions import compute_region_mask

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None, period=None, PlotRegion = False, tag =None):

  if '_HNLMass' in hist_name:
      hist_name_split = hist_name.split('_HNLMass')
      hist_name = hist_name_split[0]
      MassHNL_Hyp = int(hist_name_split[-1])
  else:
      MassHNL_Hyp = '300'

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  is_flow = hist_cfg[hist_name].get('flow', False)
  input_files = load_inputs(inputs_cfg, input_dir)

  # load Corr factor from YAML file
  corrFactor_path = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/correctionFactorsLL.yml')
  with open(corrFactor_path, 'r') as file:
      CorrFactor = yaml.safe_load(file)

  branches = load_ntuples(input_files, 'Events;1')

  #MC background
  print(f'For Monte Carlo Background ...')
  for process in branches.keys():
    if process == 'data' or process.startswith('HNL'):
        continue
    print(f'    Processing {process} ...')
    cut_region = compute_region_mask(branches[process], channel, 'MC', PlotRegion)
    hists[process] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    hists[process].fill(x=np.array(compute_var_to_plot(branches[process], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region[f'{PlotRegion}_PassTightWP']], weight=branches[process]['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])
    print(f'')

  print(f'Computing Signal ...')
  cutregion = compute_region_mask(branches[f'HNL{MassHNL_Hyp}'], channel, 'MC', PlotRegion)
  cut = cutregion[f'{PlotRegion}_PassTightWP'] 
  Signal = np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut]
  Signal_w = branches[f'HNL{MassHNL_Hyp}']['genWeight'][cut]
  hists[f'HNL{MassHNL_Hyp}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists[f'HNL{MassHNL_Hyp}'].fill(x=Signal, weight=Signal_w)

  #data
  print(f'For Data ...')
  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['data'].fill(x=np.array(compute_var_to_plot(branches['data'], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region[f'{PlotRegion}_PassTightWP']], weight=branches['data']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])
  print(f'')
  
  root_hists = { hist_name : ToRootHist_val(h, flow= is_flow) for hist_name, h in hists.items() }

  return root_hists
