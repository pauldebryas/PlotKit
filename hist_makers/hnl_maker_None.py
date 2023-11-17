import hist
from hist_makers.helpers import ToRootHist, compute_var_to_plot, load_ntuples, load_inputs, load_xbins
from hist_makers.regions.Baseline_regions import compute_region_mask
import numpy as np

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None):

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  input_files = load_inputs(inputs_cfg, input_dir)
  branches = load_ntuples(input_files)

  for process in branches.keys():
    cut_region = compute_region_mask(branches[process], channel) # no need MC info
    branch_to_plot = np.array(compute_var_to_plot(branches[process], hist_name)).flatten()
    hists[process] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    hists[process].fill(x=branch_to_plot[cut_region['ControlRegion_pass']], weight=branches[process]['genWeight'][cut_region['ControlRegion_pass']])

  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }
  return root_hists