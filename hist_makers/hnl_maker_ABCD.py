import hist
import os
import numpy as np
from hist_makers.helpers import ToRootHist, compute_var_to_plot, load_ntuples, load_inputs, load_xbins
from hist_makers.regions.ABCD_regions import compute_region_cut

def compute_QCD(inputs, hist_name, channel):
  sumw_MCevents_regionA = 0
  sumw_MCevents_regionB = 0
  mc_regionC = []
  weights_regionC = []

  for process in inputs.keys():
    if process not in ['data', 'HNL']:
      tree = inputs[process]
      cut_region = compute_region_cut(tree, channel)
      sumw_MCevents_regionA += np.sum(tree["genWeight"][cut_region['A']])
      sumw_MCevents_regionB += np.sum(tree["genWeight"][cut_region['B']])
      mc_regionC += list(compute_var_to_plot(tree, hist_name)[cut_region['C']])
      weights_regionC += list( tree["genWeight"][cut_region['C']]* (-1))
  
  tree = inputs['data']
  cut_region = compute_region_cut(tree, channel)
  NA = np.sum(tree['genWeight'][cut_region['A']]) - sumw_MCevents_regionA
  NB = np.sum(tree['genWeight'][cut_region['B']]) - sumw_MCevents_regionB
  mc_regionC += list(compute_var_to_plot(tree, hist_name)[cut_region['C']])
  weights_regionC += list( tree['genWeight'][cut_region['C']])
  TF = (NB/NA)
  return mc_regionC, (np.array(weights_regionC))*TF

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None):

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  input_files = load_inputs(inputs_cfg, input_dir)
  branches = load_ntuples(input_files)

  QCD, weights_QCD = compute_QCD(branches, hist_name, channel)

  hists['QCD'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['QCD'].fill(x=QCD, weight=weights_QCD)

  for process in branches.keys():
    cut_region = compute_region_cut(branches[process], channel)
    branch_to_plot = compute_var_to_plot(branches[process], hist_name)
    hists[process] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    hists[process].fill(x=branch_to_plot[cut_region['D']], weight=branches[process]['genWeight'][cut_region['D']])

  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }
  return root_hists

