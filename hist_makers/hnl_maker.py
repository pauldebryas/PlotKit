import hist
import os
import uproot
import ROOT
import numpy as np
from hist_makers.helpers import ToRootHist, compute_var_to_plot
from config.hnl.ABCD_regions import compute_region_cut

def compute_QCD(inputs, input_dir, hist_name, channel):
  mc_regionC = []
  weights_regionC = []
  sumw_MCevents_regionA = 0
  sumw_MCevents_regionB = 0
  MC_sample_files = []
  for MC_process in inputs.keys():
    if MC_process not in ['data', 'HNL']:
      for file in inputs[MC_process]:
        MC_sample_files.append(file)

  for file_name in MC_sample_files:
    file = uproot.open(os.path.join(input_dir, file_name))
    tree = file['Event']
    cut_region = compute_region_cut(tree, channel)
    sumw_MCevents_regionA += np.sum(tree["genWeight"].array()[cut_region['A']])
    sumw_MCevents_regionB += np.sum(tree["genWeight"].array()[cut_region['B']])
    mc_regionC += list(compute_var_to_plot(tree, hist_name)[cut_region['C']])
    weights_regionC += list( tree["genWeight"].array()[cut_region['C']]* (-1))
  
  data_file = inputs['data'][0]
  file = uproot.open(os.path.join(input_dir, data_file))
  tree = file['Event']
  cut_region = compute_region_cut(tree, channel)
  NA = np.sum(tree['genWeight'].array()[cut_region['A']]) - sumw_MCevents_regionA
  NB = np.sum(tree['genWeight'].array()[cut_region['B']]) - sumw_MCevents_regionB
  mc_regionC += list(compute_var_to_plot(tree, hist_name)[cut_region['C']])
  weights_regionC += list( tree['genWeight'].array()[cut_region['C']])

  TF = (NB/NA)
  return mc_regionC, (np.array(weights_regionC))*TF

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None):

  input_files = {}
  for input in inputs_cfg:
    if 'files' in input.keys():
      files_list = [elem for elem in input['files']] 
      for file in input['files']:
        if os.path.isfile(os.path.join(input_dir, file)) == False:
          print('WARNING: ' + file + ' is missing')
          files_list.remove(file)
      input_files[input['name']] = files_list

  hists = {}
  hist_desc = hist_cfg[hist_name]

  QCD, weights_QCD = compute_QCD(input_files, input_dir, hist_name, channel)
  x_bins = hist_desc['x_bins']
  if isinstance(x_bins, list):
    hists['QCD'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  else:
    n_bins, bin_range = x_bins.split('|')
    start,stop = bin_range.split(':')
    hists['QCD'] = hist.Hist.new.Regular(int(n_bins), float(start), float(stop), name='x').Weight()
  hists['QCD'].fill(x=QCD, weight=weights_QCD)

  for input_name, input_files_lst in input_files.items():
    x_bins = hist_desc['x_bins']
    if isinstance(x_bins, list):
      hists[input_name] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    else:
      n_bins, bin_range = x_bins.split('|')
      start,stop = bin_range.split(':')
      hists[input_name] = hist.Hist.new.Regular(int(n_bins), float(start), float(stop), name='x').Weight()

    for file_name in input_files_lst:
      file = uproot.open(os.path.join(input_dir, file_name))
      tree = file['Event']
      cut_region = compute_region_cut(tree, channel)
      branch_to_plot = compute_var_to_plot(tree, hist_name)
      hists[input_name].fill(x=branch_to_plot[cut_region['D']], weight=tree['genWeight'].array()[cut_region['D']])
  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }
  return root_hists
