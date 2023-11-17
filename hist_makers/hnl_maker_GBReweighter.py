import hist
import os
import numpy as np
import json
from hist_makers.helpers import ToRootHist, compute_var_to_plot, load_ntuples, load_inputs, load_xbins
from hist_makers.regions.GBReweighter_regions import compute_region_mask

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None):

  if channel not in ['tmm','tee']:
     raise('channel not implement yet')
    
  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  input_files = load_inputs(inputs_cfg, input_dir)
  branches = load_ntuples(input_files)

  # check if file exist first
  weights_CR_file = os.path.join(os.getenv("RUN_PATH"), f'GBReweighter/results/gb_weights_ControlRegion_{channel}.json')  
  with open(weights_CR_file, 'r') as file:
      weights_CR = json.load(file)

  #Extrapolated background
  cut_region = compute_region_mask(branches['data'], channel)
  Control_region_fail = np.array(compute_var_to_plot(branches['data'], hist_name)[cut_region['ControlRegion_fail']]).flatten()

  # check that we have the same selection
  if len(branches['data'][cut_region['ControlRegion_fail']]) != len(weights_CR['gb_weights']):
    raise 'error in regions definition (weight computation vs CR)'
  '''
  # check if we have the same events in the weights 
  extrapolated_region = branches['data'][cut_region['ControlRegion_fail']]
  for event_nb in extrapolated_region['event']:
    if event_nb not in weights_CR['event']:
      print(f'warning: {event_nb} not in weights_CR')
  '''
  sf = np.ones(len(weights_CR['gb_weights']))*weights_CR['sf']
  hists['ExtrapolatedBackground'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['ExtrapolatedBackground'].fill(x=Control_region_fail, weight=sf*weights_CR['gb_weights'])

  #signal
  cut_region = compute_region_mask(branches['HNL'], channel)
  hists['HNL'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['HNL'].fill(x=np.array(compute_var_to_plot(branches['HNL'], hist_name)).flatten()[cut_region['ControlRegion_pass']], weight=branches['HNL']['genWeight'][cut_region['ControlRegion_pass']])

  #data
  cut_region = compute_region_mask(branches['data'], channel)
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['data'].fill(x=np.array(compute_var_to_plot(branches['data'], hist_name)).flatten()[cut_region['ControlRegion_pass']], weight=branches['data']['genWeight'][cut_region['ControlRegion_pass']])

  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }

  return root_hists

