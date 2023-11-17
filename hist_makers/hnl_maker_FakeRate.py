import hist
import os
import correctionlib.convert
import numpy as np
from hist_makers.helpers import ToRootHist, compute_var_to_plot, load_ntuples, load_inputs, load_xbins
from hist_makers.regions.FakeRate_regions import compute_region_mask

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None):

  MCSample = 'DY' #MCSample for FR computation
  if channel not in ['tmm','tee', 'tem']:
     raise('channel not implement yet')

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  input_files = load_inputs(inputs_cfg, input_dir)
  branches = load_ntuples(input_files)

  # check if file exist first
  fake_rate_file =  os.path.join(os.getenv("RUN_PATH"), f'FakeRate/results/P_fake_TL_{MCSample}.json')
  #err_fake_rate_file = os.path.join(os.getenv("RUN_PATH"), f'FakeRate/results/err_P_fake_TL_{MCSample}.json')
  FakeInSideband_file = os.path.join(os.getenv("RUN_PATH"), f'FakeRate/results/P_FakeInSideband_{channel}.json')
  #err_FakeInSideband_file =  os.path.join(os.getenv("RUN_PATH"), f'FakeRate/results/P_FakeInSideband_err_{channel}.json')

  fake_rate = correctionlib.CorrectionSet.from_file(fake_rate_file)
  #err_sf = correctionlib.CorrectionSet.from_file(err_fake_rate_file)

  # no need to do this step every time: check if json file exist
  if channel in ['tmm','tee','tem']:
    FakePropInSideband = correctionlib.CorrectionSet.from_file(FakeInSideband_file)
    #err_FakePropInSideband = correctionlib.CorrectionSet.from_file(err_FakeInSideband_file)

  #True Tau MC background
  for process in branches.keys():
    if process not in ['HNL', 'data']:
       cut_region = compute_region_mask(branches[process], channel, 'MC')
       hists[process] = hist.Hist.new.Variable(x_bins, name='x').Weight()
       hists[process].fill(x=np.array(compute_var_to_plot(branches[process], hist_name)).flatten()[cut_region['ControlRegion_pass_true']], weight=branches[process]['genWeight'][cut_region['ControlRegion_pass_true']])

  #Fake Tau background
  if channel in ['tmm','tee','tem']:
    cut_region = compute_region_mask(branches['data'], channel, 'data')
    FakeRate = np.array(compute_var_to_plot(branches['data'], hist_name)[cut_region['ControlRegion_fail']]).flatten()
    weights_fakeRate = fake_rate['fake_rate'].evaluate(np.array(branches['data']['Tau_Jet_pt'][cut_region['ControlRegion_fail']]).flatten(), np.array(np.abs(branches['data']['Tau_Jet_eta'][cut_region['ControlRegion_fail']])).flatten())
    weights_FakePropInSideband = FakePropInSideband['fake_rate'].evaluate(np.array(branches['data']['Tau_pt'][cut_region['ControlRegion_fail']]), np.array(np.abs(branches['data']['Tau_eta'][cut_region['ControlRegion_fail']])))
    weights_fakeRate = weights_fakeRate/(1-weights_fakeRate)
    weights = weights_fakeRate*weights_FakePropInSideband
    hists['FakeTauBackground'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    hists['FakeTauBackground'].fill(x=FakeRate, weight=weights)
  
  #signal
  cut_region = compute_region_mask(branches['HNL'], channel, 'MC')
  hists['HNL'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['HNL'].fill(x=np.array(compute_var_to_plot(branches['HNL'], hist_name)).flatten()[cut_region['ControlRegion_pass']], weight=branches['HNL']['genWeight'][cut_region['ControlRegion_pass']])

  #data
  cut_region = compute_region_mask(branches['data'], channel, 'data')
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['data'].fill(x=np.array(compute_var_to_plot(branches['data'], hist_name)).flatten()[cut_region['ControlRegion_pass']], weight=branches['data']['genWeight'][cut_region['ControlRegion_pass']])

  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }

  return root_hists

'''
  if channel in ['ttm','tte']:
    cut_region = compute_region_mask(branches['data'], channel, 'data')
    FakeRate = np.array(compute_var_to_plot(branches['data'], hist_name)[cut_region['SidebandPF']])
    weights_fakeRate = fake_rate['fake_rate'].evaluate(np.array(branches['data']['Tau2_pt'][cut_region['SidebandPF']]), np.array(np.abs(branches['data']['Tau2_eta'][cut_region['SidebandPF']])))
    #weights_fakeRate = fake_rate.to_evaluator().evaluate(np.array(branches['data']['Tau2_pt'][cut_region['SidebandPF']]), np.array(np.abs(branches['data']['Tau2_eta'][cut_region['SidebandPF']])))
    weights_fakeRate = weights_fakeRate/(1-weights_fakeRate)
    hists['SingleFakeTauBackground'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    hists['SingleFakeTauBackground'].fill(x=FakeRate, weight=weights_fakeRate)

    FakeRate = np.array(compute_var_to_plot(branches['data'], hist_name)[cut_region['SidebandFF']])
    weights_fakeRate_Tau1 = fake_rate['fake_rate'].evaluate(np.array(branches['data']['Tau1_pt'][cut_region['SidebandFF']]), np.array(np.abs(branches['data']['Tau1_eta'][cut_region['SidebandFF']])))
    weights_fakeRate_Tau2 = fake_rate['fake_rate'].evaluate(np.array(branches['data']['Tau2_pt'][cut_region['SidebandFF']]), np.array(np.abs(branches['data']['Tau2_eta'][cut_region['SidebandFF']])))
    #weights_fakeRate_Tau1 = fake_rate.to_evaluator().evaluate(np.array(branches['data']['Tau1_pt'][cut_region['SidebandFF']]), np.array(np.abs(branches['data']['Tau1_eta'][cut_region['SidebandFF']])))
    #weights_fakeRate_Tau2 = fake_rate.to_evaluator().evaluate(np.array(branches['data']['Tau2_pt'][cut_region['SidebandFF']]), np.array(np.abs(branches['data']['Tau2_eta'][cut_region['SidebandFF']])))
    weights_fakeRate_Tau1 = weights_fakeRate_Tau1/(1-weights_fakeRate_Tau1)
    weights_fakeRate_Tau2 = weights_fakeRate_Tau2/(1-weights_fakeRate_Tau2)
    hists['DoubleFakeTauBackground'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
    hists['DoubleFakeTauBackground'].fill(x=FakeRate, weight=weights_fakeRate_Tau1*weights_fakeRate_Tau2)
'''