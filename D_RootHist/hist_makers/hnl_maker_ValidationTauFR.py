import hist
import numpy as np
import correctionlib.convert
import os

from common.helpers import load_ntuples
from D_RootHist.hist_makers.helpers import ToRootHist_val, compute_var_to_plot, load_inputs, load_xbins
from common.regions.regions import compute_region_mask

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None, period=None, PlotRegion = False, tag =None):

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  input_files = load_inputs(inputs_cfg, input_dir)

  if PlotRegion in ['ttbarRegionFR', 'ttbarRegionValidation']:
     FRregion = 'ttbarRegionFR'
  if PlotRegion in ['DYRegionFR', 'DYRegionValidation']:
     FRregion = 'DYRegionFR'

  # Load FR and err on FR
  FR = {}
  FR['Tau'] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactorsTau/{FRregion}/P_fake_TL_Tau.json'))

  branches = load_ntuples(input_files, 'Events;1')

  #True Tau MC background
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  print(f"Nevent in TrueLepton: {np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP_TauIsPromptLepton']])}\n")
  hists['TrueLepton'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['TrueLepton'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name)).flatten()[cut_region[f'{PlotRegion}_PassTightWP_TauIsPromptLepton']], weight=branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP_TauIsPromptLepton']])
  
  # FakeBackground (FR method)
  print(f'For FakeBackground ...')
  Fakes, Weights = apply_FR_tt(FR, branches, channel, hist_name, PlotRegion)
  hists['FakeBackground'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['FakeBackground'].fill(x= Fakes, weight=Weights)

  #data
  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  print(f"Nevent in data: {np.sum(branches['data']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])}\n")
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x').Weight()
  hists['data'].fill(x=np.array(compute_var_to_plot(branches['data'], hist_name)).flatten()[cut_region[f'{PlotRegion}_PassTightWP']], weight=branches['data']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])

  root_hists = { hist_name : ToRootHist_val(h) for hist_name, h in hists.items() }

  return root_hists


def apply_FR_tt(FR, branches, channel, hist_name, PlotRegion):

  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  sumw_MC = np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassLooseNotTightWP_TauIsPromptLepton']])
  print(f'sumw TrueHadronicTau events in AppRegion = {sumw_MC} \n')

  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  cut = cut_region[f'{PlotRegion}_PassLooseNotTightWP']
  N_events = np.sum(cut)
  print(f'N_events in AppRegion= {N_events} \n')

  weight = FR['Tau']['fake_rate'].evaluate(np.array(branches['data'][f'Tau_pt'][cut]).flatten(), np.array(np.abs(branches['data'][f'Tau_eta'][cut])).flatten())
  weights = weight/(1-weight)

  FakesProp = (N_events - sumw_MC)/N_events
  if FakesProp < 0:
      FakesProp = 0

  print(f'FakesProp = {round(FakesProp,1)} \n')
  Fakes = np.array(compute_var_to_plot(branches['data'], hist_name)[cut]).flatten()
  Weights = weights*FakesProp
  print(f'FFs * FakesProp = {round(np.sum(Weights)*FakesProp,1)} \n')

  return Fakes, Weights

