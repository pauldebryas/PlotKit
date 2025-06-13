import hist
import numpy as np
import correctionlib.convert
import os
import yaml

from common.helpers import load_ntuples
from D_RootHist.hist_makers.helpers import ToRootHist_val, compute_var_to_plot, load_inputs, load_xbins
from common.regions.regions import compute_region_mask

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None, period=None, PlotRegion = False, tag =None):

  name_lep = 'Muon'

  hists = {}
  x_bins = load_xbins(hist_cfg, hist_name)
  is_flow = hist_cfg[hist_name].get('flow', False)
  input_files = load_inputs(inputs_cfg, input_dir)

  if PlotRegion in ['ttbarRegionFR', 'ttbarRegionValidation']:
     FRregion = 'ttbarRegionFR'
  if PlotRegion in ['DYRegionFR', 'DYRegionValidation']:
     FRregion = 'DYRegionFR'
  
  # Load FR and err on FR
  FR = {}
  FR[name_lep] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactors{name_lep}/{FRregion}/P_fake_TL_{name_lep}.json'))
  # load Corr factor from YAML file
  corrFactor_path = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/correctionFactorsLL.yml')
  with open(corrFactor_path, 'r') as file:
      CorrFactor = yaml.safe_load(file)

  branches = load_ntuples(input_files, 'Events;1')

  #True MC background
  print(f'For Prompt Background ...')
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  print(f" - Prompt {name_lep} events in TightWP: {np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP_{name_lep}IsPromptLepton']])}\n")
  hists['TrueLepton'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['TrueLepton'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region[f'{PlotRegion}_PassTightWP_{name_lep}IsPromptLepton']], weight=branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP_{name_lep}IsPromptLepton']])
  print(f'')
  
  # FakeBackground (FR method)
  print(f'For FakeBackground ...')
  Fakes, Weights = apply_FR_1LL(FR, branches, channel, hist_name, PlotRegion, CorrFactor, name_lep)
  hists['FakeBackground'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['FakeBackground'].fill(x= Fakes, weight=Weights)
  print(f'')

  #data
  print(f'For Data ...')
  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  print(f" - Nevent in data: {np.sum(branches['data']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])}\n")
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['data'].fill(x=np.array(compute_var_to_plot(branches['data'], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region[f'{PlotRegion}_PassTightWP']], weight=branches['data']['genWeight'][cut_region[f'{PlotRegion}_PassTightWP']])
  print(f'')
  
  root_hists = { hist_name : ToRootHist_val(h, flow= is_flow) for hist_name, h in hists.items() }

  return root_hists

def apply_FR_1LL(FR, branches, channel, hist_name, PlotRegion, CorrFactor, name_lep):

  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  sumw_MC = np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassLooseNotTightWP_{name_lep}IsPromptLepton']])
  print(f' - Prompt {name_lep} events in AppRegion = {sumw_MC} \n')

  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  cut = cut_region[f'{PlotRegion}_PassLooseNotTightWP']
  N_events = np.sum(cut)
  print(f' - N_events in AppRegion= {N_events} \n')

  #ConeCorrectedPt = branches['data'][f'{name_lep}_pt']*(1+branches['data'][f'{name_lep}_pfRelIso03_all']-0.15)
  pt_parton = np.concatenate(branches['data'][f'{name_lep}_ConeCorrectedPt'])*CorrFactor[name_lep]
  ConeCorrectedPt = np.where(cut, pt_parton, branches['data'][f'{name_lep}_pt'])

  weight = FR[name_lep]['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{name_lep}_eta'][cut])).flatten())
  weights = weight/(1-weight)

  FakesProp = (N_events - sumw_MC)/N_events
  if FakesProp < 0:
      FakesProp = 0

  print(f' - FakesProp = {round(FakesProp,1)} \n')
  Fakes = np.array(compute_var_to_plot(branches['data'], hist_name, CorrFactor = CorrFactor)[cut]).flatten()
  Weights = weights*FakesProp
  print(f' - FR weights = {round(np.sum(Weights),1)} \n')

  return Fakes, Weights

'''
def compute_FRweight(FR, branches, Lepton, cut):
  #Lep_pt = branches['data'][f'{Lepton}_pt'][cut]
  ConeCorrectedPt = branches['data'][f'{Lepton}_pt']*(1+branches['data'][f'{Lepton}_pfRelIso03_all']-0.15)
  #Lep_pt*(1+branches['data'][f'{Lepton}_pfRelIso03_all'][cut])
  weights = FR[name_lep]['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{Lepton}_eta'][cut])).flatten())
  return weights

def compute_FR(weight):
  return weight/(1-weight)

def apply_FR_2Mu(FR, branches, channel, hist_name, PlotRegion):
  application_region = ['PFF','PFP','PPF']

  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  sumw_MC = {}
  N_MC = {}
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  #Set values to 0 ...
  for XXX in application_region:
    sumw_MC[f'{PlotRegion}{XXX}_TrueLepton'] = 0
    N_MC[f'{PlotRegion}{XXX}_TrueLepton'] = 0
  # ... and fill values
  for XXX in application_region:
    sumw_MC[f'{PlotRegion}{XXX}_TrueLepton'] += np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}{XXX}_TrueLepton']])
    N_MC[f'{PlotRegion}{XXX}_TrueLepton'] += np.sum(cut_region[f'{PlotRegion}{XXX}_TrueLepton'])

  #Save events in App region and corresponding FR
  Fakes = []
  Weights = []
  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  for XXX in application_region:
    cut = cut_region[f'{PlotRegion}{XXX}']
    print(f'{XXX} region: \n')
    N_events_XXX = np.sum(cut)
    print(f'N_events AppRegion= {N_events_XXX} \n')

    if XXX == 'PFF':
      weights2 = compute_FRweight(FR, branches, 'Muon1', cut)
      weights3 = compute_FRweight(FR, branches, 'Muon2', cut)
      weights = (-1)*compute_FR(weights2)*compute_FR(weights3)

    if XXX == 'PFP':
      weights2 = compute_FRweight(FR, branches, 'Muon1', cut)
      weights = compute_FR(weights2)

    if XXX == 'PPF':
      weights3 = compute_FRweight(FR, branches, 'Muon2', cut)
      weights = compute_FR(weights3)

    print(f'sumw FFs= {round(np.sum(weights),1)} \n')

    FakesProp = (N_events_XXX - sumw_MC[f'{PlotRegion}{XXX}_TrueLepton'])/N_events_XXX
    if FakesProp < 0:
      FakesProp = 0

    print(f'sumw FakesProp = {round(FakesProp,1)} \n')

    Fakes = np.concatenate( [Fakes, np.array(compute_var_to_plot(branches['data'], hist_name)[cut]).flatten()])

    Weights = np.concatenate( [Weights, weights*FakesProp])

    print(f'FFs * FakesProp = {round(np.sum(weights)*FakesProp,1)} \n')
    print(f'\n')

  return Fakes, Weights

'''