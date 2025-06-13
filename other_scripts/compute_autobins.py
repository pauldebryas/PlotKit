import os
import correctionlib.convert
import numpy as np
import awkward as ak
import yaml
import json

from common.helpers import load_ntuples, get_hnl_masses
from D_RootHist.hist_makers.helpers import compute_var_to_plot, load_inputs, load_xbins
from common.regions.regions import compute_region_mask

# ------------------------- Parameters -------------------------
period= '2018'
tag = 'AddJETcorr'
channel = 'tmm'
# --------------------------------------------------------------
def compute_adaptive_binning(inputs, x_range, nbin_max, w_threshold_bck, max_relative_bck_err):

    min_val, max_val = x_range

    mask = (inputs['signal'] >= min_val) & (inputs['signal'] <= max_val)
    if len(inputs['signal'][mask]) < nbin_max:
      if len(inputs['signal'][mask]) == 0:
        print("Warning: No signal values in the specified range. Using a single bin.")
        nbin_max = 1
        n_threshold_signal = len(inputs['signal']) + 1
      else:
        nbin_max = len(inputs['signal'][mask])
        n_threshold_signal = np.floor(len(inputs['signal'][mask]) / nbin_max).astype(int)
    else:
      n_threshold_signal = np.floor(len(inputs['signal'][mask]) / nbin_max).astype(int)


    # Granularity control
    step_int_scale = 1e4  # granularity (0.0001 bins)
    min_value_int = int(min_val * step_int_scale)
    max_value_int = int(max_val * step_int_scale)
    n_bins = max_value_int - min_value_int + 1

    # Convert values to int bin positions
    signal_vals_int = (inputs['signal'] * step_int_scale).astype(int)
    bkg1_vals_int = (inputs['FakeBackground'] * step_int_scale).astype(int)
    bkg2_vals_int = (inputs['TrueLepton'] * step_int_scale).astype(int)

    # Initialize histograms
    sig_hist = np.zeros(n_bins)
    bkg1_hist = np.zeros(n_bins)
    bkg2_hist = np.zeros(n_bins)
    bkg1_w2_hist = np.zeros(n_bins)
    bkg2_w2_hist = np.zeros(n_bins)

    # Fill signal histogram (counts only)
    for v in signal_vals_int:
        idx = v - min_value_int
        if 0 <= idx < n_bins:
            sig_hist[idx] += 1

    # Fill background1 histogram (weighted yield and sum of weights squared)
    for v, w in zip(bkg1_vals_int, inputs['FakeBackground_w']):
        idx = v - min_value_int
        if 0 <= idx < n_bins:
            bkg1_hist[idx] += w
            bkg1_w2_hist[idx] += w**2

    # Fill background2 histogram (weighted yield and sum of weights squared)
    for v, w in zip(bkg2_vals_int, inputs['TrueLepton_w']):
        idx = v - min_value_int
        if 0 <= idx < n_bins:
            bkg2_hist[idx] += w
            bkg2_w2_hist[idx] += w**2

    # Compute bin edges from top to bottom
    edges = [max_value_int]
    i = max_value_int
    while i > min_value_int:
        acc_sig = 0
        acc_bkg1 = 0
        acc_bkg2 = 0
        acc_bkg1_w2 = 0
        acc_bkg2_w2 = 0

        j = i
        while j >= min_value_int:
            idx = j - min_value_int
            acc_sig += sig_hist[idx]
            acc_bkg1 += bkg1_hist[idx]
            acc_bkg2 += bkg2_hist[idx]
            acc_bkg1_w2 += bkg1_w2_hist[idx]
            acc_bkg2_w2 += bkg2_w2_hist[idx]

            # Relative errors (avoid zero division)
            bkg1_err_ok = acc_bkg1 > 0 and np.sqrt(acc_bkg1_w2) / acc_bkg1 <= max_relative_bck_err
            bkg2_err_ok = acc_bkg2 > 0 and np.sqrt(acc_bkg2_w2) / acc_bkg2 <= max_relative_bck_err
            total_bkg_yield = acc_bkg1 + acc_bkg2

            j -= 1
            if acc_sig >= n_threshold_signal and total_bkg_yield >= w_threshold_bck and bkg1_err_ok and bkg2_err_ok:
                break
        # print('') 
        # print(f'Sig yield: {acc_sig}')
        # print(f'bkg yield: {total_bkg_yield}')
        # print(f'True bkg yield: {acc_bkg1}')
        # print(f'Fake bkg yield: {acc_bkg2}')
        # print(f'Rel. err. True bkg: {np.sqrt(acc_bkg1_w2) / acc_bkg1}')
        # print(f'Rel. err. Fake bkg: {np.sqrt(acc_bkg2_w2) / acc_bkg2}')
        # print('')
        edges.append(j)
        if j <= min_value_int:
            if (acc_sig >= n_threshold_signal*0.5 and total_bkg_yield >= w_threshold_bck and bkg1_err_ok and bkg2_err_ok) | (nbin_max == 1):
              break
            else:
              print('removing last edge')
              edges.remove(j)
              break
        i = j

    # Finalize edges
    edges = sorted(set(edges))
    edges = np.array(edges, dtype=float) / step_int_scale

    # Ensure full coverage
    if edges[0] > min_val:
        edges[0] = min_val
    if edges[-1] < max_val:
        edges[-1] = max_val

    x_bins = edges.tolist()
    print(f"Adaptive binning complete: {len(x_bins) - 1} bins.")
    return x_bins 

def load_inputs_cfg(BCKestimationMethod, period, channel):
    
    inputs = []
    config_path_BCKG = os.path.join(os.getenv("RUN_PATH"),"common","config", "all", "inputs")

    #channel config files
    if BCKestimationMethod in ['FakeRate','ValidationMuFR', 'ValidationEleFR', 'ValidationTauFR']:
        #load TrueLepton background
        with open(os.path.join(config_path_BCKG, 'inputs_TrueLepton.yaml'), 'r') as f:
            MC_inputs = yaml.safe_load(f)
            if BCKestimationMethod in ['ValidationMuFR', 'ValidationEleFR']:
                MC_inputs[0]['files'].remove('ZHToTauTau_anatuple.root')
        inputs.append(MC_inputs[0])

        #load Fake Background
        with open(os.path.join(config_path_BCKG, 'inputs_FakeBackground.yaml'), 'r') as f:
            FB_inputs = yaml.safe_load(f)
        inputs.append(FB_inputs[0])

        #load signal
        if BCKestimationMethod == 'FakeRate':
            with open(os.path.join(config_path_BCKG, 'inputs_AllSignal.yaml'), 'r') as f:
                signal_inputs = yaml.safe_load(f)
            for i in range(len(signal_inputs)):
                inputs.append(signal_inputs[i])

        #load data
        with open(f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_{channel}/inputs_data.yaml', 'r') as f:
            data_inputs = yaml.safe_load(f)
        inputs.append(data_inputs[0])
    
    if len(inputs) == 0:
        raise (f'{BCKestimationMethod} not inplemented in load_inputs')
    
    return inputs

def AddDNNscore(branches, channel, MassHNL_Hyp, period, ModelName):
  # load score:
  json_file = f'/eos/user/p/pdebryas/HNL/DNNscore/{period}/{ModelName}/DNNscore_M{MassHNL_Hyp}.json'
  with open(json_file, 'r') as file:
      DNNscore = json.load(file)
  sample_list = DNNscore[channel].keys()
  #print(sample_list)
  for process in list(branches.keys()):
    #print(process)
    process_score =  (branches[process]['event']/branches[process]['event'])*-1. #initialize with default value
    for sample in sample_list:
      if f'mask_{sample}' in branches[process].fields:
        #print(f'   {sample}')
        mask_sample = branches[process][f'mask_{sample}'] == 1
        mask_events =  ak.where(np.isin(branches[process]['event'], DNNscore[channel][sample]['event']), True, False)
        mask_total = (mask_sample & mask_events)
        # Create a dictionary mapping sel_events to scores
        event_to_score = dict(zip(DNNscore[channel][sample]['event'], DNNscore[channel][sample]['scores']))
        where_score = np.vectorize(lambda x: event_to_score.get(x, -1.))(branches[process]['event'])
        # Ensure where_score has a proper numeric dtype
        where_score = np.array(where_score, dtype=np.float32)  # Convert to float32
        where_score = ak.Array(where_score)  # Convert to Awkward Array
        process_score = ak.where(mask_total, ak.Array(where_score) , process_score)
    branches[process] = ak.with_field(branches[process], process_score, "DNNscore")
    #print(branches[process]["DNNscore"])
  return branches

def load_FR(tag, period, tag_LL):
  Leptons = ['Electron', 'Muon', 'Tau']
  pathtag = {
    'Electron': tag_LL, 
    'Muon': tag_LL, 
    'Tau': tag
  }
  # Load FR and err on FR
  FR = {}
  FR_err = {}
  for lepton in Leptons:
    FR[f'{lepton}_ttbar'] = correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{pathtag[lepton]}/{period}/FakeFactors{lepton}/ttbarRegionFR/P_fake_TL_{lepton}.json'))
    FR_err[f'{lepton}_ttbar'] = correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{pathtag[lepton]}/{period}/FakeFactors{lepton}/ttbarRegionFR/P_fake_TL_{lepton}_err.json'))
    FR[f'{lepton}_DY'] =    correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{pathtag[lepton]}/{period}/FakeFactors{lepton}/DYRegionFR/P_fake_TL_{lepton}.json'))
    FR_err[f'{lepton}_DY'] =    correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{pathtag[lepton]}/{period}/FakeFactors{lepton}/DYRegionFR/P_fake_TL_{lepton}_err.json'))
  return FR, FR_err

def load_W_DY_ttbar(tag, period, PlotRegion, channel):
  # Load DY and ttbar weights
  with open(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/Prop_DY_ttbar_in{PlotRegion}AppRegion.yml'), 'r') as stream: 
      file_content = yaml.safe_load(stream)
      W_DY_ttbar = file_content[channel]
  return W_DY_ttbar

def load_CorrFactor(tag_LL, period):
  # load Corr factor from YAML file
  with open(os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag_LL}/{period}/correctionFactorsLL.yml'), 'r') as file:
      CorrFactor = yaml.safe_load(file)
  return CorrFactor

def get_var_config(channel):
  with open(os.path.join(os.getenv("RUN_PATH"),f'common/config/all/histograms/histograms_{channel}.yaml'), 'r') as file:
      vars = yaml.safe_load(file)
  return vars

def compute_FRweight(FR, FR_err, branches, Lepton, cut, W_DY_ttbar, CorrFactor):
  Lepton_pt = np.array(branches['data'][f'{Lepton["name"]}_pt'][cut]).flatten()
  Lepton_abseta = np.array(np.abs(branches['data'][f'{Lepton["name"]}_eta'][cut])).flatten()
  if Lepton["type"] != 'Tau':
    Lepton_parton_pt = np.array(branches['data'][f'{Lepton["name"]}_ConeCorrectedPt'][cut]).flatten()*CorrFactor[Lepton["type"]]
    Lepton_iso = np.array(branches['data'][f'{Lepton["name"]}_pfRelIso03_all'][cut]).flatten()
    Lepton_pt = np.where(Lepton_iso < 0.15, Lepton_pt, Lepton_parton_pt)

  weights_DY = FR[f'{Lepton["type"]}_DY']['fake_rate'].evaluate(Lepton_pt, Lepton_abseta)
  weights_DY_err = FR_err[f'{Lepton["type"]}_DY']['fake_rate_err'].evaluate(Lepton_pt, Lepton_abseta)
  weights_ttbar = FR[f'{Lepton["type"]}_ttbar']['fake_rate'].evaluate(Lepton_pt, Lepton_abseta)
  weights_ttbar_err = FR_err[f'{Lepton["type"]}_ttbar']['fake_rate_err'].evaluate(Lepton_pt, Lepton_abseta)
  weights = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar
  weights_err = np.sqrt((W_DY_ttbar['DY'][0] * weights_DY_err) ** 2 + (W_DY_ttbar['ttbar'][0]* weights_ttbar_err) ** 2 + (weights_DY * W_DY_ttbar['DY'][1]) ** 2 + (weights_ttbar * W_DY_ttbar['ttbar'][1]) ** 2)

  return weights, weights_err

def compute_FR(weight):
  return weight/(1-weight)

def compute_ptcorr(branches, lepton, CorrFactor):
  lepton_pt = np.array(branches['data'][f'{lepton["name"]}_pt']).flatten()
  if lepton['type'] != 'Tau':
    lepton_parton_pt = np.array(branches['data'][f'{lepton["name"]}_ConeCorrectedPt']).flatten()*CorrFactor[lepton['type']]
    lepton_iso = np.array(branches['data'][f'{lepton["name"]}_pfRelIso03_all']).flatten()
    lepton_pt = np.where(lepton_iso < 0.15, lepton_pt, lepton_parton_pt)
  return lepton_pt

def apply_FR_method(FR, FR_err, branches, channel, hist_name, cutRegion, PlotRegion, W_DY_ttbar, CorrFactor, lepton, application_region):
  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  #Set values to 0 ...
  sumw_MC = {}
  N_MC = {}
  for XXX in application_region:
    sumw_MC[f'AppRegion{XXX}_TrueLepton'] = 0
    N_MC[f'AppRegion{XXX}_TrueLepton'] = 0
  # ... and fill values
  for XXX in application_region:
    cut = cutRegion['TrueLepton'][f'{PlotRegion}_AppRegion{XXX}_TrueLepton']
    sumw_MC[f'AppRegion{XXX}_TrueLepton'] += np.sum(branches['TrueLepton']['genWeight'][cut])
    N_MC[f'AppRegion{XXX}_TrueLepton'] += np.sum(cut)

  Fakes = {}
  Fakes['var'] = []
  Fakes['lepton1_ptcorr'] = []
  Fakes['lepton2_ptcorr'] = []
  Fakes['lepton3_ptcorr'] = []
  Weights = {}
  Weights['nom'] = []
  Weights['l1_errUp'] = []
  Weights['l2_errUp'] = []
  Weights['l3_errUp'] = []
  Weights['l1_errDown'] = []
  Weights['l2_errDown'] = []
  Weights['l3_errDown'] = []
  Weights['FakeProp_errUp'] = {}
  Weights['FakeProp_errDown'] = {}
  for XXX in application_region:
    Weights['FakeProp_errUp'][XXX] = []
    Weights['FakeProp_errDown'][XXX] = []

  for XXX in application_region:
    cut = cutRegion['data'][f'{PlotRegion}_AppRegion{XXX}']
    N_events_XXX = np.sum(cut)

    if XXX == 'FFF':
      weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
      weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
      weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)
      weights = compute_FR(weights1)*compute_FR(weights2)*compute_FR(weights3)
      weights1_errUp =  compute_FR(weights1 + weights1_err)*compute_FR(weights2)*compute_FR(weights3)
      weights2_errUp =  compute_FR(weights1)*compute_FR(weights2 + weights2_err)*compute_FR(weights3)
      weights3_errUp =  compute_FR(weights1)*compute_FR(weights2)*compute_FR(weights3 + weights3_err)
      weights1_errDown =  compute_FR(weights1 - weights1_err)*compute_FR(weights2)*compute_FR(weights3)
      weights2_errDown =  compute_FR(weights1)*compute_FR(weights2 - weights2_err)*compute_FR(weights3)
      weights3_errDown =  compute_FR(weights1)*compute_FR(weights2)*compute_FR(weights3 - weights3_err)

    if XXX == 'FFP':
      weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
      weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
      weights = (-1)*compute_FR(weights1)*compute_FR(weights2)
      weights1_errUp =  (-1)*compute_FR(weights1 + weights1_err)*compute_FR(weights2)
      weights2_errUp =  (-1)*compute_FR(weights1)*compute_FR(weights2 + weights2_err)
      weights3_errUp =  (-1)*compute_FR(weights1)*compute_FR(weights2)
      weights1_errDown =  (-1)*compute_FR(weights1 - weights1_err)*compute_FR(weights2)
      weights2_errDown =  (-1)*compute_FR(weights1)*compute_FR(weights2 - weights2_err)
      weights3_errDown =  (-1)*compute_FR(weights1)*compute_FR(weights2)

    if XXX == 'FPF':
      weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
      weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)
      weights = (-1)*compute_FR(weights1)*compute_FR(weights3)
      weights1_errUp =  (-1)*compute_FR(weights1 + weights1_err)*compute_FR(weights3)
      weights2_errUp = (-1)*compute_FR(weights1)*compute_FR(weights3)
      weights3_errUp =  (-1)*compute_FR(weights1)*compute_FR(weights3 + weights3_err)
      weights1_errDown =  (-1)*compute_FR(weights1 - weights1_err)*compute_FR(weights3)
      weights2_errDown = (-1)*compute_FR(weights1)*compute_FR(weights3)
      weights3_errDown =  (-1)*compute_FR(weights1)*compute_FR(weights3 - weights3_err)

    if XXX == 'PFF':
      weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
      weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)
      weights = (-1)*compute_FR(weights2)*compute_FR(weights3)
      weights1_errUp =  (-1)*compute_FR(weights2)*compute_FR(weights3)
      weights2_errUp =  (-1)*compute_FR(weights2 + weights2_err)*compute_FR(weights3)
      weights3_errUp =  (-1)*compute_FR(weights2)*compute_FR(weights3 + weights3_err)
      weights1_errDown =  (-1)*compute_FR(weights2)*compute_FR(weights3)
      weights2_errDown =  (-1)*compute_FR(weights2 - weights2_err)*compute_FR(weights3)
      weights3_errDown =  (-1)*compute_FR(weights2)*compute_FR(weights3 - weights3_err)

    if XXX == 'FPP':
      weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
      weights = compute_FR(weights1)
      weights1_errUp =  compute_FR(weights1 + weights1_err)
      weights2_errUp =  compute_FR(weights1)
      weights3_errUp =  compute_FR(weights1)
      weights1_errDown =  compute_FR(weights1 - weights1_err)
      weights2_errDown =  compute_FR(weights1)
      weights3_errDown =  compute_FR(weights1)

    if XXX == 'PFP':
      weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
      weights = compute_FR(weights2)
      weights1_errUp =  compute_FR(weights2)
      weights2_errUp =  compute_FR(weights2 + weights2_err)
      weights3_errUp =  compute_FR(weights2)
      weights1_errDown =  compute_FR(weights2)
      weights2_errDown =  compute_FR(weights2 - weights2_err)
      weights3_errDown =  compute_FR(weights2)

    if XXX == 'PPF':
      weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)
      weights = compute_FR(weights3)
      weights1_errUp =  compute_FR(weights3)
      weights2_errUp =  compute_FR(weights3)
      weights3_errUp =  compute_FR(weights3 + weights3_err)
      weights1_errDown =  compute_FR(weights3)
      weights2_errDown =  compute_FR(weights3)
      weights3_errDown =  compute_FR(weights3 - weights3_err)

    FakesProp = (N_events_XXX - sumw_MC[f'AppRegion{XXX}_TrueLepton'])/N_events_XXX
    if FakesProp < 0:
      FakesProp = 0
    FakesProp_err = np.sqrt(N_events_XXX)/N_events_XXX

    Fakes['var'] = np.concatenate( [Fakes['var'], np.array(compute_var_to_plot(branches['data'], hist_name)[cut]).flatten()])
    #np.array(branches['data'][f"{lepton[channel]['lepton1']['name']}_pt"][cut]).flatten()
    Fakes['lepton1_ptcorr'] = np.concatenate( [Fakes['lepton1_ptcorr'], compute_ptcorr(branches, lepton[channel]['lepton1'], CorrFactor)[cut]])
    Fakes['lepton2_ptcorr'] = np.concatenate( [Fakes['lepton2_ptcorr'], compute_ptcorr(branches, lepton[channel]['lepton2'], CorrFactor)[cut]])
    Fakes['lepton3_ptcorr'] = np.concatenate( [Fakes['lepton3_ptcorr'], compute_ptcorr(branches, lepton[channel]['lepton3'], CorrFactor)[cut]])
    Weights['nom'] = np.concatenate( [Weights['nom'], weights*FakesProp])
    Weights['l1_errUp'] = np.concatenate( [Weights['l1_errUp'], weights1_errUp*FakesProp])
    Weights['l2_errUp'] = np.concatenate( [Weights['l2_errUp'], weights2_errUp*FakesProp])
    Weights['l3_errUp'] = np.concatenate( [Weights['l3_errUp'], weights3_errUp*FakesProp])
    Weights['l1_errDown'] = np.concatenate( [Weights['l1_errDown'], weights1_errDown*FakesProp])
    Weights['l2_errDown'] = np.concatenate( [Weights['l2_errDown'], weights2_errDown*FakesProp])
    Weights['l3_errDown'] = np.concatenate( [Weights['l3_errDown'], weights3_errDown*FakesProp])

    for XXX_bis in application_region:
      if XXX_bis == XXX:
        Weights['FakeProp_errUp'][XXX_bis] = np.concatenate( [Weights['FakeProp_errUp'][XXX_bis], weights*(FakesProp+FakesProp_err)])
        Weights['FakeProp_errDown'][XXX_bis] = np.concatenate( [Weights['FakeProp_errDown'][XXX_bis], weights*(FakesProp-FakesProp_err)])
      else: 
        Weights['FakeProp_errUp'][XXX_bis] = np.concatenate( [Weights['FakeProp_errUp'][XXX_bis], weights*FakesProp])
        Weights['FakeProp_errDown'][XXX_bis] = np.concatenate( [Weights['FakeProp_errDown'][XXX_bis], weights*FakesProp])

  return Fakes, Weights

def load_XsecUnc(tag, period, channel):
  with open(os.path.join(os.getenv("RUN_PATH"), f'A_cutflow/results/{tag}/{period}/Xsec_unc.yml'), 'r') as f:
      XsecUnc = yaml.safe_load(f)
      XsecUnc = XsecUnc[channel]
  return XsecUnc

def main():
    eos_path = f'/eos/user/p/pdebryas/HNL/anatuple/'
    anatuple_path = os.path.join(eos_path, period, tag, channel , 'anatuple')

    output_path = f'{os.getenv("RUN_PATH")}/D_RootHist/results/'
    output_folder = os.path.join(output_path, period, tag, channel)
    os.makedirs(output_folder, exist_ok=True)

    #load global parameters
    with open(os.path.join(os.getenv("RUN_PATH"), f'common/config/all/config_FakeRate.yaml'), 'r') as f:
        GlobalParms = yaml.safe_load(f)

    ModelName = GlobalParms['ModelName']
    tag_LL = GlobalParms['tag_LL']
    lepton = GlobalParms['lepton']
    application_region = GlobalParms['application_region']
    nbin_max = GlobalParms['nbin_max']
    min_bck = GlobalParms['min_bck']
    max_rel_error = GlobalParms['max_rel_error']

    vars_config = get_var_config(channel)
    inputs_cfg = load_inputs_cfg('FakeRate', period, channel)
    input_files = load_inputs(inputs_cfg, anatuple_path)
    hnl_masses = get_hnl_masses(period)
    FR, FR_err = load_FR(tag, period, tag_LL)
    W_DY_ttbar = load_W_DY_ttbar(tag, period, 'SignalRegion', channel)
    CorrFactor = load_CorrFactor(tag_LL, period)
    XsecUnc = load_XsecUnc(tag, period, channel)

    print('Load files...')
    branches = load_ntuples(input_files, 'Events;1')

    auto_xbins = {}
    for var in vars_config.keys():
        for MassHNL_Hyp in hnl_masses:
            hist_name = f'{var}_HNLMass{MassHNL_Hyp}'
            x_bins = load_xbins(vars_config, var)
            print(f'Processing HNL mass hypothesis: {MassHNL_Hyp} GeV')
            if var.startswith('DNNscore'):
                branches = AddDNNscore(branches, channel, MassHNL_Hyp, period, ModelName)

            cutRegion = {}
            for process in ['TrueLepton', 'data', f'HNL{MassHNL_Hyp}']:
                if process == 'data':
                    cutRegion[process] = compute_region_mask(branches[process], channel, 'data', 'SignalRegion')
                else:
                    cutRegion[process] = compute_region_mask(branches[process], channel, 'MC', 'SignalRegion')

            inputs = {}
            print('--- Compute nominal ---')
            print('')
            print(f'- Computing Prompt Background ...')
            cut = cutRegion['TrueLepton'][f'SignalRegion_PassTightWP_TrueLeptons']
            inputs['TrueLepton'] = np.array(compute_var_to_plot(branches['TrueLepton'], var)).flatten()[cut]
            inputs['TrueLepton_w'] = branches['TrueLepton']['genWeight'][cut]

            print(f'- Computing FakeBackground ...')
            Fakes, Weights = apply_FR_method(FR, FR_err, branches, channel, var, cutRegion, 'SignalRegion', W_DY_ttbar, CorrFactor, lepton, application_region)
            inputs['FakeBackground'] = Fakes['var']
            inputs['FakeBackground_w'] = Weights['nom']

            print(f'- Computing Signal ...')
            cut = cutRegion[f'HNL{MassHNL_Hyp}'][f'SignalRegion_PassTightWP'] 
            inputs['signal'] = np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], var)).flatten()[cut]
            inputs['signal_w'] = branches[f'HNL{MassHNL_Hyp}']['genWeight'][cut]

            print(f'- Computing Data ...')
            cut = cutRegion['data'][f'SignalRegion_PassTightWP']
            inputs['data'] = np.array(compute_var_to_plot(branches['data'], var)).flatten()[cut]
            inputs['data_w'] = branches['data']['genWeight'][cut]

            print('--- Compute corrections ---')
            print('')
            print('- Compute corrections for TrueLepton')
            for MainMCsample in XsecUnc.keys():
                if XsecUnc[MainMCsample]['relunc'] != None:
                    cut = cutRegion['TrueLepton'][f'SignalRegion_PassTightWP_TrueLeptons']
                    mask_sample = branches['TrueLepton'][f'mask_{MainMCsample}'][cut]
                    weight_nom = branches['TrueLepton']['genWeight'][cut]
                    weight_up  = np.where(mask_sample == 1, (1. + XsecUnc[MainMCsample]['relunc'])*weight_nom, weight_nom)  
                    weight_down= np.where(mask_sample == 1, (1. - XsecUnc[MainMCsample]['relunc'])*weight_nom, weight_nom)  
                    MainMCsampleName = MainMCsample.replace('_','').replace('-','')
                    inputs[f'TrueLepton_XsecUnc{MainMCsampleName}Up'] = np.array(compute_var_to_plot(branches['TrueLepton'], hist_name)).flatten()[cut]
                    inputs[f'TrueLepton_XsecUnc{MainMCsampleName}Up'] = weight_up
                    inputs[f'TrueLepton_XsecUnc{MainMCsampleName}Down'] = np.array(compute_var_to_plot(branches['TrueLepton'], hist_name)).flatten()[cut]
                    inputs[f'TrueLepton_XsecUnc{MainMCsampleName}Down'] = weight_down
            #here
            inputs = produceNPhists(inputs, branches, 'TrueLepton', cut, x_bins, hist_name, channel)
            print('')
            print('- Compute corrections for FakeBackground')
            for XXX in application_region:
                print(f'   ... Add hist for FakeProp in region {XXX}')
                for ud in ['Up', 'Down']:
                    hists[f'FakeBackground_statFakeProp{XXX}{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
                    hists[f'FakeBackground_statFakeProp{XXX}{period}{ud}'].fill(x= Fakes['var'], weight=Weights[f'FakeProp_err{ud}'][XXX])
                Weights_syst, bins_syst = ComputeFakeFactorsSyst(channel, period, Fakes, tag, W_DY_ttbar)
                hists = AddFakeFactorsUnc(hists, channel, period, x_bins, Fakes['var'], Weights, Weights_syst, bins_syst, is_flow)

            print('')
            print(f'- Compute corrections for Signal')
            cut = cutRegion[f'HNL{MassHNL_Hyp}'][f'{PlotRegion}_PassTightWP']
            hists = produceNPhists(hists, branches, f'HNL{MassHNL_Hyp}', cut, x_bins, hist_name, channel, is_flow)
            print('')    
            auto_xbins[hist_name] = compute_adaptive_binning(inputs, x_bins, nbin_max, min_bck, max_rel_error)




    print('')
    print('- Compute ES corrections for (TrueLepton/signal)')
    for Treename in load_Tree(input_files):
      if Treename.startswith('Events_'):
        print(f'   ... Add hist for {Treename} ...')
        branches = load_ntuples(input_files, Treename)
          
        SFBranchName = Treename.replace('Events_', "").replace(';1', "")
        if SFBranchName.endswith('_up'): 
            SFBranchName = SFBranchName.replace('_up', 'Up')
        if SFBranchName.endswith('_down'):
            SFBranchName = SFBranchName.replace('_down', 'Down')
        SFBranchName = SFBranchName.replace('_', "")

        if hist_name[-7:] == 'tauhOSL':
          branches = AddOSLepton(branches, channel)

        if hist_name.startswith('DNNscore'):
          branches = AddDNNscore(branches, channel, MassHNL_Hyp, period)
 
        #True Tau MC background
        cut = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)[f'{PlotRegion}_PassTightWP_TrueLeptons']
        hists[f'TrueLepton_{SFBranchName}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f'TrueLepton_{SFBranchName}'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name)).flatten()[cut], weight=branches['TrueLepton']['genWeight'][cut])

        #signal
        cut = compute_region_mask(branches[f'HNL{MassHNL_Hyp}'], channel, 'MC', PlotRegion)[f'{PlotRegion}_PassTightWP']
        hists[f'HNL{MassHNL_Hyp}_{SFBranchName}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f'HNL{MassHNL_Hyp}_{SFBranchName}'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name)).flatten()[cut], weight=branches[f'HNL{MassHNL_Hyp}']['genWeight'][cut])



    # Save to a YAML file
    with open(f'{output_folder}/auto_xbins.yml', 'w') as file:
        yaml.dump(auto_xbins, file, default_flow_style=False)

if __name__ == "__main__":
    main()


def produceNPhists(hists, branches, process, cut_region, x_bins, hist_name, channel, is_flow):
  common_corr_list = ['weightcorr_Tau_TauID_genuineElectron_Total_Central', 
                      'weightcorr_Tau_TauID_genuineMuon_Total_Central', 
                      'weightcorr_Tau_TauID_genuineTau_Total_Central', 
                      'weightcorr_Muon_MuID_Total_Central', 
                      'weightcorr_Electron_EleID_Total_Central',
                      'weightcorr_Muon_TrgSF_singleMu_Total_Central',
                      'weightcorr_Electron_TrgSF_singleEle_stat_Total_Central',
                      'weightcorr_L1PreFiring_Total_Central',
                      'weightcorr_PileUp_Total_Central',
                      'weightcorr_bTagSF_Loose_Total_Central'
                     ]

  combine_corr_list = [
                      'weightcorr_Tau1_TauID_genuineElectron_Total_Central',
                      'weightcorr_Tau1_TauID_genuineMuon_Total_Central',
                      'weightcorr_Tau1_TauID_genuineTau_Total_Central',
                      'weightcorr_Muon1_MuID_Total_Central', 
                      'weightcorr_Electron1_EleID_Total_Central', 
                      'weightcorr_Muon1_TrgSF_singleMu_Total_Central', 
                      'weightcorr_Electron1_TrgSF_singleEle_stat_Total_Central'
                     ]
  
  print(f'For {process} ...')
  for field in branches[process].fields:
    if field in common_corr_list:
      #print(field)
      name = field.replace('_Total_Central', '')
      NPs = [element for element in branches[process].fields if element.startswith(name)]
      for NP in NPs:
        name_TH1 = NP.replace('weightcorr_', '').replace('_', '')
        name_TH1 = name_TH1.replace('MuonMuID', 'MuID')
        name_TH1 = name_TH1.replace('ElectronEleID', 'EleID')
        name_TH1 = name_TH1.replace('TauTauID', 'TauID')
        if name_TH1.endswith('Uprel') or name_TH1.endswith('Downrel'):
          name_TH1 = name_TH1.replace('Uprel', 'Up')
          name_TH1 = name_TH1.replace('Downrel', 'Down')
          hists[f'{process}_{name_TH1}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
          if len(branches[process][NP][cut_region]) > 0:
            if isinstance(branches[process][NP][cut_region][0],  ak.Array):
              corr_weights = branches[process]['genWeight'][cut_region] * ak.concatenate(branches[process][NP][cut_region], axis = -1)
            else:
              corr_weights = np.array(branches[process]['genWeight'][cut_region] * branches[process][NP][cut_region]).flatten()
          else:
            corr_weights = np.array(branches[process]['genWeight'][cut_region] * branches[process][NP][cut_region]).flatten()
          #if isinstance(branches[process][NP][cut_region][0],  ak.Array):
          #  corr_weights = branches[process]['genWeight'][cut_region] * ak.concatenate(branches[process][NP][cut_region], axis = -1)
          #else:
          #  corr_weights = np.array(branches[process]['genWeight'][cut_region] * branches[process][NP][cut_region]).flatten()
          #corr_weights = branches[process]['genWeight'][cut_region] * np.array(branches[process][NP][cut_region]).flatten()
          hists[f'{process}_{name_TH1}'].fill(x=np.array(compute_var_to_plot(branches[process], hist_name)).flatten()[cut_region], weight=corr_weights)

    if field in combine_corr_list:
      #print(field)
      name = field.replace('_Total_Central', '')
      NPs = [element for element in branches[process].fields if element.startswith(name)]
      for NP1 in NPs:
        name_TH1 = NP1.replace('weightcorr_', '').replace('_', '')
        name_TH1 = name_TH1.replace('Tau1', '')
        name_TH1 = name_TH1.replace('Muon1', '')
        name_TH1 = name_TH1.replace('Electron1', '')
        if name_TH1.endswith('Uprel') or name_TH1.endswith('Downrel'):
          name_TH1 = name_TH1.replace('Uprel', 'Up')
          name_TH1 = name_TH1.replace('Downrel', 'Down')
          NP2 = NP1.replace('Tau1', 'Tau2')
          NP2 = NP1.replace('Muon1', 'Muon2')
          NP2 = NP1.replace('Electron1', 'Electron2')
          hists[f'{process}_{name_TH1}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
          if len(branches[process][NP1][cut_region]) > 0:
            if isinstance(branches[process][NP1][cut_region][0],  ak.Array):
              corr_weights = branches[process]['genWeight'][cut_region] * ak.concatenate(branches[process][NP1][cut_region] * branches[process][NP2][cut_region], axis = -1)
            else:
              corr_weights = np.array(branches[process]['genWeight'][cut_region] * branches[process][NP1][cut_region] * branches[process][NP2][cut_region]).flatten()
          else:
            corr_weights = np.array(branches[process]['genWeight'][cut_region] * branches[process][NP1][cut_region] * branches[process][NP2][cut_region]).flatten()
          hists[f'{process}_{name_TH1}'].fill(x=np.array(compute_var_to_plot(branches[process], hist_name)).flatten()[cut_region], weight=corr_weights)

    if (field not in common_corr_list) & (field not in combine_corr_list) & (field.endswith('_Total_Central')):
      name = field.replace('_Total_Central', '')
      name = name.replace('weightcorr_', '')
      split_name  = name.split("_")
      if (name not in ['Tau_TauID', 'Tau1_TauID']) & (split_name[0] not in ['Tau2','Muon2','Electron2']):
        print(f'WARNING: field {field} not taken into account !!!!!!!!!!!!!!!!!!!!!!!!!!!!')

  return hists