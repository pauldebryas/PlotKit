import hist
import os
import correctionlib.convert
import numpy as np
import awkward as ak
import yaml
import json

from common.helpers import load_ntuples
from D_RootHist.hist_makers.helpers import ToRootHist, compute_var_to_plot, load_inputs, load_xbins, load_Tree#, load_ntuples
from common.regions.regions import compute_region_mask

with open(os.path.join(os.getenv("RUN_PATH"), f'common/config/all/config_FakeRate.yaml'), 'r') as f:
    GlobalParms = yaml.safe_load(f)

print_ARinfo = GlobalParms['print_ARinfo']
ModelName = GlobalParms['ModelName']
tag_LL = GlobalParms['tag_LL']
lepton = GlobalParms['lepton']
application_region = GlobalParms['application_region']
nbin_max = GlobalParms['nbin_max']
min_bck = GlobalParms['min_bck']
max_rel_error = GlobalParms['max_rel_error']

def make_histograms(input_dir, hist_name=None, hist_cfg=None, inputs_cfg=None, channel =None, period=None, PlotRegion = False, tag =None):

  if '_HNLMass' in hist_name:
      hist_name_split = hist_name.split('_HNLMass')
      hist_name = hist_name_split[0]
      MassHNL_Hyp = int(hist_name_split[-1])
  else:
      MassHNL_Hyp = '300'

  out = ''
  hists = {}

  if PlotRegion == 'SignalRegion':
    is_correction = True
  else:
    is_correction = False

  # is_correction = False
    
  x_bins = load_xbins(hist_cfg, hist_name)
  is_flow = hist_cfg[hist_name].get('flow', False)
  input_files = load_inputs(inputs_cfg, input_dir)
  FR, FR_err = load_FR(tag, period)
  W_DY_ttbar = load_W_DY_ttbar(tag, period, PlotRegion, channel)
  XsecUnc = load_XsecUnc(tag, period, channel)
  CorrFactor = load_CorrFactor(tag, period)

  print('Load files...')
  branches = load_ntuples(input_files, 'Events;1')
  if hist_name[-7:] == 'tauhOSL':
    branches = AddOSLepton(branches, channel)
  if hist_name.startswith('DNNscore'):
    branches = AddDNNscore(branches, channel, MassHNL_Hyp, period)

  cutRegion = {}
  for process in ['TrueLepton', 'data', f'HNL{MassHNL_Hyp}']:
    if process == 'data':
      cutRegion[process] = compute_region_mask(branches[process], channel, 'data', PlotRegion)
    else:
      cutRegion[process] = compute_region_mask(branches[process], channel, 'MC', PlotRegion)

  print(f'Computing Prompt Background ...')
  cut = cutRegion['TrueLepton'][f'{PlotRegion}_PassTightWP_TrueLeptons']
  TrueLepton = np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut]
  TrueLepton_w = branches['TrueLepton']['genWeight'][cut]

  print(f'Computing FakeBackground ...')
  Fakes, Weights, out = apply_FR_method(FR, FR_err, branches, channel, hist_name, cutRegion, PlotRegion, W_DY_ttbar, CorrFactor, out)
  FakeBackground = Fakes['var']
  FakeBackground_w = Weights['nom']

  print(f'Computing Signal ...')
  cut = cutRegion[f'HNL{MassHNL_Hyp}'][f'{PlotRegion}_PassTightWP'] 
  Signal = np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut]
  Signal_w = branches[f'HNL{MassHNL_Hyp}']['genWeight'][cut]

  print(f'Computing Data ...')
  cut = cutRegion['data'][f'{PlotRegion}_PassTightWP']
  Data = np.array(compute_var_to_plot(branches['data'], hist_name, CorrFactor = CorrFactor)).flatten()[cut]
  Data_w = branches['data']['genWeight'][cut]

  if len(x_bins) == 2:
    print('--- Auto binning ---')
    #print(f'Data: {Data}')
    print(f'Number of signal events: {len(Signal)}')
    #print(f'FakeBackground: {FakeBackground}')
    #print(f'Number of FakeBackground events: {len(FakeBackground)}')
    print(f'sum w of FakeBackground: {np.sum(FakeBackground_w)}')
    #print(f'TrueLepton : {TrueLepton}')
    #print(f'Number of TrueLepton events: {len(TrueLepton)}')
    print(f'sum w of TrueLepton: {np.sum(TrueLepton_w)}')
    x_bins = compute_adaptive_binning(Signal, TrueLepton, TrueLepton_w, FakeBackground, FakeBackground_w, x_bins, nbin_max, min_bck, max_rel_error)

  #True Tau MC background
  hists['TrueLepton'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['TrueLepton'].fill(x=TrueLepton, weight=TrueLepton_w)
  # FakeBackground (FR method)
  hists['FakeBackground'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['FakeBackground'].fill(x= FakeBackground, weight=FakeBackground_w)
  #signal
  hists[f'HNL{MassHNL_Hyp}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists[f'HNL{MassHNL_Hyp}'].fill(x=Signal, weight=Signal_w)
  #data
  hists['data'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
  hists['data'].fill(x=Data, weight=Data_w)

  #compute corrections 
  if is_correction:
    print('--- Compute corrections ---')
    print('')
    print('- Compute corrections for TrueLepton')
    # for MainMCsample in XsecUnc.keys():
    #   if XsecUnc[MainMCsample]['relunc'] != None:
    #     cut = cutRegion['TrueLepton'][f'{PlotRegion}_PassTightWP_TrueLeptons']
    #     mask_sample = branches['TrueLepton'][f'mask_{MainMCsample}'][cut]
    #     weight_nom = branches['TrueLepton']['genWeight'][cut]
    #     weight_up  = np.where(mask_sample == 1, (1. + XsecUnc[MainMCsample]['relunc'])*weight_nom, weight_nom)  
    #     weight_down= np.where(mask_sample == 1, (1. - XsecUnc[MainMCsample]['relunc'])*weight_nom, weight_nom)  
    #     MainMCsampleName = MainMCsample.replace('_','').replace('-','')
    #     hists[f'TrueLepton_XsecUnc{MainMCsampleName}Up'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #     hists[f'TrueLepton_XsecUnc{MainMCsampleName}Up'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=weight_up)
    #     hists[f'TrueLepton_XsecUnc{MainMCsampleName}Down'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #     hists[f'TrueLepton_XsecUnc{MainMCsampleName}Down'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=weight_down)


    # extract DY list of files (put in a function)
    configMainMCbackgroundFile = os.path.join(os.getenv("RUN_PATH"), f'common/config/all/inputs/inputs_MainMCbackground_{channel}.yaml')
    with open(configMainMCbackgroundFile, 'r') as f:
      configMainMCbackground = yaml.safe_load(f)
    DY_entry = next((item for item in configMainMCbackground if item.get("name") == "DY"), None)
    if DY_entry:
        DYsamplelist = [f.replace("_anatuple.root", "") for f in DY_entry.get("files", [])]
    # ------ 

    cut = cutRegion['TrueLepton'][f'{PlotRegion}_PassTightWP_TrueLeptons']
    weight_nom = branches['TrueLepton']['genWeight'][cut]
    TrueLep_events = np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut]
    # initialize "other sample" mask to all true
    #mask_dummy = branches['TrueLepton'][f'mask_WWW'][cut]
    mask_samplecovred = np.array(np.zeros(np.sum(cut)), dtype=bool)
    for MainMCsample in XsecUnc.keys():
      print(MainMCsample)
      if (MainMCsample != "DYJetsToLL"):
        mask_sample = branches['TrueLepton'][f'mask_{MainMCsample}'][cut]
        mask_samplecovred = mask_samplecovred | np.array(mask_sample, dtype=bool)

      elif MainMCsample == "DYJetsToLL":
        mask_DYsample = np.array(np.zeros(np.sum(cut)), dtype=bool)
        for DYfile in DYsamplelist:
          mask_DYfile = branches['TrueLepton'][f'mask_{DYfile}'][cut]
          mask_samplecovred = mask_samplecovred | np.array(mask_DYfile, dtype=bool)
          mask_DYsample = mask_DYsample | np.array(mask_DYfile, dtype=bool)
        mask_sample = mask_DYsample

      #print(np.sum(mask_sample))
      #print(np.sum(mask_samplecovred))

      weight_up  = np.where(mask_sample == 1, (1. + XsecUnc[MainMCsample]['RelUncUp'])*weight_nom, weight_nom)  
      weight_down= np.where(mask_sample == 1, (1. - XsecUnc[MainMCsample]['RelUncDown'])*weight_nom, weight_nom)  
      MainMCsampleName = MainMCsample.replace('_','').replace('-','')
      hists[f'TrueLepton_XsecUnc{MainMCsampleName}Up'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'TrueLepton_XsecUnc{MainMCsampleName}Up'].fill(x=TrueLep_events, weight=weight_up)
      hists[f'TrueLepton_XsecUnc{MainMCsampleName}Down'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'TrueLepton_XsecUnc{MainMCsampleName}Down'].fill(x=TrueLep_events, weight=weight_down)

    weight_up  = np.where(mask_samplecovred == 0, (1. + 0.2)*weight_nom, weight_nom)  
    weight_down= np.where(mask_samplecovred == 0, (1. - 0.2)*weight_nom, weight_nom)  
    hists[f'TrueLepton_XsecUncOthersUp'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    hists[f'TrueLepton_XsecUncOthersUp'].fill(x=TrueLep_events, weight=weight_up)
    hists[f'TrueLepton_XsecUncOthersDown'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    hists[f'TrueLepton_XsecUncOthersDown'].fill(x=TrueLep_events, weight=weight_down)

    #statweight_up, statweight_down = compute_stat_unc(branches['TrueLepton']['genWeight'][cut])
    #hists[f'TrueLepton_StatUncUp'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #hists[f'TrueLepton_StatUncUp'].fill(x=TrueLep_events, weight=statweight_up)
    #hists[f'TrueLepton_StatUncDown'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #hists[f'TrueLepton_StatUncDown'].fill(x=TrueLep_events, weight=statweight_down)

    hists = produceNPhists(hists, branches, 'TrueLepton', cut, x_bins, hist_name, channel, is_flow, CorrFactor)

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
    hists = produceNPhists(hists, branches, f'HNL{MassHNL_Hyp}', cut, x_bins, hist_name, channel, is_flow, CorrFactor)

    # Add Signal Theory Unc.
    Signal_w_up  = (1. + 0.04)*Signal_w
    Signal_w_down= (1. - 0.04)*Signal_w
    hists[f'HNL{MassHNL_Hyp}_PDFuncUp'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    hists[f'HNL{MassHNL_Hyp}_PDFuncUp'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=Signal_w_up)
    hists[f'HNL{MassHNL_Hyp}_PDFuncDown'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    hists[f'HNL{MassHNL_Hyp}_PDFuncDown'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=Signal_w_down)
  
    #statweight_up, statweight_down = compute_stat_unc(Signal_w)
    #hists[f'HNL{MassHNL_Hyp}_StatUncUp'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #hists[f'HNL{MassHNL_Hyp}_StatUncUp'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=statweight_up)
    #hists[f'HNL{MassHNL_Hyp}_StatUncDown'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
    #hists[f'HNL{MassHNL_Hyp}_StatUncDown'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=statweight_down)

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
        hists[f'TrueLepton_{SFBranchName}'].fill(x=np.array(compute_var_to_plot(branches['TrueLepton'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=branches['TrueLepton']['genWeight'][cut])

        #signal
        cut = compute_region_mask(branches[f'HNL{MassHNL_Hyp}'], channel, 'MC', PlotRegion)[f'{PlotRegion}_PassTightWP']
        hists[f'HNL{MassHNL_Hyp}_{SFBranchName}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f'HNL{MassHNL_Hyp}_{SFBranchName}'].fill(x=np.array(compute_var_to_plot(branches[f'HNL{MassHNL_Hyp}'], hist_name, CorrFactor = CorrFactor)).flatten()[cut], weight=branches[f'HNL{MassHNL_Hyp}']['genWeight'][cut])

  root_hists = { hist_name : ToRootHist(h, hist_name, flow= is_flow) for hist_name, h in hists.items() }

  if print_ARinfo:
    print(out)
    
  filename = f'AR_info_{PlotRegion}'
  file_path = os.path.join(os.getenv("RUN_PATH"), f'D_RootHist/results/{period}/{tag}/{channel}/{filename}.txt')
  with open(file_path, 'w') as file:
    file.write(out)
    file.write('\n')

  return root_hists

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
  weights = W_DY_ttbar[Lepton["type"]]['DY'][0]*weights_DY + W_DY_ttbar[Lepton["type"]]['ttbar'][0]*weights_ttbar
  weights_err = np.sqrt((W_DY_ttbar[Lepton["type"]]['DY'][0] * weights_DY_err) ** 2 + (W_DY_ttbar[Lepton["type"]]['ttbar'][0]* weights_ttbar_err) ** 2 + (weights_DY * W_DY_ttbar[Lepton["type"]]['DY'][1]) ** 2 + (weights_ttbar * W_DY_ttbar[Lepton["type"]]['ttbar'][1]) ** 2)

  return weights, weights_err

def compute_FR(weight, weight_err):
    nom = weight/(1-weight)
    err = weight_err/((1-weight)**2)
    up = np.where(nom + err<0, 0, nom+err)
    down = np.where(nom - err<0, 0, nom-err)
    return nom, up, down

def AddOSLepton(branches, channel):
  info_to_save= ['mass', 'pt', 'phi']
  for process in branches.keys():
    mask1 = branches[process]["Tau_charge"] != branches[process][f"{lepton[channel]['lepton1']['name']}_charge"]
    mask2 = branches[process]["Tau_charge"] != branches[process][f"{lepton[channel]['lepton2']['name']}_charge"]
    mask_total = mask1 | mask2
    for str in info_to_save:
      branches[process][f"OSLepton_{str}"] = np.where(mask1, branches[process][f"{lepton[channel]['lepton1']['name']}_{str}"], branches[process][f"{lepton[channel]['lepton2']['name']}_{str}"])
    print(f'filtering {len(mask_total) - np.sum(mask_total)} events out of {len(mask_total)}')
    branches[process] = branches[process][mask_total]
  return branches

def AddDNNscore(branches, channel, MassHNL_Hyp, period):
  if channel in ['tee_ss', 'tee_os']:
    channel = 'tee'
  if channel in ['tmm_ss', 'tmm_os']:
    channel = 'tmm'
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

def produceNPhists(hists, branches, process, cut_region, x_bins, hist_name, channel, is_flow, CorrFactor):
  common_corr_list = ['weightcorr_Tau_TauID_genuineElectron_Total_Central', 
                      'weightcorr_Tau_TauID_genuineMuon_Total_Central', 
                      'weightcorr_Tau_TauID_genuineTau_Total_Central', 
                      'weightcorr_Muon_MuID_Total_Central', 
                      'weightcorr_Electron_EleID_Total_Central',
                      'weightcorr_Muon_TrgSF_singleMu_Total_Central',
                      'weightcorr_Electron_TrgSF_singleEle_stat_Total_Central',
                      'weightcorr_L1PreFiring_ECAL_Total_Central',
                      'weightcorr_L1PreFiring_Muon_Total_Central',
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
          hists[f'{process}_{name_TH1}'].fill(x=np.array(compute_var_to_plot(branches[process], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region], weight=corr_weights)

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
          hists[f'{process}_{name_TH1}'].fill(x=np.array(compute_var_to_plot(branches[process], hist_name, CorrFactor = CorrFactor)).flatten()[cut_region], weight=corr_weights)

    if (field not in common_corr_list) & (field not in combine_corr_list) & (field.endswith('_Total_Central')):
      name = field.replace('_Total_Central', '')
      name = name.replace('weightcorr_', '')
      split_name  = name.split("_")
      if (name not in ['Tau_TauID', 'Tau1_TauID']) & (split_name[0] not in ['Tau2','Muon2','Electron2']):
        print(f'WARNING: field {field} not taken into account !!!!!!!!!!!!!!!!!!!!!!!!!!!!')

  return hists

def apply_FR_method(FR, FR_err, branches, channel, hist_name, cutRegion, PlotRegion, W_DY_ttbar, CorrFactor, out):
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
    out += f'{XXX} region: \n'
    N_events_XXX = np.sum(cut)
    out += f'N_events AppRegion= {N_events_XXX} \n'

    if XXX == 'FFF':
        weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
        weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
        weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)

        f1, f1_up, f1_down = compute_FR(weights1, weights1_err)
        f2, f2_up, f2_down = compute_FR(weights2, weights2_err)
        f3, f3_up, f3_down = compute_FR(weights3, weights3_err)

        weights          = f1 * f2 * f3
        weights1_errUp   = f1_up   * f2     * f3
        weights1_errDown = f1_down * f2     * f3
        weights2_errUp   = f1      * f2_up  * f3
        weights2_errDown = f1      * f2_down* f3
        weights3_errUp   = f1      * f2     * f3_up
        weights3_errDown = f1      * f2     * f3_down


    if XXX == 'FFP':
        weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
        weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)

        f1, f1_up, f1_down = compute_FR(weights1, weights1_err)
        f2, f2_up, f2_down = compute_FR(weights2, weights2_err)

        weights          = -1 * f1 * f2
        weights1_errUp   = -1 * f1_up   * f2
        weights1_errDown = -1 * f1_down * f2
        weights2_errUp   = -1 * f1      * f2_up
        weights2_errDown = -1 * f1      * f2_down
        weights3_errUp   = weights      
        weights3_errDown = weights      


    if XXX == 'FPF':
        weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)
        weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)

        f1, f1_up, f1_down = compute_FR(weights1, weights1_err)
        f3, f3_up, f3_down = compute_FR(weights3, weights3_err)

        weights          = -1 * f1 * f3
        weights1_errUp   = -1 * f1_up   * f3
        weights1_errDown = -1 * f1_down * f3
        weights3_errUp   = -1 * f1      * f3_up
        weights3_errDown = -1 * f1      * f3_down
        weights2_errUp   = weights      
        weights2_errDown = weights      


    if XXX == 'PFF':
        weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)
        weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)

        f2, f2_up, f2_down = compute_FR(weights2, weights2_err)
        f3, f3_up, f3_down = compute_FR(weights3, weights3_err)

        weights          = -1 * f2 * f3
        weights2_errUp   = -1 * f2_up   * f3
        weights2_errDown = -1 * f2_down * f3
        weights3_errUp   = -1 * f2      * f3_up
        weights3_errDown = -1 * f2      * f3_down
        weights1_errUp   = weights      
        weights1_errDown = weights      


    if XXX == 'FPP':
        weights1, weights1_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton1"], cut, W_DY_ttbar, CorrFactor)

        f1, f1_up, f1_down = compute_FR(weights1, weights1_err)

        weights          = f1
        weights1_errUp   = f1_up
        weights1_errDown = f1_down
        weights2_errUp   = weights   
        weights2_errDown = weights
        weights3_errUp   = weights
        weights3_errDown = weights


    if XXX == 'PFP':
        weights2, weights2_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton2"], cut, W_DY_ttbar, CorrFactor)

        f2, f2_up, f2_down = compute_FR(weights2, weights2_err)

        weights         = f2
        weights2_errUp   = f2_up
        weights2_errDown = f2_down
        weights1_errUp   = weights   
        weights1_errDown = weights
        weights3_errUp   = weights
        weights3_errDown = weights


    if XXX == 'PPF':
        weights3, weights3_err = compute_FRweight(FR, FR_err, branches, lepton[channel]["lepton3"], cut, W_DY_ttbar, CorrFactor)

        f3, f3_up, f3_down = compute_FR(weights3, weights3_err)

        weights         = f3
        weights3_errUp   = f3_up
        weights3_errDown = f3_down
        weights1_errUp   = weights   
        weights1_errDown = weights
        weights2_errUp   = weights
        weights2_errDown = weights


    out += f'sumw FFs= {round(np.sum(weights),1)} \n'
    FakesProp = (N_events_XXX - sumw_MC[f'AppRegion{XXX}_TrueLepton'])/N_events_XXX
    if FakesProp < 0:
      FakesProp = 0
    FakesProp_err = np.sqrt(N_events_XXX)/N_events_XXX

    out += f'sumw FakesProp = {round(FakesProp,1)} \n'
    Fakes['var'] = np.concatenate( [Fakes['var'], np.array(compute_var_to_plot(branches['data'], hist_name, CorrFactor = CorrFactor)[cut]).flatten()])
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

    out += f'FFs * FakesProp = {round(np.sum(weights)*FakesProp,1)} \n'
    out += f'\n'

  return Fakes, Weights, out

def AddFakeFactorsUnc(hists, channel, period, x_bins, Fakes, Weights, Weights_syst, bins_syst, is_flow):
  print(f' ... Add hist for FakeFactorsUnc for {channel} ...')
  nbins = len(bins_syst['lepton1']) # should be the same for all leptons

  if channel not in ['tee', 'tee_ss', 'tee_os', 'tmm', 'tmm_ss', 'tmm_os', 'tem', 'tte', 'ttm']:
    print(f'missing channel {channel}')
    raise
  
  if channel in ['tee', 'tee_ss', 'tee_os']:
    for ud in ['Up', 'Down']:
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'].fill(x= Fakes, weight=Weights[f'l1_err{ud}'])

      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'].fill(x= Fakes, weight=(Weights[f'l2_err{ud}']*Weights[f'l3_err{ud}'])/Weights['nom'])

      for i in range(nbins-1):
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton1'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton2'][i][ud]*Weights_syst['lepton3'][i][ud]*Weights['nom'])
          
  if channel in ['tmm', 'tmm_ss', 'tmm_os']:
    for ud in ['Up', 'Down']:
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'].fill(x= Fakes, weight=Weights[f'l1_err{ud}'])

      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'].fill(x= Fakes, weight=(Weights[f'l2_err{ud}']*Weights[f'l3_err{ud}'])/Weights['nom'])

      for i in range(nbins-1):
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton1'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton2'][i][ud]*Weights_syst['lepton3'][i][ud]*Weights['nom'])

  if channel == 'tem':
    for ud in ['Up', 'Down']:
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'].fill(x= Fakes, weight=Weights[f'l1_err{ud}'])

      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'].fill(x= Fakes, weight=Weights[f'l2_err{ud}'])

      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'].fill(x= Fakes,  weight=Weights[f'l3_err{ud}'])

      for i in range(nbins-1):
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton1'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton2'][i]}to{bins_syst['lepton2'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton2'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton3'][i][ud]*Weights['nom'])

  if channel == 'tte':
    for ud in ['Up', 'Down']:
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'].fill(x= Fakes, weight=(Weights[f'l1_err{ud}']*Weights[f'l2_err{ud}'])/Weights['nom'])

      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsElectron{period}{ud}'].fill(x= Fakes, weight=Weights[f'l3_err{ud}'])

      for i in range(nbins-1):
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton1'][i][ud]*Weights_syst['lepton2'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsElectronPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton3'][i][ud]*Weights['nom'])

  if channel == 'ttm':
    for ud in ['Up', 'Down']:
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsTau{period}{ud}'].fill(x= Fakes, weight=(Weights[f'l1_err{ud}']*Weights[f'l2_err{ud}'])/Weights['nom'])

      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
      hists[f'FakeBackground_statFakeFactorsMuon{period}{ud}'].fill(x= Fakes,  weight=Weights[f'l3_err{ud}'])

      for i in range(nbins-1):
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsTauPt{bins_syst['lepton1'][i]}to{bins_syst['lepton1'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton1'][i][ud]*Weights_syst['lepton2'][i][ud]*Weights['nom'])

        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"] = hist.Hist.new.Variable(x_bins, name='x', flow=is_flow).Weight()
        hists[f"FakeBackground_systFakeFactorsMuonPt{bins_syst['lepton3'][i]}to{bins_syst['lepton3'][i+1]}Bin{i}{period}{ud}"].fill(x= Fakes, weight=Weights_syst['lepton3'][i][ud]*Weights['nom'])

  return hists 

def ComputeFakeFactorsSyst(channel, period, Fakes, tag, W_DY_ttbar):
  # need to be updated: ptcorr for LL. Use DY and ttbar treated as independent source ans using DY/ttbar relative fraction
  TauSystRelErrorFile = os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/Syst_Unc_TauFF.yaml')
  LLSystRelErrorFile =  os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag_LL}/{period}/Syst_Unc_LLFF.yaml')
  with open(TauSystRelErrorFile, 'r') as f:
    TauSystRelError = yaml.safe_load(f)
  with open(LLSystRelErrorFile, 'r') as f:
    LLSystRelError = yaml.safe_load(f)
  SystRelError = {**TauSystRelError, **LLSystRelError}
  Weights_lep = {}
  Bins_lep = {}
  for leptonNb in ['lepton1', 'lepton2', 'lepton3']:
    Weights_lep[leptonNb] = {}
    lepton_type = lepton[channel][leptonNb]['type']
    errLepDY = SystRelError[lepton_type]['DY'] 
    errLepttbar = SystRelError[lepton_type]['ttbar']
    #binning should be the same for DY and ttbar
    bins = errLepDY['bin'] 
    bins[0] = 0
    bins[-1] = 1000000
    Bins_lep[leptonNb] = bins
    for i in range(len(bins) -1):
      #combine_err = (W_DY_ttbar[lepton_type]['DY'][0]*abs(errLepDY['ratio'][i])+W_DY_ttbar[lepton_type]['ttbar'][0]*abs(errLepttbar['ratio'][i]))
      combine_err = np.sqrt((W_DY_ttbar[lepton_type]['DY'][0]*errLepDY['ratio'][i])**2 + (W_DY_ttbar[lepton_type]['ttbar'][0] * errLepttbar['ratio'][i])**2)
      if combine_err > 0.3:
        print(f'WARNING: {lepton[channel][leptonNb]["name"]} bin {i} has a large error: {combine_err} | {errLepDY["ratio"][i]} (DY) and {errLepttbar["ratio"][i]} (ttbar)')
      Weights_lep[leptonNb][i] = {}
      mask_bin_i =  (bins[i] <= Fakes[f'{leptonNb}_ptcorr']) & (bins[i+1] > Fakes[f'{leptonNb}_ptcorr'])
      corr_weight_up =  np.ones(len(mask_bin_i))
      corr_weight_up[mask_bin_i] = (1. + combine_err)
      Weights_lep[leptonNb][i]['Up'] = corr_weight_up
      corr_weight_down =  np.ones(len(mask_bin_i))
      corr_weight_down[mask_bin_i] = (1. - combine_err)
      Weights_lep[leptonNb][i]['Down'] = corr_weight_down
  return Weights_lep, Bins_lep 

def load_FR(tag, period):
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
  with open(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FRweights_{PlotRegion}AppRegion.yml'), 'r') as stream: 
      file_content = yaml.safe_load(stream)
      # Collect only leptons that contain the requested channel
      W_DY_ttbar = {
          lepton: channels[channel]
          for lepton, channels in file_content.items()
          if channel in channels
      }
  return W_DY_ttbar

# def load_XsecUnc(tag, period, channel):
#   with open(os.path.join(os.getenv("RUN_PATH"), f'A_cutflow/results/{tag}/{period}/Xsec_unc.yml'), 'r') as f:
#       XsecUnc = yaml.safe_load(f)
#       XsecUnc = XsecUnc[channel]
#   return XsecUnc

def load_XsecUnc(tag, period, channel):
  with open(os.path.join(os.getenv("RUN_PATH"), f'A_cutflow/results/Xsec_unc.yml'), 'r') as f:
      XsecUnc = yaml.safe_load(f)
  return XsecUnc

def load_CorrFactor(tag, period):
  # load Corr factor from YAML file
  with open(os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag_LL}/{period}/correctionFactorsLL.yml'), 'r') as file:
      CorrFactor = yaml.safe_load(file)
  return CorrFactor

def compute_ptcorr(branches, lepton, CorrFactor):
  lepton_pt = np.array(branches['data'][f'{lepton["name"]}_pt']).flatten()
  if lepton['type'] != 'Tau':
    lepton_parton_pt = np.array(branches['data'][f'{lepton["name"]}_ConeCorrectedPt']).flatten()*CorrFactor[lepton['type']]
    lepton_iso = np.array(branches['data'][f'{lepton["name"]}_pfRelIso03_all']).flatten()
    lepton_pt = np.where(lepton_iso < 0.15, lepton_pt, lepton_parton_pt)
  return lepton_pt

def compute_adaptive_binning(signal_vals, bkg1_vals, bkg1_weights, bkg2_vals, bkg2_weights, x_range, nbin_max, w_threshold_bck, max_relative_bck_err):

    min_val, max_val = x_range

    mask = (signal_vals >= min_val) & (signal_vals <= max_val)
    if len(signal_vals[mask]) < nbin_max:
      if len(signal_vals[mask]) == 0:
        print("Warning: No signal values in the specified range. Using a single bin.")
        nbin_max = 1
        n_threshold_signal = len(signal_vals) + 1
      else:
        nbin_max = len(signal_vals[mask])
        n_threshold_signal = np.floor(len(signal_vals[mask]) / nbin_max).astype(int)
    else:
      n_threshold_signal = np.floor(len(signal_vals[mask]) / nbin_max).astype(int)

    # Granularity control
    step_int_scale = 1e4  # granularity (0.0001 bins)
    min_value_int = int(min_val * step_int_scale)
    max_value_int = int(max_val * step_int_scale)
    n_bins = max_value_int - min_value_int + 1

    # Convert values to int bin positions
    signal_vals_int = (signal_vals * step_int_scale).astype(int)
    bkg1_vals_int = (bkg1_vals * step_int_scale).astype(int)
    bkg2_vals_int = (bkg2_vals * step_int_scale).astype(int)

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
    for v, w in zip(bkg1_vals_int, bkg1_weights):
        idx = v - min_value_int
        if 0 <= idx < n_bins:
            bkg1_hist[idx] += w
            bkg1_w2_hist[idx] += w**2

    # Fill background2 histogram (weighted yield and sum of weights squared)
    for v, w in zip(bkg2_vals_int, bkg2_weights):
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
        print('') 
        print(f'Sig yield: {acc_sig}')
        print(f'bkg yield: {total_bkg_yield}')
        print(f'True bkg yield: {acc_bkg1}')
        print(f'Fake bkg yield: {acc_bkg2}')
        print(f'Rel. err. True bkg: {np.sqrt(acc_bkg1_w2) / acc_bkg1}')
        print(f'Rel. err. Fake bkg: {np.sqrt(acc_bkg2_w2) / acc_bkg2}')
        print('')
        edges.append(j)

        if j <= min_value_int:
            if (acc_sig >= n_threshold_signal*0.5 and total_bkg_yield >= w_threshold_bck and bkg1_err_ok and bkg2_err_ok) or ((len(edges) - 1) <= 1):
                break
            else:
                print('removing last edge')
                edges.remove(j)
                # NEW safeguard: keep removing edges until merged bin yields are valid
                while len(edges) > 1:
                    # convert to array indices
                    edge_int = edges[-2]
                    end_idx   = edge_int - min_value_int

                    # recompute merged background yields for this last bin
                    acc_bkg1_merged = np.sum(bkg1_hist[:end_idx+1])
                    acc_bkg2_merged = np.sum(bkg2_hist[:end_idx+1])
                    total_bkg_yield_merged = acc_bkg1_merged + acc_bkg2_merged

                    print("Check merged:")
                    print("  True bkg:", acc_bkg1_merged)
                    print("  Fake bkg:", acc_bkg2_merged)
                    print("  Total:", total_bkg_yield_merged)

                    if acc_bkg1_merged <= 0 or acc_bkg2_merged <= 0 or total_bkg_yield_merged <= 0:
                        print("Removing problematic edge due to negative merged yield")
                        edges.pop(-2)  # drop the second-to-last edge
                    else:
                        break
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

def compute_stat_unc(weights):
    sumw = np.sum(weights)
    sumw2 = np.sum(weights**2)
    delta = np.sqrt(sumw2)
    # relative uncertainty
    rel_unc = delta / sumw if sumw != 0 else 0.0
    weight_up = weights * (1.0 + rel_unc)
    weight_down = weights * (1.0 - rel_unc)
    return weight_up, weight_down