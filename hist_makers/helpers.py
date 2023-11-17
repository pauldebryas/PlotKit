import numpy as np
import ROOT
import uproot
import os

def compute_var_to_plot(tree, hist_name):
  
  if hist_name == 'mT_tautau':
    var_to_plot = mT_l1l2(tree, 'Tau1', 'Tau2')
    return var_to_plot

  if hist_name == 'mT_total_tautau':
    var_to_plot = mT_total_l1l2(tree, 'Tau1', 'Tau2')
    return var_to_plot
  
  if hist_name == 'dr_etau':
    var_to_plot = delta_r(tree, 'Electron', 'Tau')
    return var_to_plot

  if hist_name == 'dr_mutau':
    var_to_plot = delta_r(tree, 'Muon', 'Tau1')
    return var_to_plot
    
  if hist_name == 'pt_sum_tee':
    var_to_plot = pt_sum_l1l2l3(tree, 'Tau', 'Electron1', 'Electron2')
    return var_to_plot

  if hist_name == 'pt_sum_tem':
    var_to_plot = pt_sum_l1l2l3(tree, 'Electron', 'Muon', 'Tau')
    return var_to_plot

  if hist_name == 'pt_sum_tmm':
    var_to_plot = pt_sum_l1l2l3(tree, 'Tau', 'Muon1', 'Muon2')
    return var_to_plot

  if hist_name == 'pt_sum_tte':
    var_to_plot = pt_sum_l1l2l3(tree, 'Electron', 'Tau1', 'Tau2')
    return var_to_plot

  if hist_name == 'pt_sum_ttm':
    var_to_plot = pt_sum_l1l2l3(tree, 'Muon', 'Tau1', 'Tau2')
    return var_to_plot

  if hist_name[-7:] == '_Jet_pt':
    var_to_plot = lepton_Jet_pt(tree, hist_name[:-7])
    return var_to_plot
  
  if hist_name[-7:] == '_charge':
    var_to_plot = lepton_charge(tree, hist_name[:-7])
    return var_to_plot

  if hist_name[-3:] == '_pt':
    var_to_plot = lepton_pt(tree, hist_name[:-3])
    return var_to_plot

  if hist_name[-4:] == '_phi':
    var_to_plot = lepton_phi(tree, hist_name[:-4])
    return var_to_plot

  if hist_name[-5:] == '_mass':
    var_to_plot = lepton_mass(tree, hist_name[:-5])
    return var_to_plot

  if hist_name[-4:] == '_eta':
    var_to_plot = lepton_eta(tree, hist_name[:-4])
    return var_to_plot

  if hist_name[-3:] == '_dz':
    var_to_plot = lepton_dz(tree, hist_name[:-3])
    return var_to_plot

  if hist_name[-4:] == '_dxy':
    var_to_plot = lepton_dxy(tree, hist_name[:-4])
    return var_to_plot

  if hist_name[-10:] == '_decayMode':
    var_to_plot = lepton_decayMode(tree, hist_name[:-10])
    return var_to_plot

  print('Error compute_var_to_plot function do not find the variable')
  return 

def lepton_charge(tree, Lepton):
    return tree[f'{Lepton}_charge'] 

def lepton_pt(tree, Lepton):
    return tree[f'{Lepton}_pt']

def lepton_Jet_pt(tree, Lepton):
    return tree[f'{Lepton}_Jet_pt']

def lepton_phi(tree, Lepton):
    return tree[f'{Lepton}_phi'] 

def lepton_mass(tree, Lepton):
    return tree[f'{Lepton}_mass'] 

def lepton_eta(tree, Lepton):
    return tree[f'{Lepton}_eta'] 

def lepton_dz(tree, Lepton):
    return tree[f'{Lepton}_dz'] 

def lepton_dxy(tree, Lepton):
    return tree[f'{Lepton}_dxy'] 

def lepton_decayMode(tree, Lepton):
    return tree[f'{Lepton}_decayMode'] 

def delta_r2(tree, Lepton1, Lepton2):
  '''
  Calculates deltaR squared between two particles v1, v2 
  ''' 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v1_eta = tree[f'{Lepton1}_eta'] 
  v2_phi = tree[f'{Lepton2}_phi'] 
  v2_eta = tree[f'{Lepton2}_eta'] 

  dphi = (v1_phi - v2_phi + np.pi) % (2 * np.pi) - np.pi
  deta = v1_eta - v2_eta
  dr2 = dphi**2 + deta**2
  return dr2

def delta_r(tree, Lepton1, Lepton2):
  '''
  Calculates deltaR between two particles v1, v2
  '''
  return np.sqrt(delta_r2(tree, Lepton1, Lepton2))

def pt_sum_l1l2l3(tree, Lepton1, Lepton2, Lepton3):
  return tree[f'{Lepton1}_pt'] +tree[f'{Lepton2}_pt'] +tree[f'{Lepton3}_pt'] 

def mT_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates Tranverse mass between two particles v1, v2 
  ''' 
  v1_mass = tree[f'{Lepton1}_mass'] 
  v2_mass = tree[f'{Lepton2}_mass'] 
  v1_pt = tree[f'{Lepton1}_pt'] 
  v2_pt = tree[f'{Lepton2}_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v2_phi = tree[f'{Lepton2}_phi'] 
  return np.sqrt( v1_mass**2 + v2_mass**2 + 2*(np.sqrt(v1_mass**2 + v1_pt**2)*np.sqrt(v2_mass**2 + v2_pt**2) - v1_pt*v2_pt*np.cos(abs(v1_phi - v2_phi))))

def Z_cut(tree, Lepton1, Lepton2):
  mass_z = 91.2 #GeV
  interval = 15.
  ZCut = (mT_l1l2(tree, Lepton1, Lepton2) < (mass_z - interval)) |  (mT_l1l2(tree, Lepton1, Lepton2) > (mass_z + interval))
  return ZCut

def mT_total_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates total Tranverse mass between two particles v1, v2 (observable from doi:10.1007/JHEP11(2014)056)
  ''' 
  v1_pt = tree[f'{Lepton1}_pt'] 
  v2_pt = tree[f'{Lepton2}_pt'] 
  MET_pt= tree['MET_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v2_phi = tree[f'{Lepton2}_phi'] 
  MET_phi= tree['MET_phi'] 

  mT2_l1l2 = 2.*v1_pt*v2_pt*(1.-np.cos(abs(v1_phi - v2_phi)))
  mT2_l1MET = 2.*v1_pt*MET_pt*(1.-np.cos(abs(v1_phi - MET_phi)))
  mT2_l2MET = 2.*MET_pt*v2_pt*(1.-np.cos(abs(MET_phi - v2_phi)))
  return np.sqrt(mT2_l1l2 + mT2_l1MET +  mT2_l2MET)

def ListToVector(list, type="string"):
  vec = ROOT.std.vector(f"{type}")()
  for item in list:
    vec.push_back(item)
  return vec

def ToRootHist(h):
  x_bins = ListToVector(h.axes[0].edges, type="double")
  root_hist = ROOT.TH1D('', '', len(x_bins) - 1, x_bins.data())
  for i in range(len(x_bins) - 1):
    root_hist.SetBinContent(i + 1, h.values()[i])
    root_hist.SetBinError(i + 1, h.variances()[i] ** 0.5)
  root_hist.SetDirectory(0)
  return root_hist

def load_ntuples(samples):
  processes = samples.keys()
  branches = {}
  for p in processes:
      list = []
      for file_path in samples[p]:
          # Load the dataframes
          DataUproot = uproot.open(file_path)
          list.append(DataUproot['Event;1'].arrays())
      branches[p] = np.concatenate(list)
  return branches

def equalObs(x, nbin):
    #calculate equal-frequency bins 
    nlen = len(x)
    return np.array(np.interp(np.linspace(0, nlen, nbin + 1), np.arange(nlen), np.sort(x)))

def load_inputs(inputs_cfg, input_dir):
  input_files = {}
  for input in inputs_cfg:
    if 'files' in input.keys():
      files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
      for file in files_list:
        if os.path.isfile(file) == False:
          print('WARNING: ' + file + ' is missing')
          files_list.remove(file)
      input_files[input['name']] = files_list
  return input_files

def load_xbins(hist_cfg, hist_name):
  hist_desc = hist_cfg[hist_name]
  x_bins = hist_desc['x_bins']
  if not isinstance(x_bins, list):
    n_bins, bin_range = x_bins.split('|')
    start,stop = bin_range.split(':')
    x_bins = np.arange(int(start),int(stop), (int(stop)-int(start))/int(n_bins))
  return x_bins