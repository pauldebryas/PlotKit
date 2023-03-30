import numpy as np
import ROOT

def compute_var_to_plot(tree, hist_name):
  
  if hist_name == 'mT_tautau':
    var_to_plot = mT_l1l2(tree, 'Tau1', 'Tau2')

  if hist_name == 'mT_total_tautau':
    var_to_plot = mT_total_l1l2(tree, 'Tau1', 'Tau2')
  
  if hist_name == 'dr_etau':
    var_to_plot = delta_r(tree, 'Electron', 'Tau')

  if hist_name == 'dr_mutau':
    var_to_plot = delta_r(tree, 'Muon', 'Tau1')
    
  if hist_name == 'pt_sum_tee':
    var_to_plot = pt_sum_l1l2l3(tree, 'Tau', 'Electron1', 'Electron2')

  if hist_name == 'pt_sum_tem':
    var_to_plot = pt_sum_l1l2l3(tree, 'Electron', 'Muon', 'Tau')

  if hist_name == 'pt_sum_tmm':
    var_to_plot = pt_sum_l1l2l3(tree, 'Tau', 'Muon1', 'Muon2')

  if hist_name == 'pt_sum_tte':
    var_to_plot = pt_sum_l1l2l3(tree, 'Electron', 'Tau1', 'Tau2')

  if hist_name == 'pt_sum_ttm':
    var_to_plot = pt_sum_l1l2l3(tree, 'Muon', 'Tau1', 'Tau2')

  if hist_name[-7:] == '_charge':
    var_to_plot = lepton_charge(tree, hist_name[:-7])

  if hist_name[-3:] == '_pt':
    var_to_plot = lepton_pt(tree, hist_name[:-3])

  if hist_name[-4:] == '_phi':
    var_to_plot = lepton_phi(tree, hist_name[:-4])

  if hist_name[-5:] == '_mass':
    var_to_plot = lepton_mass(tree, hist_name[:-5])

  if hist_name[-4:] == '_eta':
    var_to_plot = lepton_eta(tree, hist_name[:-4])

  if hist_name[-3:] == '_dz':
    var_to_plot = lepton_dz(tree, hist_name[:-3])

  if hist_name[-4:] == '_dxy':
    var_to_plot = lepton_dxy(tree, hist_name[:-4])

  return var_to_plot

def lepton_charge(tree, Lepton):
    return tree[f'{Lepton}_charge'].array()

def lepton_pt(tree, Lepton):
    return tree[f'{Lepton}_pt'].array()

def lepton_phi(tree, Lepton):
    return tree[f'{Lepton}_phi'].array()

def lepton_mass(tree, Lepton):
    return tree[f'{Lepton}_mass'].array()

def lepton_eta(tree, Lepton):
    return tree[f'{Lepton}_eta'].array()

def lepton_dz(tree, Lepton):
    return tree[f'{Lepton}_dz'].array()

def lepton_dxy(tree, Lepton):
    return tree[f'{Lepton}_dxy'].array()

def delta_r2(tree, Lepton1, Lepton2):
  '''
  Calculates deltaR squared between two particles v1, v2 
  ''' 
  v1_phi = tree[f'{Lepton1}_phi'].array()
  v1_eta = tree[f'{Lepton1}_eta'].array()
  v2_phi = tree[f'{Lepton2}_phi'].array()
  v2_eta = tree[f'{Lepton2}_eta'].array()

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
  return tree[f'{Lepton1}_pt'].array()+tree[f'{Lepton2}_pt'].array()+tree[f'{Lepton3}_pt'].array()

def mT_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates Tranverse mass between two particles v1, v2 
  ''' 
  v1_mass = tree[f'{Lepton1}_mass'].array()
  v2_mass = tree[f'{Lepton2}_mass'].array()
  v1_pt = tree[f'{Lepton1}_pt'].array()
  v2_pt = tree[f'{Lepton2}_pt'].array()
  v1_phi = tree[f'{Lepton1}_phi'].array()
  v2_phi = tree[f'{Lepton2}_phi'].array()
  return np.sqrt( v1_mass**2 + v2_mass**2 + 2*(np.sqrt(v1_mass**2 + v1_pt**2)*np.sqrt(v2_mass**2 + v2_pt**2) - v1_pt*v2_pt*np.cos(abs(v1_phi - v2_phi))))

def mT_total_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates total Tranverse mass between two particles v1, v2 (observable from doi:10.1007/JHEP11(2014)056)
  ''' 
  v1_pt = tree[f'{Lepton1}_pt'].array()
  v2_pt = tree[f'{Lepton2}_pt'].array()
  MET_pt= tree['MET_pt'].array()
  v1_phi = tree[f'{Lepton1}_phi'].array()
  v2_phi = tree[f'{Lepton2}_phi'].array()
  MET_phi= tree['MET_phi'].array()

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