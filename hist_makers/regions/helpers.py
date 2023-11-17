import numpy as np

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
