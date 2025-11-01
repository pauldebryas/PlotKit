import numpy as np

def mT_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates tranverse mass between two particles v1, v2 
  ''' 
  v1_mass = tree[f'{Lepton1}_mass'] 
  v2_mass = tree[f'{Lepton2}_mass'] 
  v1_pt = tree[f'{Lepton1}_pt'] 
  v2_pt = tree[f'{Lepton2}_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v2_phi = tree[f'{Lepton2}_phi'] 
  return np.sqrt( v1_mass**2 + v2_mass**2 + 2*(np.sqrt(v1_mass**2 + v1_pt**2)*np.sqrt(v2_mass**2 + v2_pt**2) - v1_pt*v2_pt*np.cos(abs(v1_phi - v2_phi))))

def mass_l1l2(tree, Lepton1, Lepton2):
  '''
  Calculates invariant mass between two particles v1, v2 
  ''' 
  v1_eta = tree[f'{Lepton1}_eta'] 
  v2_eta = tree[f'{Lepton2}_eta'] 
  v1_pt = tree[f'{Lepton1}_pt'] 
  v2_pt = tree[f'{Lepton2}_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v2_phi = tree[f'{Lepton2}_phi'] 
  return np.sqrt(2*v1_pt*v2_pt*(np.cosh(v1_eta-v2_eta) - np.cos(v1_phi-v2_phi)))

def mass_l1l2l3(tree, Lepton1, Lepton2, Lepton3):
    '''
    Calculates invariant mass between three particles v1, v2, v3
    '''
    v1_eta = tree[f'{Lepton1}_eta']
    v2_eta = tree[f'{Lepton2}_eta']
    v3_eta = tree[f'{Lepton3}_eta']
    
    v1_pt = tree[f'{Lepton1}_pt']
    v2_pt = tree[f'{Lepton2}_pt']
    v3_pt = tree[f'{Lepton3}_pt']
    
    v1_phi = tree[f'{Lepton1}_phi']
    v2_phi = tree[f'{Lepton2}_phi']
    v3_phi = tree[f'{Lepton3}_phi']
    
    v1_m = tree[f'{Lepton1}_mass'] 
    v2_m = tree[f'{Lepton2}_mass']
    v3_m = tree[f'{Lepton3}_mass']

    # Compute transverse energy
    v1_E = np.sqrt(v1_m**2 + v1_pt**2 * np.cosh(v1_eta)**2)
    v2_E = np.sqrt(v2_m**2 + v2_pt**2 * np.cosh(v2_eta)**2)
    v3_E = np.sqrt(v3_m**2 + v3_pt**2 * np.cosh(v3_eta)**2)

    # Compute momentum components
    v1_px, v1_py, v1_pz = v1_pt * np.cos(v1_phi), v1_pt * np.sin(v1_phi), v1_pt * np.sinh(v1_eta)
    v2_px, v2_py, v2_pz = v2_pt * np.cos(v2_phi), v2_pt * np.sin(v2_phi), v2_pt * np.sinh(v2_eta)
    v3_px, v3_py, v3_pz = v3_pt * np.cos(v3_phi), v3_pt * np.sin(v3_phi), v3_pt * np.sinh(v3_eta)

    # Compute invariant mass
    total_E = v1_E + v2_E + v3_E
    total_px = v1_px + v2_px + v3_px
    total_py = v1_py + v2_py + v3_py
    total_pz = v1_pz + v2_pz + v3_pz

    return np.sqrt(total_E**2 - (total_px**2 + total_py**2 + total_pz**2))

def Z_cut_OSSF(tree, Lepton1, Lepton2, interval = 10.):
  if Lepton1[:-1] != Lepton2[:-1]:
    print('!!!!! Warning: lepton1 and lepton2 must have the same flavor.')
  SSLepton_mask = (tree[f'{Lepton1}_charge'] == tree[f'{Lepton2}_charge'])

  mass_z = 91.2 #GeV
  ZCut = (mass_l1l2(tree, Lepton1, Lepton2) < (mass_z - interval)) |  (mass_l1l2(tree, Lepton1, Lepton2) > (mass_z + interval))
  
  return np.where(SSLepton_mask, SSLepton_mask, ZCut)

def Z_cut_l1l2l3(tree, Lepton1, Lepton2, Lepton3, interval = 15.):
  mass_z = 91.2 #GeV
  ZCut = (mass_l1l2l3(tree, Lepton1, Lepton2, Lepton3) < (mass_z - interval)) |  (mass_l1l2l3(tree, Lepton1, Lepton2, Lepton3) > (mass_z + interval))
  return ZCut

def LowMass_cut_OSSF(tree, Lepton1, Lepton2):
  if Lepton1[:-1] != Lepton2[:-1]:
    print('!!!!! Warning: lepton1 and lepton2 must have the same flavor.')
  SSLepton_mask = (tree[f'{Lepton1}_charge'] == tree[f'{Lepton2}_charge'])

  LowerCut = 5. #GeV
  LowMassCut = (mass_l1l2(tree, Lepton1, Lepton2) > LowerCut)
  
  return np.where(SSLepton_mask, SSLepton_mask, LowMassCut)

def MET_cut(tree, cut_nb, lowerThan = False):
  met_pt = tree[f'MET_pt'] 
  met_phi = tree[f'MET_phi'] 

  met_x = met_pt * np.cos(met_phi)
  met_y = met_pt * np.sin(met_phi)

  MET = np.sqrt(met_x**2 + met_y**2)

  if lowerThan == True:
    return MET <= cut_nb
  else:
    return MET >= cut_nb
  
def Z_cut(tree, Lepton1, Lepton2, interval = 15.):
  mass_z = 91.2 #GeV
  ZCut = (mass_l1l2(tree, Lepton1, Lepton2) < (mass_z - interval)) |  (mass_l1l2(tree, Lepton1, Lepton2) > (mass_z + interval))
  return ZCut

def delta_R(tree, Lepton1, Lepton2):
    '''
    Calculates Delta_R between two particles v1, v2
    '''
    v1_eta = tree[f'{Lepton1}_eta']
    v2_eta = tree[f'{Lepton2}_eta']
    v1_phi = tree[f'{Lepton1}_phi']
    v2_phi = tree[f'{Lepton2}_phi']
    
    delta_eta = v1_eta - v2_eta
    delta_phi = np.abs(v1_phi - v2_phi)
    delta_phi = np.where(delta_phi > np.pi, 2*np.pi - delta_phi, delta_phi)
    
    return np.sqrt(delta_eta**2 + delta_phi**2)

def compute_invariant_mass_MET_lepton(tree, Lepton1):
    # Extract kinematic variables
    v1_pt = tree[f'{Lepton1}_pt']
    v1_phi = tree[f'{Lepton1}_phi']
    MET_pt = tree['MET_pt']
    MET_phi = tree['MET_phi']

    # Compute Δφ (angle between lepton and MET)
    delta_phi = np.abs(v1_phi - MET_phi)
    delta_phi = np.arctan2(np.sin(delta_phi), np.cos(delta_phi))  # Keep Δφ in [-π, π]

    # Compute invariant mass using transverse quantities
    M_inv = np.sqrt(2 * v1_pt * MET_pt * (1 - np.cos(delta_phi)))

    return M_inv