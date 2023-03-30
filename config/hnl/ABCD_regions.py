import uproot
import numpy as np

def compute_region_cut(tree, channel):
  if channel == 'ttm':
    #leptons
    mu_iso = tree["Muon_pfRelIso03_all"].array() <= 0.15
    inv_mu_iso = tree["Muon_pfRelIso03_all"].array() > 0.17

    Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    inv_Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"].array() < 0.95
    inv_Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"].array() < 0.95

    #common cuts
    no_bjets = tree["nbjets"].array() == 0
    matchingIsoMu24 = tree["matchingIsoMu24"].array()
    no_samesign = (tree["Muon_charge"].array() != tree["Tau1_charge"].array()) | (tree["Muon_charge"].array() != tree["Tau2_charge"].array())
    met_cut = tree["MET_pt"].array() >= 25.

    cut_region = {}
    common_cut = no_bjets & matchingIsoMu24 & no_samesign & met_cut
    cut_region['A'] =  common_cut & inv_Tau1_VSJet & inv_Tau2_VSJet & inv_mu_iso
    cut_region['B'] =  common_cut & Tau1_VSJet     & Tau2_VSJet     & inv_mu_iso
    cut_region['C'] =  common_cut & inv_Tau1_VSJet & inv_Tau2_VSJet & mu_iso
    cut_region['D'] =  common_cut & Tau1_VSJet     & Tau2_VSJet     & mu_iso
    
  if channel == 'tte':
    #leptons
    e_iso = tree["Electron_pfRelIso03_all"].array() <= 0.15
    Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
  
    inv_e_iso = tree["Electron_pfRelIso03_all"].array() > 0.19
    inv_Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"].array() < 0.96
    inv_Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"].array() < 0.96

    #common cuts
    no_bjets = tree["nbjets"].array() == 0
    matchingHLTDoubleTau =  tree["matchingHLTDoubleTau"].array()
    no_samesign = (tree["Electron_charge"].array() != tree["Tau1_charge"].array()) | (tree["Electron_charge"].array() != tree["Tau2_charge"].array())
    met_cut = tree["MET_pt"].array() >= 25.

    cut_region = {}
    common_cut = no_bjets  & no_samesign & met_cut & matchingHLTDoubleTau
    cut_region['A'] = common_cut & inv_e_iso & inv_Tau1_VSJet & inv_Tau2_VSJet  
    cut_region['B'] = common_cut & inv_e_iso & Tau1_VSJet     & Tau2_VSJet   
    cut_region['C'] = common_cut & e_iso     & inv_Tau1_VSJet & inv_Tau2_VSJet  
    cut_region['D'] = common_cut & e_iso     & Tau1_VSJet     & Tau2_VSJet   

  if channel == 'tmm':
    #leptons
    mu1_iso = tree["Muon1_pfRelIso03_all"].array() <= 0.15
    mu2_iso = tree["Muon2_pfRelIso03_all"].array() <= 0.15
    inv_mu2_iso = tree["Muon2_pfRelIso03_all"].array() > 0.17
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() < 0.96

    #common cuts
    no_bjets = tree["nbjets"].array() == 0
    matchingIsoMu24 = tree["matchingIsoMu24"].array()
    no_samesign = (tree["Muon1_charge"].array() != tree["Tau_charge"].array()) | (tree["Muon1_charge"].array() != tree["Muon2_charge"].array())
    met_cut = tree["MET_pt"].array() >= 25.
    Zcut = Z_cut(tree, 'Muon1', 'Muon2')

    cut_region = {}
    common_cut = no_bjets & matchingIsoMu24 & no_samesign & met_cut & Zcut & mu1_iso
    cut_region['A'] = common_cut & inv_mu2_iso & inv_Tau_VSJet 
    cut_region['B'] = common_cut & inv_mu2_iso & Tau_VSJet
    cut_region['C'] = common_cut & mu2_iso     & inv_Tau_VSJet 
    cut_region['D'] = common_cut & mu2_iso     & Tau_VSJet  

  if channel == 'tee':
    #leptons
    e1_iso = tree["Electron1_pfRelIso03_all"].array() <= 0.15
    e2_iso = tree["Electron2_pfRelIso03_all"].array() <= 0.15
    inv_e2_iso = tree["Electron2_pfRelIso03_all"].array() > 0.18
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() < 0.93
    #common cuts
    no_bjets = tree["nbjets"].array() == 0
    no_samesign = (tree["Electron1_charge"].array() != tree["Electron2_charge"].array()) | (tree["Electron1_charge"].array() != tree["Tau_charge"].array())
    met_cut = tree["MET_pt"].array() >= 25.
    matchingEle32 = tree["matchingEle32"].array()
    Zcut = Z_cut(tree, 'Electron1', 'Electron2')

    cut_region = {}
    common_cut = no_bjets & no_samesign & met_cut & matchingEle32 & Zcut & e1_iso
    cut_region['A'] = common_cut & inv_e2_iso & inv_Tau_VSJet
    cut_region['B'] = common_cut & inv_e2_iso & Tau_VSJet 
    cut_region['C'] = common_cut & e2_iso     & inv_Tau_VSJet
    cut_region['D'] = common_cut & e2_iso     & Tau_VSJet  

  if channel == 'tem':
    #leptons
    mu_iso = tree["Muon_pfRelIso03_all"].array() <= 0.15
    inv_mu_iso = tree["Muon_pfRelIso03_all"].array() > 0.19
    e_iso = tree["Electron_pfRelIso03_all"].array() <= 0.15
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"].array() < 0.95

    #common cuts
    no_bjets = tree["nbjets"].array() == 0
    no_samesign = (tree["Muon_charge"].array() != tree["Electron_charge"].array()) | (tree["Muon_charge"].array() != tree["Tau_charge"].array())
    matchingIsoMu24 = tree["matchingIsoMu24"].array()
    met_cut = tree["MET_pt"].array() >= 25.

    cut_region = {}
    common_cut = no_bjets & no_samesign &  matchingIsoMu24 & met_cut & e_iso 
    cut_region['A'] = common_cut & inv_mu_iso & inv_Tau_VSJet
    cut_region['B'] = common_cut & inv_mu_iso & Tau_VSJet 
    cut_region['C'] = common_cut & mu_iso     & inv_Tau_VSJet
    cut_region['D'] = common_cut & mu_iso     & Tau_VSJet 

  return cut_region

def Z_cut(tree, Lepton1, Lepton2):
  mass_z = 91.2 #GeV
  interval = 10.
  ZCut = (mT_l1l2(tree, Lepton1, Lepton2) < (mass_z - interval)) |  (mT_l1l2(tree, Lepton1, Lepton2) > (mass_z + interval))
  return ZCut

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