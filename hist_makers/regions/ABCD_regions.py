from hist_makers.regions.helpers import Z_cut

def compute_region_cut(tree, channel):
  if channel == 'ttm':
    #leptons
    mu_iso = tree["Muon_pfRelIso03_all"] <= 0.15
    inv_mu_iso = tree["Muon_pfRelIso03_all"]  > 0.17

    Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    inv_Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"]  < 0.95
    inv_Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"]  < 0.95

    #common cuts
    no_bjets = tree["nbjets"]  == 0
    no_samesign = (tree["Muon_charge"]  != tree["Tau1_charge"] ) | (tree["Muon_charge"]  != tree["Tau2_charge"] )
    met_cut = tree["MET_pt"]  >= 25.
    Zcut = Z_cut(tree, 'Tau1', 'Tau2')

    cut_region = {}
    common_cut = no_bjets & no_samesign & met_cut & Zcut
    cut_region['A'] =  common_cut & inv_Tau1_VSJet & inv_Tau2_VSJet & inv_mu_iso
    cut_region['B'] =  common_cut & Tau1_VSJet     & Tau2_VSJet     & inv_mu_iso
    cut_region['C'] =  common_cut & inv_Tau1_VSJet & inv_Tau2_VSJet & mu_iso
    cut_region['D'] =  common_cut & Tau1_VSJet     & Tau2_VSJet     & mu_iso
    
  if channel == 'tte':
    #leptons
    e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
    Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"]  >= 0.9632
  
    inv_e_iso = tree["Electron_pfRelIso03_all"]  > 0.19
    inv_Tau1_VSJet = tree[f"Tau1_rawDeepTau2018v2p5VSjet"]  < 0.96
    inv_Tau2_VSJet = tree[f"Tau2_rawDeepTau2018v2p5VSjet"]  < 0.96

    #common cuts
    no_bjets = tree["nbjets"]  == 0
    #matchingHLTDoubleTau =  tree["matchingHLTDoubleTau"] 
    no_samesign = (tree["Electron_charge"]  != tree["Tau1_charge"] ) | (tree["Electron_charge"]  != tree["Tau2_charge"] )
    met_cut = tree["MET_pt"]  >= 25.
    Zcut = Z_cut(tree, 'Tau1', 'Tau2')

    cut_region = {}
    common_cut = no_bjets  & no_samesign & met_cut & Zcut #matchingHLTDoubleTau
    cut_region['A'] = common_cut & inv_e_iso & inv_Tau1_VSJet & inv_Tau2_VSJet
    cut_region['B'] = common_cut & inv_e_iso & Tau1_VSJet     & Tau2_VSJet   
    cut_region['C'] = common_cut & e_iso     & inv_Tau1_VSJet & inv_Tau2_VSJet  
    cut_region['D'] = common_cut & e_iso     & Tau1_VSJet     & Tau2_VSJet   

  if channel == 'tmm':
    #leptons
    mu1_iso = tree["Muon1_pfRelIso03_all"]  <= 0.15
    mu2_iso = tree["Muon2_pfRelIso03_all"]  <= 0.15
    inv_mu2_iso = tree["Muon2_pfRelIso03_all"]  > 0.17
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  < 0.96

    #common cuts
    no_bjets = tree["nbjets"]  == 0
    no_samesign = (tree["Muon1_charge"]  != tree["Tau_charge"] ) | (tree["Muon1_charge"]  != tree["Muon2_charge"] )
    met_cut = tree["MET_pt"]  >= 25.
    Zcut = Z_cut(tree, 'Muon1', 'Muon2')

    cut_region = {}
    common_cut = no_bjets & no_samesign & met_cut & Zcut & mu1_iso
    cut_region['A'] = common_cut & inv_mu2_iso & inv_Tau_VSJet 
    cut_region['B'] = common_cut & inv_mu2_iso & Tau_VSJet
    cut_region['C'] = common_cut & mu2_iso     & inv_Tau_VSJet 
    cut_region['D'] = common_cut & mu2_iso     & Tau_VSJet  

  if channel == 'tee':
    #leptons
    e1_iso = tree["Electron1_pfRelIso03_all"]  <= 0.15
    e2_iso = tree["Electron2_pfRelIso03_all"]  <= 0.15
    inv_e2_iso = tree["Electron2_pfRelIso03_all"]  > 0.18
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  < 0.93
    #common cuts
    no_bjets = tree["nbjets"]  == 0
    no_samesign = (tree["Electron1_charge"]  != tree["Electron2_charge"] ) | (tree["Electron1_charge"]  != tree["Tau_charge"] )
    met_cut = tree["MET_pt"]  >= 25.
    Zcut = Z_cut(tree, 'Electron1', 'Electron2')

    cut_region = {}
    common_cut = no_bjets & no_samesign & met_cut & Zcut & e1_iso
    cut_region['A'] = common_cut & inv_e2_iso & inv_Tau_VSJet
    cut_region['B'] = common_cut & inv_e2_iso & Tau_VSJet 
    cut_region['C'] = common_cut & e2_iso     & inv_Tau_VSJet
    cut_region['D'] = common_cut & e2_iso     & Tau_VSJet  

  if channel == 'tem':
    #leptons
    mu_iso = tree["Muon_pfRelIso03_all"]  <= 0.15
    inv_mu_iso = tree["Muon_pfRelIso03_all"]  > 0.19
    e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
    Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  >= 0.9632
    inv_Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  < 0.95

    #common cuts
    no_bjets = tree["nbjets"]  == 0
    no_samesign = (tree["Muon_charge"]  != tree["Electron_charge"] ) | (tree["Muon_charge"]  != tree["Tau_charge"] )
    met_cut = tree["MET_pt"]  >= 25.

    cut_region = {}
    common_cut = no_bjets & no_samesign & met_cut & e_iso 
    cut_region['A'] = common_cut & inv_mu_iso & inv_Tau_VSJet
    cut_region['B'] = common_cut & inv_mu_iso & Tau_VSJet 
    cut_region['C'] = common_cut & mu_iso     & inv_Tau_VSJet
    cut_region['D'] = common_cut & mu_iso     & Tau_VSJet 

  return cut_region
