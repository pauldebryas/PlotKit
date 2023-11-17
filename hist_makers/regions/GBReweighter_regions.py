from hist_makers.regions.helpers import Z_cut

def compute_region_mask(tree, channel):

  if channel == 'tmm':
      #leptons
      mu1_iso = tree["Muon1_pfRelIso03_all"] <= 0.15 
      mu2_iso = tree["Muon2_pfRelIso03_all"] <= 0.15
      Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"] >= 0.9632

      #common cuts
      no_samesign = (tree["Muon1_charge"] != tree["Tau_charge"]) | (tree["Muon1_charge"] != tree["Muon2_charge"])
      met_cut = tree["MET_pt"]  >= 25.
      Zcut = Z_cut(tree, 'Muon1', 'Muon2')
      no_bjets = tree["nbjets"] == 0

      Loose_selection = no_samesign & met_cut & mu1_iso & mu2_iso 

  if channel == 'tee':
      #leptons
      e1_iso = tree["Electron1_pfRelIso03_all"]  <= 0.15
      e2_iso = tree["Electron2_pfRelIso03_all"]  <= 0.15
      Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  >= 0.9632

      #common cuts
      no_samesign = (tree["Electron1_charge"]  != tree["Electron2_charge"] ) | (tree["Electron1_charge"]  != tree["Tau_charge"] )
      met_cut = tree["MET_pt"]  >= 25.
      Zcut = Z_cut(tree, 'Electron1', 'Electron2')
      no_bjets = tree["nbjets"]  == 0

      Loose_selection = no_samesign & met_cut & e1_iso & e2_iso 

  # control region: invert b jet requirement 
  # Ff determined in invert Z region (Sideband_Zveto)
  Sideband_Zveto = Loose_selection & ~Zcut
  ControlRegion = Loose_selection & Zcut & ~no_bjets
  SignalRegion = Loose_selection & Zcut & no_bjets

  cut_region = {}
  # Fail tau_h ID
  cut_region['SidebandZveto_fail'] = Sideband_Zveto & ~Tau_VSJet
  cut_region['ControlRegion_fail'] = ControlRegion & ~Tau_VSJet
  cut_region['SignalRegion_fail'] = SignalRegion & ~Tau_VSJet
  # Pass tau_h ID
  cut_region['SidebandZveto_pass'] = Sideband_Zveto & Tau_VSJet
  cut_region['ControlRegion_pass'] = ControlRegion & Tau_VSJet
  #cut_region['SignalRegion_pass'] = SignalRegion & Tau_VSJet --> Blind

  return cut_region

