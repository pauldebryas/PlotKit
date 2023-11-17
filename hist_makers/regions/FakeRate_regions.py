from hist_makers.regions.helpers import Z_cut

def compute_region_mask(tree, channel, type):

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

      Loose_selection = no_samesign & met_cut & Zcut & mu1_iso & mu2_iso 

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

      Loose_selection = no_samesign & met_cut & Zcut & e1_iso & e2_iso 

  if channel == 'tem':
      #leptons
      mu_iso = tree["Muon_pfRelIso03_all"]  <= 0.15
      e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
      Tau_VSJet = tree[f"Tau_rawDeepTau2018v2p5VSjet"]  >= 0.9632

      #common cuts
      no_samesign = (tree["Muon_charge"]  != tree["Electron_charge"] ) | (tree["Muon_charge"]  != tree["Tau_charge"] )
      met_cut = tree["MET_pt"]  >= 25.
      no_bjets = tree["nbjets"]  == 0

      Loose_selection = no_samesign & met_cut & mu_iso & e_iso 

  ControlRegion = Loose_selection & ~no_bjets
  SignalRegion = Loose_selection  & no_bjets

  cut_region = {}
  cut_region['ControlRegion'] = ControlRegion
  # Fail tau_h ID
  cut_region['ControlRegion_fail'] = ControlRegion & ~Tau_VSJet
  cut_region['SignalRegion_fail'] = SignalRegion & ~Tau_VSJet
  # Pass tau_h ID
  cut_region['ControlRegion_pass'] = ControlRegion & Tau_VSJet
  #cut_region['SignalRegion_pass'] = SignalRegion & Tau_VSJet --> Blind

  if type == 'MC':
    #MC information: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched
    fake_tau =  (tree["Tau_genPartFlav"]  == 0)

    cut_region["ControlRegion_fake"] = ( ControlRegion & fake_tau)
    cut_region["ControlRegion_pass_true"] = ( ControlRegion & Tau_VSJet & ~fake_tau)
    cut_region["ControlRegion_pass_fake"] = ( ControlRegion & Tau_VSJet & fake_tau)
    cut_region["ControlRegion_fail_true"] = ( ControlRegion & ~Tau_VSJet & ~fake_tau)
    cut_region["ControlRegion_fail_fake"] = ( ControlRegion & ~Tau_VSJet & fake_tau)

  return cut_region
