from common.regions.helpers import Z_cut_OSSF, LowMass_cut_OSSF, MET_cut, mass_l1l2, Z_cut, delta_R, compute_invariant_mass_MET_lepton, Z_cut_l1l2l3

DeepTauVersion = 'idDeepTau2018v2p5VSjet' 
# Tau_genPartFlav:      Flavour of genParticle for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched
# Muon_genPartFlav:     Flavour of genParticle for MC matching to status==1 muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched
# Electron_genPartFlav: Flavour of genParticle for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched

def compute_region_mask(tree, channel, type, RegionName):

  cut_region = {}

  Tau_VSJet_cut = 6 #DeepTau2018v2p5VSjet ID working points (deepTau2018v2p5): 1 = VVVLoose, 2 = VVLoose, 3 = VLoose, 4 = Loose, 5 = Medium, 6 = Tight, 7 = VTight, 8 = VVTight

  if channel == 'tmm':
      #leptons ID
      mu1_iso = tree["Muon1_pfRelIso03_all"] <= 0.15 
      mu2_iso = tree["Muon2_pfRelIso03_all"] <= 0.15
      Tau_VSJet = tree[f"Tau_{DeepTauVersion}"] >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Muon1_charge"] != tree["Tau_charge"]) | (tree["Muon1_charge"] != tree["Muon2_charge"])
      met_cut = MET_cut(tree,  25.)
      Zcut = Z_cut_OSSF(tree, 'Muon1', 'Muon2')
      LowMasscut = LowMass_cut_OSSF(tree, 'Muon1', 'Muon2')
      no_bjets = tree["nbjetsLoose"] == 0
      OSlightLepton = tree["Muon1_charge"] != tree["Muon2_charge"]
      MET_lt50 = MET_cut(tree,  50., True)
      IMlightLepton_gt20 = mass_l1l2(tree, 'Muon1', 'Muon2') >= 20.
      MET_gt70 = MET_cut(tree,  70., False)

      if RegionName == 'SignalRegion':
        passLooseWP = no_samesign & met_cut & Zcut & LowMasscut & no_bjets 
        passTightWP = passLooseWP & Tau_VSJet & mu1_iso & mu2_iso
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & mu1_iso & mu2_iso)

      if RegionName == 'InvertedBjetsVetoRegion':
        passLooseWP = no_samesign & met_cut & Zcut & LowMasscut & ~no_bjets 
        passTightWP = passLooseWP & Tau_VSJet & mu1_iso & mu2_iso
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & mu1_iso & mu2_iso)

      if RegionName == 'ControlRegion':
        passLooseWP = no_samesign & LowMasscut & ~no_bjets 
        passTightWP = passLooseWP & Tau_VSJet & mu1_iso & mu2_iso
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & mu1_iso & mu2_iso)

      if RegionName == 'DYRegion':
        passLooseWP = mu1_iso & mu2_iso & OSlightLepton & ~Zcut & no_bjets
        passTightWP = passLooseWP & Tau_VSJet
        passLooseNotTightWP = passLooseWP & ~Tau_VSJet
 
      if RegionName == 'DYClosureRegion':
        passLooseWP = mu1_iso & mu2_iso & OSlightLepton & ~Zcut & ~MET_lt50
        passTightWP = passLooseWP & Tau_VSJet
        passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'ttClosureRegion':
        passLooseWP = mu1_iso & mu2_iso & OSlightLepton & Zcut & MET_gt70 & ~no_bjets
        passTightWP = passLooseWP & Tau_VSJet
        passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'ValidationRegionMuFR':
        passLooseWP = Zcut & ~no_bjets & Tau_VSJet 
        passTightWP = passLooseWP & mu1_iso & mu2_iso
        passLooseNotTightWP = passLooseWP & ~(mu1_iso & mu2_iso)
        cut_region[f'{RegionName}PFF'] = passLooseWP & ~mu1_iso & ~mu2_iso
        cut_region[f'{RegionName}PFF'] = passLooseWP & ~mu1_iso &  mu2_iso
        cut_region[f'{RegionName}PPF'] = passLooseWP &  mu1_iso & ~mu2_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
      if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
        cut_region[f'{RegionName}_AppRegionFFF'] = passLooseWP & ~Tau_VSJet & ~mu1_iso & ~mu2_iso
        cut_region[f'{RegionName}_AppRegionFFP'] = passLooseWP & ~Tau_VSJet & ~mu1_iso &  mu2_iso
        cut_region[f'{RegionName}_AppRegionFPF'] = passLooseWP & ~Tau_VSJet &  mu1_iso & ~mu2_iso
        cut_region[f'{RegionName}_AppRegionFPP'] = passLooseWP & ~Tau_VSJet &  mu1_iso &  mu2_iso
        cut_region[f'{RegionName}_AppRegionPFF'] = passLooseWP &  Tau_VSJet & ~mu1_iso & ~mu2_iso
        cut_region[f'{RegionName}_AppRegionPFP'] = passLooseWP &  Tau_VSJet & ~mu1_iso &  mu2_iso
        cut_region[f'{RegionName}_AppRegionPPF'] = passLooseWP &  Tau_VSJet &  mu1_iso & ~mu2_iso
        cut_region[f'{RegionName}_AppRegionPPP'] = passLooseWP &  Tau_VSJet &  mu1_iso &  mu2_iso 

      if type == 'MC':  
          TauIsPromptLepton =  (tree["Tau_genPartFlav"]  == 1) | (tree["Tau_genPartFlav"]  == 2)  | ((tree["Tau_genPartFlav"]  == 3) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 4) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 5) & (tree["Tau_isPrompt"] == 1))
          Muon1IsPromptLepton =  (tree["Muon1_genPartFlav"]  == 1) | (tree["Muon1_genPartFlav"]  == 15)
          Muon2IsPromptLepton =  (tree["Muon2_genPartFlav"]  == 1) | (tree["Muon2_genPartFlav"]  == 15)

          #used in XsecFile
          cut_region[f'{RegionName}_PassTightWP_TrueLeptons'] = passTightWP & TauIsPromptLepton & Muon1IsPromptLepton & Muon2IsPromptLepton

          #used for TauFRwightsDYttbar
          cut_region[f'{RegionName}_PassLooseNotTightWP_NotTrueLeptons'] = passLooseNotTightWP & ~(TauIsPromptLepton & Muon1IsPromptLepton & Muon2IsPromptLepton)

          #used in MCFakeFactor computation
          cut_region[f'{RegionName}_PassLooseWP_FakeTau'] = passLooseWP & ~TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau_passTight'] = passLooseWP & ~TauIsPromptLepton & Tau_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon1'] = passLooseWP & ~Muon1IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon1_passTight'] = passLooseWP & ~Muon1IsPromptLepton  & mu1_iso
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon2'] = passLooseWP & ~Muon2IsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon2_passTight'] = passLooseWP & ~Muon2IsPromptLepton  & mu2_iso

          #used in TauFakeFactorMC computation
          cut_region[f'{RegionName}_PassLooseWP_TauIsPromptLepton'] = passLooseWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_TauIsPromptLepton'] = passTightWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_TauIsPromptLepton'] = passLooseNotTightWP & TauIsPromptLepton

          if RegionName == 'ValidationRegionMuFR': 
            cut_region[f'{RegionName}_PassLooseWP_MuonsArePromptLepton'] = passLooseWP & Muon1IsPromptLepton & Muon2IsPromptLepton
            cut_region[f'{RegionName}_PassTightWP_MuonsArePromptLepton'] = passTightWP & Muon1IsPromptLepton & Muon2IsPromptLepton
            cut_region[f'{RegionName}_PassLooseNotTightWP_MuonsArePromptLepton'] = passLooseNotTightWP & Muon1IsPromptLepton & Muon2IsPromptLepton
            cut_region[f'{RegionName}PFF_TrueLepton'] = passLooseWP & ~mu1_iso & ~mu2_iso & (Muon1IsPromptLepton | Muon2IsPromptLepton)
            cut_region[f'{RegionName}PFF_TrueLepton'] = passLooseWP & ~mu1_iso &  mu2_iso & Muon1IsPromptLepton
            cut_region[f'{RegionName}PPF_TrueLepton'] = passLooseWP &  mu1_iso & ~mu2_iso & Muon2IsPromptLepton

          if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
            # cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~mu1_iso & ~mu2_iso & (TauIsPromptLepton | Muon1IsPromptLepton | Muon2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~mu1_iso &  mu2_iso & (TauIsPromptLepton | Muon1IsPromptLepton )
            # cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  mu1_iso & ~mu2_iso & (TauIsPromptLepton | Muon2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  mu1_iso &  mu2_iso &  TauIsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~mu1_iso & ~mu2_iso & (Muon1IsPromptLepton | Muon2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~mu1_iso &  mu2_iso & Muon1IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  mu1_iso & ~mu2_iso & Muon2IsPromptLepton

            cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~mu1_iso & ~mu2_iso & (TauIsPromptLepton & Muon1IsPromptLepton & Muon2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~mu1_iso &  mu2_iso & (TauIsPromptLepton & Muon1IsPromptLepton )
            cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  mu1_iso & ~mu2_iso & (TauIsPromptLepton & Muon2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  mu1_iso &  mu2_iso &  TauIsPromptLepton
            cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~mu1_iso & ~mu2_iso & (Muon1IsPromptLepton & Muon2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~mu1_iso &  mu2_iso & Muon1IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  mu1_iso & ~mu2_iso & Muon2IsPromptLepton

  if channel == 'tee':
      #leptons ID
      e1_iso = tree["Electron1_pfRelIso03_all"]  <= 0.15
      e2_iso = tree["Electron2_pfRelIso03_all"]  <= 0.15
      Tau_VSJet = tree[f"Tau_{DeepTauVersion}"]  >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Electron1_charge"]  != tree["Electron2_charge"] ) | (tree["Electron1_charge"]  != tree["Tau_charge"] )
      met_cut = MET_cut(tree,  25.)
      Zcut = Z_cut_OSSF(tree, 'Electron1', 'Electron2')
      LowMasscut = LowMass_cut_OSSF(tree, 'Electron1', 'Electron2')
      no_bjets = tree["nbjetsLoose"]  == 0
      OSlightLepton = tree["Electron1_charge"] != tree["Electron2_charge"]
      MET_lt50 = MET_cut(tree,  50., True)
      MET_gt70 = MET_cut(tree,  70., False)

      if RegionName == 'DYRegion':
          passLooseWP = e1_iso & e2_iso & OSlightLepton & ~Zcut & no_bjets
          passTightWP = passLooseWP & Tau_VSJet
          passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'DYClosureRegion':
          passLooseWP = e1_iso & e2_iso & OSlightLepton & ~Zcut & ~MET_lt50
          passTightWP = passLooseWP & Tau_VSJet
          passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'ttClosureRegion':
        passLooseWP = e1_iso & e2_iso & OSlightLepton & Zcut & MET_gt70 & ~no_bjets
        passTightWP = passLooseWP & Tau_VSJet
        passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'SignalRegion':
          passLooseWP = no_samesign & met_cut & Zcut & LowMasscut & no_bjets 
          passTightWP = passLooseWP & Tau_VSJet & e1_iso & e2_iso
          passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e1_iso & e2_iso)

      if RegionName == 'InvertedBjetsVetoRegion':
          passLooseWP = no_samesign & met_cut & Zcut & LowMasscut & ~no_bjets
          passTightWP = passLooseWP & Tau_VSJet & e1_iso & e2_iso
          passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e1_iso & e2_iso)

      if RegionName == 'ControlRegion':
          passLooseWP = no_samesign & LowMasscut & ~no_bjets
          passTightWP = passLooseWP & Tau_VSJet & e1_iso & e2_iso
          passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e1_iso & e2_iso)

      if RegionName == 'ValidationRegionEleFR':
        passLooseWP = Zcut & ~no_bjets & Tau_VSJet 
        passTightWP = passLooseWP & e1_iso & e2_iso
        passLooseNotTightWP = passLooseWP & ~(e1_iso & e2_iso)
        cut_region[f'{RegionName}PFF'] = passLooseWP & ~e1_iso & ~e2_iso
        cut_region[f'{RegionName}PFF'] = passLooseWP & ~e1_iso &  e2_iso
        cut_region[f'{RegionName}PPF'] = passLooseWP &  e1_iso & ~e2_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
      if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
        cut_region[f'{RegionName}_AppRegionFFF'] = passLooseWP & ~Tau_VSJet & ~e1_iso & ~e2_iso
        cut_region[f'{RegionName}_AppRegionFFP'] = passLooseWP & ~Tau_VSJet & ~e1_iso &  e2_iso
        cut_region[f'{RegionName}_AppRegionFPF'] = passLooseWP & ~Tau_VSJet &  e1_iso & ~e2_iso
        cut_region[f'{RegionName}_AppRegionFPP'] = passLooseWP & ~Tau_VSJet &  e1_iso &  e2_iso
        cut_region[f'{RegionName}_AppRegionPFF'] = passLooseWP &  Tau_VSJet & ~e1_iso & ~e2_iso
        cut_region[f'{RegionName}_AppRegionPFP'] = passLooseWP &  Tau_VSJet & ~e1_iso &  e2_iso
        cut_region[f'{RegionName}_AppRegionPPF'] = passLooseWP &  Tau_VSJet &  e1_iso & ~e2_iso
        cut_region[f'{RegionName}_AppRegionPPP'] = passLooseWP &  Tau_VSJet &  e1_iso &  e2_iso 

      if type == 'MC':  
          TauIsPromptLepton =  (tree["Tau_genPartFlav"]  == 1) | (tree["Tau_genPartFlav"]  == 2)  | ((tree["Tau_genPartFlav"]  == 3) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 4) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 5) & (tree["Tau_isPrompt"] == 1))
          Electron1IsPromptLepton =  (tree["Electron1_genPartFlav"]  == 1) | (tree["Electron1_genPartFlav"]  == 15)
          Electron2IsPromptLepton =  (tree["Electron2_genPartFlav"]  == 1) | (tree["Electron2_genPartFlav"]  == 15)

          cut_region[f'{RegionName}_PassTightWP_TrueLeptons'] = passTightWP & TauIsPromptLepton & Electron1IsPromptLepton & Electron2IsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_NotTrueLeptons'] = passLooseNotTightWP & ~(TauIsPromptLepton & Electron1IsPromptLepton & Electron2IsPromptLepton)

          cut_region[f'{RegionName}_PassLooseWP_FakeTau'] = passLooseWP & ~TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau_passTight'] = passLooseWP & ~TauIsPromptLepton & Tau_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron1'] = passLooseWP & ~Electron1IsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron1_passTight'] = passLooseWP & ~Electron1IsPromptLepton & e1_iso
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron2'] = passLooseWP & ~Electron2IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron2_passTight'] = passLooseWP & ~Electron2IsPromptLepton & e2_iso

          cut_region[f'{RegionName}_PassLooseWP_TauIsPromptLepton'] = passLooseWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_TauIsPromptLepton'] = passTightWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_TauIsPromptLepton'] = passLooseNotTightWP & TauIsPromptLepton

          if RegionName == 'ValidationRegionEleFR': 
            cut_region[f'{RegionName}_PassLooseWP_ElectronsArePromptLepton'] = passLooseWP & Electron1IsPromptLepton & Electron2IsPromptLepton
            cut_region[f'{RegionName}_PassTightWP_ElectronsArePromptLepton'] = passTightWP & Electron1IsPromptLepton & Electron2IsPromptLepton
            cut_region[f'{RegionName}_PassLooseNotTightWP_ElectronsArePromptLepton'] = passLooseNotTightWP & Electron1IsPromptLepton & Electron2IsPromptLepton
            cut_region[f'{RegionName}PFF_TrueLepton'] = passLooseWP & ~e1_iso & ~e2_iso & (Electron1IsPromptLepton | Electron2IsPromptLepton)
            cut_region[f'{RegionName}PFP_TrueLepton'] = passLooseWP & ~e1_iso &  e2_iso & Electron1IsPromptLepton
            cut_region[f'{RegionName}PPF_TrueLepton'] = passLooseWP &  e1_iso & ~e2_iso & Electron2IsPromptLepton

          if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
            # cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e1_iso & ~e2_iso & (TauIsPromptLepton | Electron1IsPromptLepton | Electron2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e1_iso &  e2_iso & (TauIsPromptLepton | Electron1IsPromptLepton )
            # cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e1_iso & ~e2_iso & (TauIsPromptLepton | Electron2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e1_iso &  e2_iso &  TauIsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e1_iso & ~e2_iso & (Electron1IsPromptLepton | Electron2IsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e1_iso &  e2_iso & Electron1IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  e1_iso & ~e2_iso & Electron2IsPromptLepton

            cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e1_iso & ~e2_iso & (TauIsPromptLepton & Electron1IsPromptLepton & Electron2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e1_iso &  e2_iso & (TauIsPromptLepton & Electron1IsPromptLepton )
            cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e1_iso & ~e2_iso & (TauIsPromptLepton & Electron2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e1_iso &  e2_iso &  TauIsPromptLepton
            cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e1_iso & ~e2_iso & (Electron1IsPromptLepton & Electron2IsPromptLepton)
            cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e1_iso &  e2_iso & Electron1IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  e1_iso & ~e2_iso & Electron2IsPromptLepton

            
  if channel == 'tem':
      #leptons ID
      mu_iso = tree["Muon_pfRelIso03_all"]  <= 0.15
      e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
      Tau_VSJet = tree[f"Tau_{DeepTauVersion}"]  >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Muon_charge"]  != tree["Electron_charge"] ) | (tree["Muon_charge"]  != tree["Tau_charge"] )
      met_cut = MET_cut(tree,  25.)
      no_bjets = tree["nbjetsLoose"]  == 0
      OSlightLepton = tree["Electron_charge"] != tree["Muon_charge"]
      IMlightLepton_gt20 = mass_l1l2(tree, 'Electron', 'Muon') >= 20.
      
      if RegionName == 'ttbarRegionFR':
        passLooseWP = mu_iso & e_iso & OSlightLepton & ~no_bjets
        passTightWP = passLooseWP & Tau_VSJet
        passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      # if RegionName == 'ttbarRegionValidationEleFR':
      #   passLooseWP = Tau_VSJet & mu_iso & OSlightLepton & ~no_bjets 
      #   passTightWP = passLooseWP & e_iso
      #   passLooseNotTightWP = passLooseWP & ~e_iso

      # if RegionName == 'ttValidationRegion':
      #   passLooseWP = mu_iso & e_iso & OSlightLepton & ~IMlightLepton_gt20 & ~no_bjets
      #   passTightWP = passLooseWP & Tau_VSJet
      #   passLooseNotTightWP = passLooseWP & ~Tau_VSJet

      if RegionName == 'SignalRegion':
        passLooseWP = no_samesign & met_cut & no_bjets
        passTightWP = passLooseWP & Tau_VSJet & e_iso & mu_iso 
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e_iso & mu_iso)
      
      if RegionName == 'InvertedBjetsVetoRegion':
        passLooseWP = no_samesign & met_cut & ~no_bjets
        passTightWP = passLooseWP & Tau_VSJet & e_iso & mu_iso 
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e_iso & mu_iso)

      if RegionName == 'ControlRegion':
        passLooseWP = no_samesign & ~no_bjets
        passTightWP = passLooseWP & Tau_VSJet & e_iso & mu_iso 
        passLooseNotTightWP = passLooseWP & ~(Tau_VSJet & e_iso & mu_iso)

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
      if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
        cut_region[f'{RegionName}_AppRegionFFF'] = passLooseWP & ~Tau_VSJet & ~e_iso & ~mu_iso
        cut_region[f'{RegionName}_AppRegionFFP'] = passLooseWP & ~Tau_VSJet & ~e_iso &  mu_iso
        cut_region[f'{RegionName}_AppRegionFPF'] = passLooseWP & ~Tau_VSJet &  e_iso & ~mu_iso
        cut_region[f'{RegionName}_AppRegionFPP'] = passLooseWP & ~Tau_VSJet &  e_iso &  mu_iso
        cut_region[f'{RegionName}_AppRegionPFF'] = passLooseWP &  Tau_VSJet & ~e_iso & ~mu_iso
        cut_region[f'{RegionName}_AppRegionPFP'] = passLooseWP &  Tau_VSJet & ~e_iso &  mu_iso
        cut_region[f'{RegionName}_AppRegionPPF'] = passLooseWP &  Tau_VSJet &  e_iso & ~mu_iso
        cut_region[f'{RegionName}_AppRegionPPP'] = passLooseWP &  Tau_VSJet &  e_iso &  mu_iso 
    
      if type == 'MC':  
          TauIsPromptLepton =  (tree["Tau_genPartFlav"]  == 1) | (tree["Tau_genPartFlav"]  == 2)  | ((tree["Tau_genPartFlav"]  == 3) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 4) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 5) & (tree["Tau_isPrompt"] == 1))
          ElectronIsPromptLepton = (tree["Electron_genPartFlav"]  == 1) | (tree["Electron_genPartFlav"]  == 15)
          MuonIsPromptLepton =     (tree["Muon_genPartFlav"]  == 1) | (tree["Muon_genPartFlav"]  == 15)
        
          cut_region[f'{RegionName}_PassTightWP_TrueLeptons'] = passTightWP & TauIsPromptLepton & ElectronIsPromptLepton & MuonIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_NotTrueLeptons'] = passLooseNotTightWP & ~(TauIsPromptLepton & ElectronIsPromptLepton & MuonIsPromptLepton)

          cut_region[f'{RegionName}_PassLooseWP_FakeTau'] = passLooseWP & ~TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau_passTight'] = passLooseWP & ~TauIsPromptLepton & Tau_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron'] = passLooseWP & ~ElectronIsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron_passTight'] = passLooseWP & ~ElectronIsPromptLepton & e_iso
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon'] = passLooseWP & ~MuonIsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon_passTight'] = passLooseWP & ~MuonIsPromptLepton  & mu_iso

          cut_region[f'{RegionName}_PassLooseWP_TauIsPromptLepton'] = passLooseWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_TauIsPromptLepton'] = passTightWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_TauIsPromptLepton'] = passLooseNotTightWP & TauIsPromptLepton

          cut_region[f'{RegionName}_PassLooseWP_ElectronIsPromptLepton'] = passLooseWP & ElectronIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_ElectronIsPromptLepton'] = passTightWP & ElectronIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_ElectronIsPromptLepton'] = passLooseNotTightWP & ElectronIsPromptLepton

          cut_region[f'{RegionName}_PassLooseWP_MuonIsPromptLepton'] = passLooseWP & MuonIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_MuonIsPromptLepton'] = passTightWP & MuonIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_MuonIsPromptLepton'] = passLooseNotTightWP & MuonIsPromptLepton

          if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
            # cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e_iso & ~mu_iso & (TauIsPromptLepton | ElectronIsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e_iso &  mu_iso & (TauIsPromptLepton | ElectronIsPromptLepton )
            # cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e_iso & ~mu_iso & (TauIsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e_iso &  mu_iso &  TauIsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e_iso & ~mu_iso & (ElectronIsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e_iso &  mu_iso & ElectronIsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  e_iso & ~mu_iso & MuonIsPromptLepton

            cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e_iso & ~mu_iso & (TauIsPromptLepton & ElectronIsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau_VSJet & ~e_iso &  mu_iso & (TauIsPromptLepton & ElectronIsPromptLepton )
            cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e_iso & ~mu_iso & (TauIsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau_VSJet &  e_iso &  mu_iso &  TauIsPromptLepton
            cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e_iso & ~mu_iso & (ElectronIsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau_VSJet & ~e_iso &  mu_iso & ElectronIsPromptLepton
            cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau_VSJet &  e_iso & ~mu_iso & MuonIsPromptLepton

  if channel == 'tte':
      #leptons ID
      e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
      Tau1_VSJet = tree[f"Tau1_{DeepTauVersion}"]  >= Tau_VSJet_cut
      Tau2_VSJet = tree[f"Tau2_{DeepTauVersion}"]  >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Electron_charge"]  != tree["Tau1_charge"] ) | (tree["Electron_charge"]  != tree["Tau2_charge"] )
      met_cut = MET_cut(tree,  25.)
      Zcut = Z_cut_OSSF(tree, 'Tau1', 'Tau2')
      LowMasscut = LowMass_cut_OSSF(tree, 'Tau1', 'Tau2')
      no_bjets = tree["nbjetsLoose"]  == 0

      if RegionName == 'SignalRegion':
        passLooseWP = no_samesign & met_cut & Zcut  & LowMasscut & no_bjets 
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & e_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & e_iso)

      if RegionName == 'InvertedBjetsVetoRegion':
        passLooseWP = no_samesign & met_cut & Zcut  & LowMasscut & ~no_bjets 
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & e_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & e_iso)

      if RegionName == 'ControlRegion':
        passLooseWP = no_samesign  & LowMasscut & ~no_bjets 
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & e_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & e_iso)

      if RegionName == 'ValidationRegionEleFR':
        passLooseWP = ~Zcut & Tau1_VSJet & Tau2_VSJet 
        passTightWP = passLooseWP & e_iso 
        passLooseNotTightWP = passLooseWP & ~e_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
      if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
        cut_region[f'{RegionName}_AppRegionFFF'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~e_iso
        cut_region[f'{RegionName}_AppRegionFFP'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  e_iso
        cut_region[f'{RegionName}_AppRegionFPF'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~e_iso
        cut_region[f'{RegionName}_AppRegionFPP'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  e_iso
        cut_region[f'{RegionName}_AppRegionPFF'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~e_iso
        cut_region[f'{RegionName}_AppRegionPFP'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  e_iso
        cut_region[f'{RegionName}_AppRegionPPF'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~e_iso
        cut_region[f'{RegionName}_AppRegionPPP'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet &  e_iso 

      if type == 'MC':  
          Tau1IsPromptLepton  =    (tree["Tau1_genPartFlav"]  == 1) | (tree["Tau1_genPartFlav"]  == 2)  | ((tree["Tau1_genPartFlav"]  == 3) & (tree["Tau1_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau1_genPartFlav"]  == 4) & (tree["Tau1_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau1_genPartFlav"]  == 5) & (tree["Tau1_isPrompt"] == 1))
          Tau2IsPromptLepton  =    (tree["Tau2_genPartFlav"]  == 1) | (tree["Tau2_genPartFlav"]  == 2)  | ((tree["Tau2_genPartFlav"]  == 3) & (tree["Tau2_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau2_genPartFlav"]  == 4) & (tree["Tau2_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau2_genPartFlav"]  == 5) & (tree["Tau2_isPrompt"] == 1))
          ElectronIsPromptLepton = (tree["Electron_genPartFlav"]  == 1) | (tree["Electron_genPartFlav"]  == 15)

          cut_region[f'{RegionName}_PassTightWP_TrueLeptons'] = passTightWP & Tau1IsPromptLepton & Tau2IsPromptLepton & ElectronIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_NotTrueLeptons'] = passLooseNotTightWP & ~(Tau1IsPromptLepton & Tau2IsPromptLepton & ElectronIsPromptLepton)

          cut_region[f'{RegionName}_PassLooseWP_FakeTau1'] = passLooseWP & ~Tau1IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau1_passTight'] = passLooseWP & ~Tau1IsPromptLepton & Tau1_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeTau2'] = passLooseWP & ~Tau2IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau2_passTight'] = passLooseWP & ~Tau2IsPromptLepton & Tau2_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron'] = passLooseWP & ~ElectronIsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeElectron_passTight'] = passLooseWP & ~ElectronIsPromptLepton & e_iso

          if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
            # cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~e_iso & (Tau1IsPromptLepton | Tau2IsPromptLepton | ElectronIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  e_iso & (Tau1IsPromptLepton | Tau2IsPromptLepton )
            # cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~e_iso & (Tau1IsPromptLepton | ElectronIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  e_iso &  Tau1IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~e_iso & (Tau2IsPromptLepton | ElectronIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  e_iso & Tau2IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~e_iso & ElectronIsPromptLepton

            cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~e_iso & (Tau1IsPromptLepton & Tau2IsPromptLepton & ElectronIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  e_iso & (Tau1IsPromptLepton & Tau2IsPromptLepton )
            cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~e_iso & (Tau1IsPromptLepton & ElectronIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  e_iso &  Tau1IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~e_iso & (Tau2IsPromptLepton & ElectronIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  e_iso & Tau2IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~e_iso & ElectronIsPromptLepton
      
  if channel == 'ttm':
      #leptons ID
      mu_iso = tree["Muon_pfRelIso03_all"] <= 0.15
      Tau1_VSJet = tree[f"Tau1_{DeepTauVersion}"]  >= Tau_VSJet_cut
      Tau2_VSJet = tree[f"Tau2_{DeepTauVersion}"]  >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Muon_charge"]  != tree["Tau1_charge"] ) | (tree["Muon_charge"]  != tree["Tau2_charge"] )
      met_cut = MET_cut(tree,  25.)
      Zcut = Z_cut_OSSF(tree, 'Tau1', 'Tau2')
      LowMasscut = LowMass_cut_OSSF(tree, 'Tau1', 'Tau2')
      no_bjets = tree["nbjetsLoose"]  == 0

      if RegionName == 'SignalRegion':
        passLooseWP = no_samesign & met_cut & Zcut  & LowMasscut & no_bjets 
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & mu_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & mu_iso)

      if RegionName == 'InvertedBjetsVetoRegion':
        passLooseWP = no_samesign & met_cut & Zcut  & LowMasscut & ~no_bjets
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & mu_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & mu_iso)

      if RegionName == 'ControlRegion':
        passLooseWP = no_samesign  & LowMasscut & ~no_bjets
        passTightWP = passLooseWP & Tau1_VSJet & Tau2_VSJet & mu_iso
        passLooseNotTightWP = passLooseWP & ~(Tau1_VSJet & Tau2_VSJet & mu_iso)

      if RegionName == 'ValidationRegionMuFR':
        passLooseWP = ~Zcut & Tau1_VSJet & Tau2_VSJet
        passTightWP = passLooseWP & mu_iso 
        passLooseNotTightWP = passLooseWP & ~mu_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
      if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
        cut_region[f'{RegionName}_AppRegionFFF'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~mu_iso
        cut_region[f'{RegionName}_AppRegionFFP'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  mu_iso
        cut_region[f'{RegionName}_AppRegionFPF'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~mu_iso
        cut_region[f'{RegionName}_AppRegionFPP'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  mu_iso
        cut_region[f'{RegionName}_AppRegionPFF'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~mu_iso
        cut_region[f'{RegionName}_AppRegionPFP'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  mu_iso
        cut_region[f'{RegionName}_AppRegionPPF'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~mu_iso
        cut_region[f'{RegionName}_AppRegionPPP'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet &  mu_iso 

      if type == 'MC':  
          Tau1IsPromptLepton =  (tree["Tau1_genPartFlav"]  == 1) | (tree["Tau1_genPartFlav"]  == 2)  | ((tree["Tau1_genPartFlav"]  == 3) & (tree["Tau1_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau1_genPartFlav"]  == 4) & (tree["Tau1_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau1_genPartFlav"]  == 5) & (tree["Tau1_isPrompt"] == 1))
          Tau2IsPromptLepton =  (tree["Tau2_genPartFlav"]  == 1) | (tree["Tau2_genPartFlav"]  == 2)  | ((tree["Tau2_genPartFlav"]  == 3) & (tree["Tau2_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau2_genPartFlav"]  == 4) & (tree["Tau2_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau2_genPartFlav"]  == 5) & (tree["Tau2_isPrompt"] == 1))
          MuonIsPromptLepton =   (tree["Muon_genPartFlav"]  == 1) | (tree["Muon_genPartFlav"]  == 15)
          
          cut_region[f'{RegionName}_PassTightWP_TrueLeptons'] = passTightWP & Tau1IsPromptLepton & Tau2IsPromptLepton & MuonIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_NotTrueLeptons'] = passLooseNotTightWP & ~(Tau1IsPromptLepton & Tau2IsPromptLepton & MuonIsPromptLepton)

          cut_region[f'{RegionName}_PassLooseWP_FakeTau1'] = passLooseWP & ~Tau1IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau1_passTight'] = passLooseWP & ~Tau1IsPromptLepton & Tau1_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeTau2'] = passLooseWP & ~Tau2IsPromptLepton
          cut_region[f'{RegionName}_PassLooseWP_FakeTau2_passTight'] = passLooseWP & ~Tau2IsPromptLepton & Tau2_VSJet
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon'] = passLooseWP & ~MuonIsPromptLepton 
          cut_region[f'{RegionName}_PassLooseWP_FakeMuon_passTight'] = passLooseWP & ~MuonIsPromptLepton  & mu_iso

          if RegionName == 'ValidationRegionMuFR':
            cut_region[f'{RegionName}_PassLooseWP_MuonIsPromptLepton'] = passLooseWP & MuonIsPromptLepton
            cut_region[f'{RegionName}_PassTightWP_MuonIsPromptLepton'] = passTightWP & MuonIsPromptLepton
            cut_region[f'{RegionName}_PassLooseNotTightWP_MuonIsPromptLepton'] = passLooseNotTightWP & MuonIsPromptLepton

          if RegionName in ['SignalRegion', 'InvertedBjetsVetoRegion', 'ControlRegion']:
            # cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~mu_iso & (Tau1IsPromptLepton | Tau2IsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  mu_iso & (Tau1IsPromptLepton | Tau2IsPromptLepton )
            # cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~mu_iso & (Tau1IsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  mu_iso &  Tau1IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~mu_iso & (Tau2IsPromptLepton | MuonIsPromptLepton)
            # cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  mu_iso & Tau2IsPromptLepton
            # cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~mu_iso & MuonIsPromptLepton

            cut_region[f'{RegionName}_AppRegionFFF_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet & ~mu_iso & (Tau1IsPromptLepton & Tau2IsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFFP_TrueLepton'] = passLooseWP & ~Tau1_VSJet & ~Tau2_VSJet &  mu_iso & (Tau1IsPromptLepton & Tau2IsPromptLepton )
            cut_region[f'{RegionName}_AppRegionFPF_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet & ~mu_iso & (Tau1IsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionFPP_TrueLepton'] = passLooseWP & ~Tau1_VSJet &  Tau2_VSJet &  mu_iso &  Tau1IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPFF_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet & ~mu_iso & (Tau2IsPromptLepton & MuonIsPromptLepton)
            cut_region[f'{RegionName}_AppRegionPFP_TrueLepton'] = passLooseWP &  Tau1_VSJet & ~Tau2_VSJet &  mu_iso & Tau2IsPromptLepton
            cut_region[f'{RegionName}_AppRegionPPF_TrueLepton'] = passLooseWP &  Tau1_VSJet &  Tau2_VSJet & ~mu_iso & MuonIsPromptLepton
  
  if channel == 'Zmu':
      #leptons ID
      mu_iso = tree["Muon_pfRelIso03_all"] <= 0.15
      lep1_iso = tree["Lepton1_pfRelIso03_all"] <= 0.15
      lep2_iso = tree["Lepton2_pfRelIso03_all"] <= 0.15
      dr_l1_l2 = delta_R(tree, 'Lepton1', 'Lepton2') > 0.5
      m_MET_mu_cut = compute_invariant_mass_MET_lepton(tree, 'Muon') <= 50
      m_l1_l2_l3_cut = Z_cut_l1l2l3(tree, 'Lepton1', 'Lepton2', 'Muon', interval = 10.)
      no_bjets = tree["nbjetsLoose"] == 0
      met_cut = MET_cut(tree,  50., lowerThan = True)  
      #pairEvents = tree["event"]%2 == 0

      if RegionName == 'All':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) 

      # if RegionName == 'DYRegion':
      #   LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_mu_cut & m_l1_l2_l3_cut 

      if RegionName == 'DYRegionFR':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_mu_cut & m_l1_l2_l3_cut & no_bjets

      # if RegionName == 'DYRegionValidationBis':
      #   LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_mu_cut & m_l1_l2_l3_cut & ~no_bjets

      if RegionName == 'DYRegionValidation':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & met_cut & ~no_bjets

      passLooseWP = LooseCut
      passTightWP = passLooseWP & mu_iso 
      passLooseNotTightWP = passLooseWP & ~mu_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP

      if type == 'MC':
        MuonIsPromptLepton =   (tree["Muon_genPartFlav"]  == 1) | (tree["Muon_genPartFlav"]  == 15)
        cut_region[f'{RegionName}_PassLooseWP_MuonIsPromptLepton'] = passLooseWP & MuonIsPromptLepton
        cut_region[f'{RegionName}_PassTightWP_MuonIsPromptLepton'] = passTightWP & MuonIsPromptLepton
        cut_region[f'{RegionName}_PassLooseNotTightWP_MuonIsPromptLepton'] = passLooseNotTightWP & MuonIsPromptLepton  

  if channel == 'Ze':
      #leptons ID
      e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
      lep1_iso = tree["Lepton1_pfRelIso03_all"] <= 0.15
      lep2_iso = tree["Lepton2_pfRelIso03_all"] <= 0.15
      dr_l1_l2 = delta_R(tree, 'Lepton1', 'Lepton2') > 0.5
      m_MET_e_cut = compute_invariant_mass_MET_lepton(tree, 'Electron') <= 50
      m_l1_l2_l3_cut = Z_cut_l1l2l3(tree, 'Lepton1', 'Lepton2', 'Electron', interval = 10.)
      no_bjets = tree["nbjetsLoose"] == 0
      met_cut = MET_cut(tree,  50., lowerThan = True)  
      #pairEvents = tree["event"]%2 == 0

      if RegionName == 'All':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) 

      # if RegionName == 'DYRegion':
      #   LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_e_cut & m_l1_l2_l3_cut

      if RegionName == 'DYRegionFR':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_e_cut & m_l1_l2_l3_cut & no_bjets

      # if RegionName == 'DYRegionValidationBis':
      #   LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & m_MET_e_cut & m_l1_l2_l3_cut & ~no_bjets

      if RegionName == 'DYRegionValidation':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & met_cut & ~no_bjets

      passLooseWP = LooseCut
      passTightWP = passLooseWP & e_iso 
      passLooseNotTightWP = passLooseWP & ~e_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP

      if type == 'MC':
        ElectronIsPromptLepton =   (tree["Electron_genPartFlav"]  == 1) | (tree["Electron_genPartFlav"]  == 15)
        cut_region[f'{RegionName}_PassLooseWP_ElectronIsPromptLepton'] = passLooseWP & ElectronIsPromptLepton
        cut_region[f'{RegionName}_PassTightWP_ElectronIsPromptLepton'] = passTightWP & ElectronIsPromptLepton
        cut_region[f'{RegionName}_PassLooseNotTightWP_ElectronIsPromptLepton'] = passLooseNotTightWP & ElectronIsPromptLepton

  if channel == 'tll':
      #leptons ID
      l1_iso = tree["Lepton1_pfRelIso03_all"] <= 0.15 
      l2_iso = tree["Lepton2_pfRelIso03_all"] <= 0.15
      Tau_VSJet = tree[f"Tau_{DeepTauVersion}"] >= Tau_VSJet_cut

      #common cuts
      no_samesign = (tree["Lepton1_charge"] != tree["Tau_charge"]) | (tree["Lepton1_charge"] != tree["Lepton2_charge"])
      met_cut = MET_cut(tree,  25.)
      Zcut = Z_cut_OSSF(tree, 'Lepton1', 'Lepton2')
      no_bjets = tree["nbjetsLoose"] == 0
      OSlightLepton = tree["Lepton1_charge"] != tree["Lepton2_charge"]
      MET_lt50 = MET_cut(tree,  50., True)
      MET_gt70 = MET_cut(tree,  70., False)

      if RegionName == 'DYRegionFR': 
        LooseCut = l1_iso & l2_iso & OSlightLepton & ~Zcut & no_bjets

      if RegionName == 'DYRegionValidation':
        LooseCut = l1_iso & l2_iso & OSlightLepton & ~Zcut & ~no_bjets 

      if RegionName == 'ttbarRegionValidation':
        LooseCut = l1_iso & l2_iso & OSlightLepton & Zcut & ~no_bjets & MET_gt70

      # if RegionName == 'ttbarRegionFR':
      #   LooseCut = l1_iso & l2_iso & OSlightLepton & Zcut & ~no_bjets & ~MET_gt70

      # if RegionName == 'ttbarRegionValidation':
      #   LooseCut = l1_iso & l2_iso & OSlightLepton & Zcut & ~no_bjets & MET_gt70

      passLooseWP = LooseCut
      passTightWP = passLooseWP & Tau_VSJet
      passLooseNotTightWP = passLooseWP & ~Tau_VSJet
      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP
  
      if type == 'MC':  
          TauIsPromptLepton =  (tree["Tau_genPartFlav"]  == 1) | (tree["Tau_genPartFlav"]  == 2)  | ((tree["Tau_genPartFlav"]  == 3) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 4) & (tree["Tau_isDirectPromptTauDecayProduct"] == 1)) | ((tree["Tau_genPartFlav"]  == 5) & (tree["Tau_isPrompt"] == 1))

          #used in TauFakeFactorMC computation
          cut_region[f'{RegionName}_PassLooseWP_TauIsPromptLepton'] = passLooseWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassTightWP_TauIsPromptLepton'] = passTightWP & TauIsPromptLepton
          cut_region[f'{RegionName}_PassLooseNotTightWP_TauIsPromptLepton'] = passLooseNotTightWP & TauIsPromptLepton

  if channel == 'llmu':
      #leptons ID
      mu_iso = tree["Muon_pfRelIso03_all"] <= 0.15
      lep1_iso = tree["Lepton1_pfRelIso03_all"] <= 0.15
      lep2_iso = tree["Lepton2_pfRelIso03_all"] <= 0.15
      dr_l1_l2 = delta_R(tree, 'Lepton1', 'Lepton2') > 0.5
      no_bjets = tree["nbjetsLoose"] == 0
      MET_gt50 = MET_cut(tree,  50., False)
      Z_cut_l2Muon = Z_cut(tree, 'Lepton2', 'Muon', interval=10.)
      Z_cut_l1l2Muon = Z_cut_l1l2l3(tree, 'Lepton1', 'Lepton2', 'Muon', interval = 10.)

      if RegionName == 'All':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4)

      if RegionName == 'ttbarRegionFR':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & ~no_bjets & Z_cut_l2Muon & MET_gt50

      if RegionName == 'ttbarRegionValidation':
        LooseCut   = (tree["Muon_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & ~no_bjets & Z_cut_l2Muon & ~MET_gt50 & Z_cut_l1l2Muon

      passLooseWP = LooseCut
      passTightWP = passLooseWP & mu_iso 
      passLooseNotTightWP = passLooseWP & ~mu_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP

      if type == 'MC':
        MuonIsPromptLepton =   (tree["Muon_genPartFlav"]  == 1) | (tree["Muon_genPartFlav"]  == 15)
        cut_region[f'{RegionName}_PassLooseWP_MuonIsPromptLepton'] = passLooseWP & MuonIsPromptLepton
        cut_region[f'{RegionName}_PassTightWP_MuonIsPromptLepton'] = passTightWP & MuonIsPromptLepton
        cut_region[f'{RegionName}_PassLooseNotTightWP_MuonIsPromptLepton'] = passLooseNotTightWP & MuonIsPromptLepton 

  if channel == 'lle':
      #leptons ID
      e_iso = tree["Electron_pfRelIso03_all"]  <= 0.15
      lep1_iso = tree["Lepton1_pfRelIso03_all"] <= 0.15
      lep2_iso = tree["Lepton2_pfRelIso03_all"] <= 0.15
      dr_l1_l2 = delta_R(tree, 'Lepton1', 'Lepton2') > 0.5
      no_bjets = tree["nbjetsLoose"] == 0
      MET_gt50 = MET_cut(tree,  50., False)
      Z_cut_l1Electron = Z_cut(tree, 'Lepton1', 'Electron', interval=10.)
      Z_cut_l1l2Electron = Z_cut_l1l2l3(tree, 'Lepton1', 'Lepton2', 'Electron', interval = 10.)

      if RegionName == 'All':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4)

      if RegionName == 'ttbarRegionFR':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & ~no_bjets & Z_cut_l1Electron & MET_gt50

      if RegionName == 'ttbarRegionValidation':
        LooseCut   = (tree["Electron_pfRelIso03_all"]  <= 0.4) & lep1_iso & lep2_iso & dr_l1_l2 & ~no_bjets & Z_cut_l1Electron & ~MET_gt50 & Z_cut_l1l2Electron
                
      passLooseWP = LooseCut
      passTightWP = passLooseWP & e_iso 
      passLooseNotTightWP = passLooseWP & ~e_iso

      cut_region[f'{RegionName}_PassLooseWP'] = passLooseWP
      cut_region[f'{RegionName}_PassTightWP'] = passTightWP
      cut_region[f'{RegionName}_PassLooseNotTightWP'] = passLooseNotTightWP

      if type == 'MC':
        ElectronIsPromptLepton =   (tree["Electron_genPartFlav"]  == 1) | (tree["Electron_genPartFlav"]  == 15)

        cut_region[f'{RegionName}_PassLooseWP_ElectronIsPromptLepton'] = passLooseWP & ElectronIsPromptLepton
        cut_region[f'{RegionName}_PassTightWP_ElectronIsPromptLepton'] = passTightWP & ElectronIsPromptLepton
        cut_region[f'{RegionName}_PassLooseNotTightWP_ElectronIsPromptLepton'] = passLooseNotTightWP & ElectronIsPromptLepton

  return cut_region
