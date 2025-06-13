import numpy as np
import ROOT
import uproot
import os
import awkward as ak
import yaml

def compute_var_to_plot(tree, hist_name, CorrFactor = None):

  if hist_name == 'mT_tee_tauhOSL':
    var_to_plot = mT_l1l2(tree, 'Tau', 'OSLepton')
    return var_to_plot

  if hist_name == 'mT_tmm_tauhOSL':
    var_to_plot = mT_l1l2(tree, 'Tau', 'OSLepton')
    return var_to_plot

  if hist_name == 'mT_tem_tauhOSL':
    var_to_plot = mT_l1l2(tree, 'Tau', 'OSLepton')
    return var_to_plot
   
  if hist_name == 'mT_tautau':
    var_to_plot = mT_l1l2(tree, 'Tau1', 'Tau2')
    return var_to_plot

  if hist_name == 'mT_etau':
    var_to_plot = mT_l1l2(tree, 'Electron', 'Tau')
    return var_to_plot

  if hist_name == 'mT_ee':
    var_to_plot = mT_l1l2(tree, 'Electron1', 'Electron2')
    return var_to_plot
  
  if hist_name == 'mT_mumu':
    var_to_plot = mT_l1l2(tree, 'Muon1', 'Muon2')
    return var_to_plot

  if hist_name == 'dr_etau':
    var_to_plot = delta_r(tree, 'Electron', 'Tau1')
    return var_to_plot

  if hist_name == 'dr_mutau':
    var_to_plot = delta_r(tree, 'Muon', 'Tau1')
    return var_to_plot

  if hist_name == 'dr_emu':
    var_to_plot = delta_r(tree, 'Electron', 'Muon')
    return var_to_plot

  if hist_name == 'dr_ee':
    var_to_plot = delta_r(tree, 'Electron1', 'Electron2')
    return var_to_plot

  if hist_name == 'dr_mumu':
    var_to_plot = delta_r(tree, 'Muon1', 'Muon2')
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

  if hist_name[-7:] == '_ptcorr':
    var_to_plot = lepton_pt_corr(tree, hist_name[:-7], CorrFactor)
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

  if hist_name[-7:] == '_abseta':
    var_to_plot = np.abs(lepton_eta(tree, hist_name[:-7]))
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
  

  if hist_name == 'mT_total_ttm':
    var_to_plot = mT_total_l1l2l3(tree, 'Muon', 'Tau1','Tau2')
    return var_to_plot

  if hist_name == 'mT_total_tte':
    var_to_plot = mT_total_l1l2l3(tree, 'Electron', 'Tau1','Tau2')
    return var_to_plot

  if hist_name == 'mT_total_tem':
    var_to_plot = mT_total_l1l2l3(tree, 'Tau', 'Electron','Muon')
    return var_to_plot
  
  if hist_name == 'mT_total_tee':
    var_to_plot = mT_total_l1l2l3(tree, 'Tau', 'Electron1','Electron2')
    return var_to_plot

  if hist_name == 'mT_total_tmm':
    var_to_plot = mT_total_l1l2l3(tree, 'Tau', 'Muon1','Muon2')
    return var_to_plot
  
  
  if hist_name == 'mT_total_tautau':
    var_to_plot = mT_total_l1l2(tree, 'Tau1', 'Tau2')
    return var_to_plot
  

  if hist_name[-8:] == 'eta_sumw':
    var_to_plot = eta_sumw(tree, hist_name[0:3])
    return var_to_plot
  
  if hist_name.startswith('DNNscore'):
    return tree['DNNscore']
  
  print('Error compute_var_to_plot function do not find the variable')
  return 

def eta_sumw(tree, channel):
  if channel == 'tte':
    return (tree[f'Tau1_eta'] +tree[f'Tau2_eta'] +tree[f'Electron_eta'])/3
  if channel == 'ttm':
    return (tree[f'Tau1_eta'] +tree[f'Tau2_eta'] +tree[f'Muon_eta'])/3
  if channel == 'tee':
    return (tree[f'Tau_eta'] +tree[f'Electron1_eta'] +tree[f'Electron2_eta'])/3
  if channel == 'tmm':
    return (tree[f'Tau_eta'] +tree[f'Muon1_eta'] +tree[f'Muon2_eta'])/3
  if channel == 'tem':
    return (tree[f'Tau_eta'] +tree[f'Electron_eta'] +tree[f'Muon_eta'])/3  
  
def lepton_charge(tree, Lepton):
    return tree[f'{Lepton}_charge'] 

def lepton_pt(tree, Lepton):
    return tree[f'{Lepton}_pt']

def lepton_pt_corr(tree, Lepton, CorrFactor):
    if CorrFactor == None:
       print('Error: Missing CorrFactor')
    lepton_pt = tree[f'{Lepton}_pt']
    parton_pt = np.concatenate(tree[f'{Lepton}_ConeCorrectedPt'])* CorrFactor[Lepton]
    pt_corr = np.where(tree[f'{Lepton}_pfRelIso03_all'] < 0.15, lepton_pt, parton_pt)
    return pt_corr

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

def ToRootHist_val(h, flow = False):
    x_bins = ListToVector(h.axes[0].edges, type="double")

    # ROOT histogram without explicit overflow bin
    root_hist = ROOT.TH1D('', '', len(x_bins) - 1, x_bins.data())

    # Fill normal bins
    for i in range(len(x_bins) - 1):
        root_hist.SetBinContent(i + 1, h.values()[i])
        root_hist.SetBinError(i + 1, h.variances()[i] ** 0.5)

    if flow == True:
      last_bin_idx = len(x_bins) - 1  # ROOT bins are 1-based
      # Add overflow content to the last bin
      overflow_content = h.values(flow=True)[-1]  # Overflow bin content
      overflow_variance = h.variances(flow=True)[-1]  # Overflow variance
      root_hist.SetBinContent(last_bin_idx+1, overflow_content)
      root_hist.SetBinError(last_bin_idx+1, overflow_variance ** 0.5)
      # Add underflow content to the last bin
      underflow_content = h.values(flow=True)[0]  # Underflow bin content
      underflow_variance = h.variances(flow=True)[0]  # Underflow variance
      root_hist.SetBinContent(0, underflow_content)
      root_hist.SetBinError(0, underflow_variance ** 0.5)

    root_hist.SetDirectory(0)  # Prevent ROOT from managing memory
    return root_hist

def ToRootHist(h, hist_name, flow=False):
    
    if hist_name == 'FakeBackground':
       sys_fraction = 0.2
    else:
       sys_fraction = 0.

    x_bins = ListToVector(h.axes[0].edges, type="double")

    root_hist = ROOT.TH1D('', '', len(x_bins) - 1, x_bins.data())

    for i in range(len(x_bins) - 1):
        stat_unc = h.variances()[i] ** 0.5  # Statistical error
        sys_unc = sys_fraction * h.values()[i]  # 20% systematic error
        total_unc = (stat_unc**2 + sys_unc**2) ** 0.5  # Combine in quadrature

        root_hist.SetBinContent(i + 1, h.values()[i])
        root_hist.SetBinError(i + 1, total_unc)

    if flow:
        last_bin_idx = len(x_bins) - 1
        
        # Overflow bin
        overflow_content = h.values(flow=True)[-1]
        overflow_stat_unc = h.variances(flow=True)[-1] ** 0.5
        overflow_sys_unc = sys_fraction * overflow_content
        root_hist.SetBinContent(last_bin_idx + 1, overflow_content)
        root_hist.SetBinError(last_bin_idx + 1, (overflow_stat_unc**2 + overflow_sys_unc**2) ** 0.5)

        # Underflow bin
        underflow_content = h.values(flow=True)[0]
        underflow_stat_unc = h.variances(flow=True)[0] ** 0.5
        underflow_sys_unc = sys_fraction * underflow_content
        root_hist.SetBinContent(0, underflow_content)
        root_hist.SetBinError(0, (underflow_stat_unc**2 + underflow_sys_unc**2) ** 0.5)

    root_hist.SetDirectory(0)
    return root_hist

def load_ntuples(samples, TreeName = 'Events;1'):
  processes = samples.keys()
  branches = {}
  for p in processes:
      list = []
      namefile_list = []
      for file_path in samples[p]:
          file_path_split = file_path.split("/")
          namefile_list.append(file_path_split[-1].replace('_anatuple.root', ''))

      for file_path in samples[p]:
          file_path_split = file_path.split("/")
          namefile = file_path_split[-1].replace('_anatuple.root', '')
          # Load the dataframes
          #print(namefile)
          with uproot.open(file_path) as DataUproot:
            if p != 'data':
              myarray = DataUproot[TreeName].arrays()
            else:
              myarray = DataUproot['Events;1'].arrays()
            for file_name in namefile_list:
              if file_name == namefile:
                myarray[f'mask_{file_name}'] = ak.ones_like(myarray['genWeight'])
              else:
                myarray[f'mask_{file_name}'] = ak.zeros_like(myarray['genWeight'])
            list.append(myarray)
            #test_list = np.concatenate(list)
      if len(list) != 0:
        branches[p] = np.concatenate(list)
      else:
        branches[p] = 0
  return branches

def load_Tree(input_files):
  file_path = input_files['TrueLepton'][0]
  DataUproot = uproot.open(file_path)
  return DataUproot.keys()

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
  if x_bins == 'config':
      namelightlepton = hist_name.split('_')[0]
      var = hist_name.split('_')[1]
      with open(f'{os.getenv("RUN_PATH")}/common/config/all/bin_FR.yaml', 'r') as f:
          binsconfig = yaml.safe_load(f)
      pt_bin = binsconfig[namelightlepton]['pt']
      pt_bin[-1] = 100
      eta_bin = binsconfig[namelightlepton]['abseta']
      if var in ['pt', 'ptcorr']:
        x_bins = pt_bin
      if var in ['abseta']:
         x_bins = eta_bin
  elif not isinstance(x_bins, list):
    n_bins, bin_range = x_bins.split('|')
    start,stop = bin_range.split(':')
    if n_bins == 'auto':
      x_bins = [int(start),int(stop)]
    else:
      x_bins = np.arange(int(start),int(stop), (int(stop)-int(start))/int(n_bins))
  return x_bins

def mT_total_l1l2l3(tree, Lepton1, Lepton2, Lepton3):
  '''
  Calculates total Tranverse mass between three particles v1, v2, v3
  ''' 
  MT_total = np.sqrt(MT(tree,Lepton1,Lepton3)**2 + MT(tree,Lepton1,Lepton2)**2 + MT(tree,Lepton2,Lepton3)**2 + MT_MET(tree,Lepton1)**2 + MT_MET(tree,Lepton2)**2 + MT_MET(tree,Lepton3)**2)
  return MT_total

def MT(tree, Lepton1, Lepton2):
  v1_pt = tree[f'{Lepton1}_pt'] 
  v2_pt = tree[f'{Lepton2}_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  v2_phi = tree[f'{Lepton2}_phi']
  v1_mass = tree[f'{Lepton1}_mass'] 
  v2_mass = tree[f'{Lepton2}_mass'] 
  v1_energy = np.sqrt(v1_mass**2 + v1_pt**2)
  v2_energy = np.sqrt(v2_mass**2 + v2_pt**2)
  cos_v1_v2 = np.cos(abs(v1_phi - v2_phi))

  MT = np.sqrt(v1_mass**2 + v2_mass**2 + 2*(v1_energy*v2_energy-v1_pt*v2_pt*cos_v1_v2))
  return MT

def MT_MET(tree, Lepton1):
  v1_pt = tree[f'{Lepton1}_pt'] 
  MET_pt = tree[f'MET_pt'] 
  v1_phi = tree[f'{Lepton1}_phi'] 
  MET_phi = tree[f'MET_phi']
  v1_mass = tree[f'{Lepton1}_mass'] 
  v1_energy = np.sqrt(v1_mass**2 + v1_pt**2)
  cos_v1_MET = np.cos(abs(v1_phi - MET_phi))

  MT = np.sqrt(v1_mass**2 + 2*(v1_energy*MET_pt-v1_pt*MET_pt*cos_v1_MET))
  return MT
   
def make_root_file(path_root_file, nameDir, hists):
  myfile = ROOT.TFile(path_root_file, 'RECREATE' )
  ROOT.gDirectory.mkdir(nameDir)
  ROOT.gDirectory.cd(nameDir)
  for elem in hists:    
      if elem == 'data':
          hists[elem].Write(f'data_obs')
      else:
          hists[elem].Write(f'{elem}')
  myfile.Close()   
  return

import ROOT

def load_root_file(path_root_file, nameDir):
        
  myFile = ROOT.TFile.Open(path_root_file)

  if not myFile or myFile.IsZombie():
      print("Error: Unable to open file", path_root_file)
      return None

  hists = {}
  file_dir = myFile.Get(nameDir)

  if not file_dir:
      print("Error: Unable to find directory in file", nameDir)
      myFile.Close()
      return None

  # Iterate over all keys in the directory
  for key in file_dir.GetListOfKeys():
      obj = key.ReadObj()
      # Check if the object is a histogram
      if isinstance(obj, ROOT.TH1):
          hist_clone = obj.Clone()  # Clone the histogram
          hist_clone.SetDirectory(0)  # Detach from the TFile memory
          ROOT.SetOwnership(hist_clone, True)  # Give ownership to Python

          if key.GetName() == 'data_obs':
              hists['data'] = hist_clone
          else:
              if '_' not in key.GetName():
                  hists[key.GetName()] = hist_clone
  
  myFile.Close()
  
  return hists

'''
def load_root_file(path_root_file, nameDir):
        
  myFile = ROOT.TFile.Open(path_root_file)

  if not myFile or myFile.IsZombie():
      print("Error: Unable to open file", myFile)
      return

  hists = {}
  file_dir = myFile.Get(nameDir)

  if not file_dir:
      print("Error: Unable to find directory in file", myFile)
      myFile.Close()
      return 

  # Iterate over all keys in the directory
  for key in file_dir.GetListOfKeys():
      obj = key.ReadObj()
      # Check if the object is a histogram
      if isinstance(obj, ROOT.TH1):
          # Store the histogram in the dictionary
          if key.GetName() == 'data_obs':
              hists['data'] = obj.Clone()
          else:
              if '_' not in key.GetName():
                  hists[key.GetName()] = obj.Clone()
  
  myFile.Close()
  
  return hists
'''