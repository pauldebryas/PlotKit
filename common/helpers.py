import numpy as np
import uproot
import os
import awkward as ak
import yaml

def get_hnl_masses(period):
    file_path = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/AnatupleProd/config/samples_{period}.yaml'
    
    # Open and parse YAML file
    with open(file_path, 'r') as f:
        data = yaml.safe_load(f)

    # Extract keys that start with 'HNL_tau_M-'
    var_HNL = [key for key in data if key.startswith('HNL_tau_M-')]

    # Collect all 'mass' values
    mass_list = []
    for key in var_HNL:
        mass = data[key].get('mass')
        if isinstance(mass, list):
            mass_list.extend(mass)
        elif mass is not None:
            mass_list.append(mass)
    return mass_list

def dict_files(directory):
    files = {}
    for entry in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, entry)):
            files[entry.replace('_anatuple.root','')] = [os.path.join(directory, entry)]
    return files

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
          #print(f'... loading {namefile}')
          with uproot.open(file_path) as DataUproot:
            if p != 'data':
              try:
                myarray = DataUproot[TreeName].arrays()
                EmptyFile = False
              except KeyError:
                 print(f"Warning: Tree {TreeName} not found in file {file_path}.")
                 myarray = []
                 EmptyFile = True
            else:
              try:
                myarray = DataUproot['Events;1'].arrays()
                EmptyFile = False
              except KeyError:
                 print(f"Warning: Tree {TreeName} not found in file {file_path}.")
                 myarray = []
                 EmptyFile = True
            if not EmptyFile:
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
    x_sorted = np.sort(x[~np.isnan(x)])
    nlen = len(x_sorted)
    return np.array(np.interp(np.linspace(0, nlen, nbin + 1), np.arange(nlen), np.sort(x_sorted)))

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

def have_common_elements(list1, list2):
    return bool(set(list1) & set(list2))
