import os
import yaml
import json
import numpy as np

from common.helpers import load_ntuples
from B_FakeRate.helpers import compute_FakeRateData, compute_FakeRateMC, PrintContributionMC

#computing P_fake(T|L) as function of pt/eta for hadronically decayed tau leptons

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'FinalProd'
#bins = either int and bin size such that statistic is the same for all bins (same for pt and eta) or 'config' and pt and eta bins from config file 
bins = 'config'
#----------------------------------------------------------------------------------------------------------------

RegionNames = ['DYRegionFR', 'ttbarRegionFR']
for RegionName in RegionNames:
    print(f'Computing Fake Factors for Tau in {RegionName} Region:')

    if RegionName == 'DYRegionFR':
        channels = ['tll']
    if RegionName == 'ttbarRegionFR':
        channels = ['tem']
    
    if bins == 'config':
        nbins = []
        with open(f'{os.getenv("RUN_PATH")}/common/config/all/bin_FR.yaml', 'r') as f:
            binsconfig = yaml.safe_load(f)
        nbins.append(binsconfig['Tau']['pt'])
        nbins.append(binsconfig['Tau']['abseta'])
    else:
        nbins = bins

    # load input files
    inputs = {}
    for channel in channels:
        print(f'Load inputs for channel {channel}')
        inputs_cfg_file = f'{os.getenv("RUN_PATH")}/common/config/all/inputs/inputs_AllMCbackground.yaml'
        inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_{channel}/inputs_data.yaml' 
        input_dir = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')

        with open(inputs_cfg_file, 'r') as f:
            inputs_cfg = yaml.safe_load(f)

        with open(inputs_cfg_file_data, 'r') as f:
            inputs_cfg_data = yaml.safe_load(f)

        input_files = {}
        for input in inputs_cfg:
            if 'files' not in input.keys():
                raise f'Missing files dict for {input}'
            files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
            files_list_new = list(files_list)
            for file in files_list:
                if os.path.isfile(file) == False:
                    print('WARNING: ' + file + ' is missing')
                    files_list_new.remove(file)
            input_files[input['name']] = files_list_new

        for input in inputs_cfg_data:
            if input['name'] == 'data':
                if 'files' not in input.keys():
                    raise f'Missing files dict for {input}'
                files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
                files_list_new = list(files_list)
                for file in files_list:
                    if os.path.isfile(file) == False:
                        print('WARNING: ' + file + ' is missing')
                        files_list_new.remove(file)
                input_files[input['name']] = files_list_new

        branches = {}
        branches = load_ntuples(input_files)
        inputs[channel] = branches

    output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/FakeFactorsTau/{RegionName}/')
    if not os.path.exists(output_results):
        os.makedirs(output_results)        

    output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactorsTau/{RegionName}/')
    if not os.path.exists(output_figures):
        os.makedirs(output_figures)   
        
    print(f'Computing Fake Factors in Data:')
    compute_FakeRateData('Tau', inputs, nbins, output_results, output_figures, channels, RegionName)
    print(f'Computing Fake Factors in MC:')
    # use the same binning as in data
    # Load the JSON file directly
    with open(output_results + 'P_fake_TL_Tau.json') as f:
        data_tl_tau = json.load(f)
    # Extract the bin edges and contents from the JSON files
    pt_bins = np.array(data_tl_tau['corrections'][0]['data']['edges'][0])  # X-axis bins (pt bins)
    eta_bins = np.array(data_tl_tau['corrections'][0]['data']['edges'][1])  # Y-axis bins (eta bins)
    compute_FakeRateMC('Tau', inputs, pt_bins, eta_bins, output_results, output_figures, channels, RegionName)
    print('')
    PrintContributionMC('Tau',channels, inputs, RegionName)
