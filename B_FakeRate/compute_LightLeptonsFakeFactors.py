import os
import yaml
import json
import numpy as np

from common.helpers import load_ntuples
from B_FakeRate.helpers import compute_FakeRateDataLL, compute_FakeRateMCLL, PrintContributionMC

#computing P_fake(T|L) as function of pt/eta for light leptons

#parameters -----------------------------------------------------------------------------------------------------
#Light lepton type (Electron or Muon)
namelightlepton = 'Electron'
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'LightLeptFFV2'
#bins = either int and bin size such that statistic is the same for all bins or 'config' and pt and eta bins from config file 
bins = 'config'
#----------------------------------------------------------------------------------------------------------------
RegionNames = ['ttbarRegionFR', 'DYRegionFR']
for RegionName in RegionNames:
    print(f'Computing Fake Factors for {namelightlepton} in {RegionName} region')

    #channel = channel where region is
    if namelightlepton == 'Electron':
        if RegionName in ['ttbarRegionFR']:
            channel = 'lle' 
        if RegionName in ['DYRegionFR']:
            channel = 'Ze' 
    if namelightlepton == 'Muon':
        if RegionName in ['ttbarRegionFR']:
            channel = 'llmu' 
        if RegionName in ['DYRegionFR']:
            channel = 'Zmu'

    # load Corr factor from YAML file
    corrFactor_path = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/correctionFactorsLL.yml')
    with open(corrFactor_path, 'r') as file:
        CorrFactor = yaml.safe_load(file)

    if bins == 'config':
        nbins = []
        with open(f'{os.getenv("RUN_PATH")}/common/config/all/bin_FR.yaml', 'r') as f:
            binsconfig = yaml.safe_load(f)
        nbins.append(binsconfig[namelightlepton]['pt'])
        nbins.append(binsconfig[namelightlepton]['abseta'])
    else:
        nbins = bins

    # load input files
    inputs_cfg_file = f'{os.getenv("RUN_PATH")}/common/config/all/inputs/inputs_AllMCbackground.yaml'
    input_dir = os.path.join('/eos/user/p/pdebryas/HNL_LLFF/anatuple', period, tag, channel, 'anatuple')
    if namelightlepton == 'Electron':
        inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_tte/inputs_data.yaml' 
    elif namelightlepton == 'Muon':
        inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_ttm/inputs_data.yaml' 

    with open(inputs_cfg_file, 'r') as f:
        inputs_cfg = yaml.safe_load(f)

    with open(inputs_cfg_file_data, 'r') as f:
        inputs_cfg_data = yaml.safe_load(f)

    input_files = {}

    exclude_list = [] #['WJetsToLNu_HT-100To200_anatuple.root', 'WJetsToLNu_HT-1200To2500_anatuple.root']
    exclude_list_files = [os.path.join(input_dir, elem) for elem in exclude_list] 

    for input in inputs_cfg:
        if 'files' not in input.keys():
            raise f'Missing files dict for {input}'
        files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
        files_list_new = list(files_list)
        for file in files_list:
            if (os.path.isfile(file) == False) or (file in exclude_list_files):
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

    inputs = load_ntuples(input_files)

    output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/FakeFactors{namelightlepton}/{RegionName}/')
    if not os.path.exists(output_results):
        os.makedirs(output_results)        

    output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactors{namelightlepton}/{RegionName}/')
    if not os.path.exists(output_figures):
        os.makedirs(output_figures)   

    print(f'Computing Fake Factors in Data:')
    #compute_FakeRateData('Tau', inputs, nbins, output_results, output_figures, channels, RegionName)
    compute_FakeRateDataLL(namelightlepton, inputs, nbins, output_results, output_figures, channel, RegionName, CorrFactor)
    print(f'Computing Fake Factors in MC:')
    # use the same binning as in data
    with open(output_results + f'P_fake_TL_{namelightlepton}.json') as f:
        data_tl = json.load(f)
    # Extract the bin edges and contents from the JSON files
    pt_bins = np.array(data_tl['corrections'][0]['data']['edges'][0])  # X-axis bins (pt bins)
    eta_bins = np.array(data_tl['corrections'][0]['data']['edges'][1])  # Y-axis bins (eta bins)
    compute_FakeRateMCLL(namelightlepton, inputs, pt_bins, eta_bins, output_results, output_figures, channel, RegionName, CorrFactor)
    print('')
    inputs_ch = {}
    inputs_ch[channel] = inputs
    PrintContributionMC(namelightlepton, [channel], inputs_ch, RegionName)
