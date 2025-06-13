import os
import yaml

from common.helpers import load_ntuples
from B_FakeRate.helpers import compute_FakeRate, print_nFake

#computing P_fake(T|L) as function of pt/eta of the lepton in all MCsample

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2017'
#tag = tag used for anatuple production
tag = 'AddJETcorr'
# sanity check: compute FFs by channel or not
computeFFsByChannel = False
#bins = either int and bin size such that statistic is the same for all bins (same for pt and eta) 
# or 'config' and pt and eta bins from config file 
bins = 'config'
#----------------------------------------------------------------------------------------------------------------
Leptons = ['Muon','Electron','Tau']

for Lepton in Leptons:

    print(f'Computing MC Fake Factors for {Lepton} in Signal Region:')

    if Lepton == 'Tau':
        channels = ['tem','tee','tmm','tte','ttm']
    if Lepton == 'Electron':
        channels = ['tem','tte', 'tee']
    if Lepton == 'Muon':
        channels = ['tem','ttm','tmm']
    if channels == None:
        print('Unidentified Lepton name: can be Tau, Electron or Muon')
        raise

    if bins == 'config':
        nbins = []
        with open(f'{os.getenv("RUN_PATH")}/common/config/all/bin_FR.yaml', 'r') as f:
            binsconfig = yaml.safe_load(f)
        nbins.append(binsconfig[Lepton]['pt'])
        nbins.append(binsconfig[Lepton]['abseta'])
    else:
        nbins = bins

    # load input files
    inputs = {}
    for channel in channels:
        print(f'Load inputs for channel {channel}')
        inputs_cfg_file = f'{os.getenv("RUN_PATH")}/common/config/all/inputs/inputs_AllMCbackground.yaml'
        input_dir = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')

        with open(inputs_cfg_file, 'r') as f:
            inputs_cfg = yaml.safe_load(f)

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
        
        branches = {}
        branches = load_ntuples(input_files)
        inputs[channel] = branches

    output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/FakeFactors{Lepton}/SignalRegion/')
    if not os.path.exists(output_results):
        os.makedirs(output_results)        

    output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactors{Lepton}/SignalRegion/')
    if not os.path.exists(output_figures):
        os.makedirs(output_figures)   
        
    # sanity check: compute FFs by channel
    if computeFFsByChannel:
        for channel in channels:
            #print nb of fake
            print_nFake(inputs, channel, Lepton, 'SignalRegion')
            #compute FakeFactors
            compute_FakeRate(Lepton, inputs, nbins, output_results, output_figures, [channel], 'SignalRegion')

    # Compute FR for all channel combined
    compute_FakeRate(Lepton, inputs, nbins, output_results, output_figures, channels, 'SignalRegion')

    print(f'Computing MC Fake Factors for {Lepton} in SideBand Region (invert Bjets Veto):')

    output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/FakeFactors{Lepton}/InvertedBjetsVetoRegion/')
    if not os.path.exists(output_results):
        os.makedirs(output_results)        

    output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactors{Lepton}/InvertedBjetsVetoRegion/')
    if not os.path.exists(output_figures):
        os.makedirs(output_figures)      

    # sanity check: compute FFs by channel
    if computeFFsByChannel:
        for channel in channels:
            #print nb of fake
            print_nFake(inputs, channel, Lepton, 'InvertedBjetsVetoRegion')
            #compute FakeFactors
            compute_FakeRate(Lepton, inputs, nbins, output_results, output_figures, [channel], 'InvertedBjetsVetoRegion')

    # Compute FR for all channel combined
    compute_FakeRate(Lepton, inputs, nbins, output_results, output_figures, channels, 'InvertedBjetsVetoRegion')
