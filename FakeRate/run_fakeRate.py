import os
import yaml
from hist_makers.helpers import load_ntuples
from FakeRate.helpers import compute_FakeRate

#computing P_fake(T|L) as function of matching Jet pt/eta of the tau in the MCsample

#parameters
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'AllHNLSamples'
#MCsample = MC sample used for FR computation (can be DY or TT for example)
MCsample = 'DY'
#nbins = number of bins (same for pt and eta) for FR estimation.
nbins = 5

output_results = os.path.join(os.getenv("RUN_PATH"),f'FakeRate/results/')
if not os.path.exists(output_results):
    os.makedirs(output_results)        

output_figures = os.path.join(os.getenv("RUN_PATH"),f'FakeRate/figures/')
if not os.path.exists(output_figures):
    os.makedirs(output_figures)        

inputs = {}
for channel in ['tem','tmm','tee']:
    inputs_cfg_file = f'{os.getenv("RUN_PATH")}/config/shared/inputs_MCbackground_{period}.yaml'
    input_dir = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')

    with open(inputs_cfg_file, 'r') as f:
        inputs_cfg = yaml.safe_load(f)

    input_files = {}
    for input in inputs_cfg:
        if 'files' not in input.keys():
            raise f'Missing files dict for {input}'
        files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
        for file in files_list:
            if os.path.isfile(file) == False:
                print('WARNING: ' + file + ' is missing')
                files_list.remove(file)
        input_files[input['name']] = files_list

    if input_files[MCsample] == None:
        raise(f'missing {MCsample} samples to compute fake rate')
    
    branches = load_ntuples(input_files)
    inputs[channel] = branches

# sanity check: plot FR by channel
for channel in ['tmm','tee', 'tem']:
    compute_FakeRate(inputs, MCsample, channel ,n_bins = nbins)

# Compute FR for all channel combined
compute_FakeRate(inputs, MCsample, 'all' ,n_bins = nbins)
