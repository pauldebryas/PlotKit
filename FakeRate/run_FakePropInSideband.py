import os
import yaml
from hist_makers.helpers import load_ntuples
from FakeRate.helpers import compute_FakePropInSideband

#computing P_fake(Sideband) as function of pt/eta of the tau in the MCsample

#parameters
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'AllHNLSamples'
#nbins = number of bins (same for pt and eta) for P_fake(Sideband) estimation.
nbins = 6

output_results = os.path.join(os.getenv("RUN_PATH"),f'FakeRate/results/')
if not os.path.exists(output_results):
    os.makedirs(output_results)

output_figures = os.path.join(os.getenv("RUN_PATH"),f'FakeRate/figures/')
if not os.path.exists(output_figures):
    os.makedirs(output_figures)

inputs = {}
for channel in ['tem','tmm','tee']:
    inputs_cfg_file = f'{os.getenv("RUN_PATH")}/config/shared/inputs_MCbackground_{period}.yaml'
    inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/config/hnl/hnl_{channel}/inputs.yaml'
    input_dir = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')

    with open(inputs_cfg_file, 'r') as f:
        inputs_cfg = yaml.safe_load(f)

    with open(inputs_cfg_file_data, 'r') as f:
        inputs_cfg_data = yaml.safe_load(f)

    input_files = {}
    #add MC
    for input in inputs_cfg:
        if 'files' not in input.keys():
            raise f'Missing files dict for {input}'
        files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
        for file in files_list:
            if os.path.isfile(file) == False:
                print('WARNING: ' + file + ' is missing')
                files_list.remove(file)
        input_files[input['name']] = files_list
    #add data
    for input in inputs_cfg_data:
        if input.get('name') == 'data':
            if 'files' not in input.keys():
                raise f'Missing files dict for data key'
            files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
            for file in files_list:
                if os.path.isfile(file) == False:
                    print('WARNING: ' + file + ' is missing')
                    files_list.remove(file)
            input_files['data'] = files_list

    if input_files['data'] == None:
        raise(f'missing data samples to compute fake prop in Sideband')
    
    branches = load_ntuples(input_files)
    inputs[channel] = branches

# compute fake proportion in sideband 
for channel in ['tmm','tee', 'tem']:
    compute_FakePropInSideband(inputs, channel, n_bins = nbins)