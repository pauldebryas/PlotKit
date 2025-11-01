import os
import yaml

from B_FakeRate.helpers import compute_corrFactor
from common.helpers import load_ntuples

#computing C_mu and C_e in the different channels for light leptons

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'FinalProd'
#RegionName = region where to compute corr factor
RegionName = 'All'
#----------------------------------------------------------------------------------------------------------------

Leptons = ['Muon', 'Electron']
output_factors = {}

for leptonname in Leptons:
    if leptonname == 'Muon':
        channels = ['llmu', 'Zmu']
        inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_ttm/inputs_data.yaml' 
    if leptonname == 'Electron':
        channels = ['lle', 'Ze']
        inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_tte/inputs_data.yaml' 

    input_files = {}
    for channel in channels:
        # load input files
        input_dir = os.path.join('/eos/user/p/pdebryas/HNL_LLFF/anatuple', period, tag, channel, 'anatuple')

        with open(inputs_cfg_file_data, 'r') as f:
            inputs_cfg_data = yaml.safe_load(f)

        files_Data = []
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
                files_Data += files_list_new
        input_files[channel] = files_list_new

    inputs = load_ntuples(input_files)

    output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/')
    if not os.path.exists(output_figures):
        os.makedirs(output_figures)   

    corrFactor = compute_corrFactor(inputs, channels, leptonname, RegionName, 'PassLooseWP', 30, output_figures)
    output_factors[leptonname] = corrFactor

output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/')
if not os.path.exists(output_results):
    os.makedirs(output_results)        

# Save to a YAML file
with open(f'{output_results}correctionFactorsLL.yml', 'w') as file:
    yaml.dump(output_factors, file, default_flow_style=False)

