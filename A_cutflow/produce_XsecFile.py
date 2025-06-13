import os
from A_cutflow.helpers import save_sumwEvents, save_info_XsecUnc

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2016_HIPM'
#tag = tag used for anatuple production
tag = 'AddJETcorr'
#region_name = region where MC background is studied
region_name = 'SignalRegion'
#----------------------------------------------------------------------------------------------------------------

channels = ['ttm', 'tte', 'tee','tmm','tem']

output_path_results = os.path.join(os.getenv("RUN_PATH"),f'A_cutflow/results/{tag}/{period}/')
if not os.path.exists(output_path_results):
    os.makedirs(output_path_results) 

output_path_results_file = output_path_results + f'{region_name}Sumw.yml'
output_path_Xsec_file = output_path_results + f'Xsec_unc.yml'

if not os.path.exists(output_path_results_file):
    sumwEvents = save_sumwEvents(channels, tag, period, output_path_results_file, region_name)
else:
    print(f'A MC sumw file already exist: {output_path_results_file}')

save_info_XsecUnc(channels, output_path_results_file, output_path_Xsec_file)
