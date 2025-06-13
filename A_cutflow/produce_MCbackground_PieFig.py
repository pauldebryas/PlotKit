import os
from A_cutflow.helpers import save_sumwEvents, plot_pieFig

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2016_HIPM'
#tag = tag used for anatuple production
tag = 'AddJETcorr'
#region_name = region where MC background is studied
region_name = 'SignalRegion'
#----------------------------------------------------------------------------------------------------------------

channels = ['tee','ttm', 'tte','tmm','tem']

output_path_results = os.path.join(os.getenv("RUN_PATH"),f'A_cutflow/results/{tag}/{period}/')
if not os.path.exists(output_path_results):
    os.makedirs(output_path_results) 

output_path_fig = os.path.join(os.getenv("RUN_PATH"),f'A_cutflow/figures/{tag}/{period}/')
if not os.path.exists(output_path_fig):
    os.makedirs(output_path_fig) 

output_path_results_file = output_path_results + f'{region_name}Sumw.yml'

if not os.path.exists(output_path_results_file):
    sumwEvents = save_sumwEvents(channels, tag, period, output_path_results_file, region_name)
else:
    print(f'A MC sumw file already exist: {output_path_results_file}')

print(f'Plotting all MC bacgrounds in {region_name}')
output_path_fig_file = output_path_fig + f'MCin{region_name}'
plot_pieFig(channels, output_path_results_file, output_path_fig_file)

print(f'Plotting only TrueLepton contribution in {region_name}')
output_path_fig_file = output_path_fig + f'MCin{region_name}_TrueLepton'
plot_pieFig(channels, output_path_results_file, output_path_fig_file, onlyTL = True)