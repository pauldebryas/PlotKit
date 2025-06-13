import os
import yaml
import numpy as np

from common.helpers import load_ntuples, dict_files
from common.regions.regions import compute_region_mask

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'AddJETcorr'
#region_name = region where proportion of DY/ttbar in MC background is studied (in application region where leptons are not prompt)
region_name = 'InvertedBjetsVetoRegion'# SignalRegion or InvertedBjetsVetoRegion
#----------------------------------------------------------------------------------------------------------------

channels = ['tee','tmm','tem', 'ttm', 'tte']

output_path_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/')
if not os.path.exists(output_path_results):
    os.makedirs(output_path_results) 

def load_sumwEvents(channels, tag, period, region_name):
    my_dict = {}
    my_dict_sq = {}
    for channel in channels:
        my_dict[channel] = {}
        my_dict_sq[channel] = {}
        print(f'processing channel {channel}\n')
        path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/{channel}/anatuple/'
        dictfiles = dict_files(path_to_files)
        branches = load_ntuples(dictfiles)
        for FileName in branches.keys():
            if ~FileName.endswith(f'_{period}') & (FileName[:4] != 'HNL_'):
                if len(branches[FileName]) != 0:
                    cut_region = compute_region_mask(branches[FileName], channel, 'MC', region_name)
                    my_dict[channel][FileName] = sum(((branches[FileName]['genWeight'])[cut_region[f'{region_name}_PassLooseNotTightWP_NotTrueLeptons']]).tolist())
                    my_dict_sq[channel][FileName] = np.sum(np.array((branches[FileName]['genWeight'])[cut_region[f'{region_name}_PassLooseNotTightWP_NotTrueLeptons']]) ** 2)

    return my_dict, my_dict_sq

output_path_results_file = output_path_results + f'Prop_DY_ttbar_in{region_name}AppRegion.yml'

sumwEventsFile, sumw2EventsFile = load_sumwEvents(channels, tag, period, region_name)

print(f'Computing weights: \n')

out_dict = {}
for channel in channels:
    print(f'For channel {channel}')

    out_dict[channel] = {}

    sumwEventsFile_grouped = {}
    sumwEventsFile_grouped['DYJetsToLL'] = 0
    sumwEventsFile_grouped['TT'] = 0

    sumw2EventsFile_grouped = {}
    sumw2EventsFile_grouped['DYJetsToLL'] = 0
    sumw2EventsFile_grouped['TT'] = 0

    for MCfile in sumwEventsFile[channel].keys():
        if MCfile.startswith('DYJetsToLL_'):
            sumwEventsFile_grouped['DYJetsToLL'] += sumwEventsFile[channel][MCfile]
            sumw2EventsFile_grouped['DYJetsToLL'] += sumw2EventsFile[channel][MCfile]
            continue
        if MCfile.startswith('TTTo'):
            sumwEventsFile_grouped['TT'] += sumwEventsFile[channel][MCfile]
            sumw2EventsFile_grouped['TT'] += sumw2EventsFile[channel][MCfile]
            continue
    
    S1 = sumwEventsFile_grouped['DYJetsToLL']
    S2 = sumwEventsFile_grouped['TT']
    S_tot = S1 + S2

    sigma_S1 = np.sqrt(sumw2EventsFile_grouped['DYJetsToLL'])
    sigma_S2 = np.sqrt(sumw2EventsFile_grouped['TT'])

    f1 = S1 / S_tot
    f2 = S2 / S_tot

    sigma_f1 = np.sqrt((S2 / S_tot**2 * sigma_S1) ** 2 + (S1 / S_tot**2 * sigma_S2) ** 2)
    sigma_f2 = sigma_f1  # Since f2 = 1 - f1, its uncertainty is the same

    print(f"DY = {f1} ± {sigma_f1}")
    print(f"ttbar = {f2} ± {sigma_f2}")

    out_dict[channel]['DY'] = np.array([f1, sigma_f1]).tolist()
    out_dict[channel]['ttbar'] = np.array([f2, sigma_f2]).tolist() 

with open(output_path_results_file, 'w') as outfile:
    yaml.dump(out_dict, outfile)

def print_stat_unc_table(data):
    print("{:<10} {:<10} {:<10} {:<10}".format("Category", "Process", "Weight", "Stat. Unc."))
    print("-" * 40)
    
    for category, processes in data.items():
        for process, values in processes.items():
            weight = round(values[0], 3)
            stat_unc = round(values[1], 3)
            print(f"{category:<10} {process:<10} {weight:<10} {stat_unc:<10}")

print(f'Summary: \n')
print_stat_unc_table(out_dict)

