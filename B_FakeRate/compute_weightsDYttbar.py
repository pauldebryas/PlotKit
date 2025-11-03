import os
import yaml
import numpy as np

from common.helpers import load_ntuples, dict_files
from common.regions.regions import compute_region_mask

# parameters -----------------------------------------------------------------------------------------------------
period = '2017'
tag = 'FinalProd'
region_name = 'SignalRegion'  # SignalRegion or InvertedBjetsVetoRegion or ControlRegion
# ----------------------------------------------------------------------------------------------------------------

#channels = ['tee', 'tmm', 'tem', 'ttm', 'tte', 'tee_ss', 'tee_os', 'tmm_ss', 'tmm_os']
channels = ['tee_ss','tee_os', 'tmm_ss','tmm_os']

# Define flavors, masks, and valid channels
# flavors = {
#     "Tau": {
#         "mask": f"{region_name}_FtauAR_NotTrueLeptons",
#         "channels": ['tee', 'tmm', 'tem', 'ttm', 'tte', 'tee_ss', 'tee_os', 'tmm_ss', 'tmm_os']
#     },
#     "Muon": {
#         "mask": f"{region_name}_FmuAR_NotTrueLeptons",
#         "channels": ['tmm', 'tem', 'ttm', 'tmm_ss', 'tmm_os']
#     },
#     "Electron": {
#         "mask": f"{region_name}_FeAR_NotTrueLeptons",
#         "channels": ['tee', 'tem', 'tte', 'tee_ss', 'tee_os']
#     }
# }

flavors = {
    "Tau": {
        "mask": f"{region_name}_FtauAR_NotTrueLeptons",
        "channels": ['tee_ss','tee_os', 'tmm_ss','tmm_os']
    },
    "Muon": {
        "mask": f"{region_name}_FmuAR_NotTrueLeptons",
        "channels": ['tmm_ss','tmm_os']
    },
    "Electron": {
        "mask": f"{region_name}_FeAR_NotTrueLeptons",
        "channels": ['tee_ss','tee_os']
    }
}

output_path_results = os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/')
if not os.path.exists(output_path_results):
    os.makedirs(output_path_results)

def load_sumwEvents(flavors, tag, period, region_name):
    my_dict = {flavor: {ch: {} for ch in data["channels"]} for flavor, data in flavors.items()}
    my_dict_sq = {flavor: {ch: {} for ch in data["channels"]} for flavor, data in flavors.items()}

    for flavor, data in flavors.items():
        mask_name = data["mask"]
        for channel in data["channels"]:
            print(f'processing flavor {flavor}, channel {channel}\n')
            path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/{channel}/anatuple/'
            dictfiles = dict_files(path_to_files)
            branches = load_ntuples(dictfiles)

            for FileName in branches.keys():
                if ~FileName.endswith(f'_{period}') & (FileName[:4] != 'HNL_'):
                    if len(branches[FileName]) != 0:
                        cut_region = compute_region_mask(branches[FileName], channel, 'MC', region_name)
                        mask = cut_region[mask_name]
                        weights = branches[FileName]['genWeight'][mask]
                        if len(weights) > 0:
                            my_dict[flavor][channel][FileName] = np.sum(weights)
                            my_dict_sq[flavor][channel][FileName] = np.sum(weights**2)

    return my_dict, my_dict_sq


output_path_results_file = output_path_results + f'FRweights_{region_name}AppRegion.yml'

sumwEventsFile, sumw2EventsFile = load_sumwEvents(flavors, tag, period, region_name)

print(f'Computing weights: \n')

out_dict = {}
for flavor, data in flavors.items():
    out_dict[flavor] = {}
    for channel in data["channels"]:
        #print(f'For flavor {flavor}, channel {channel}')

        out_dict[flavor][channel] = {}

        sumwEventsFile_grouped = {'DY': 0, 'TT': 0}
        sumw2EventsFile_grouped = {'DY': 0, 'TT': 0}

        for MCfile in sumwEventsFile[flavor][channel].keys():
            if MCfile.startswith('DYJetsToLL_') or MCfile.startswith('W1JetsToLNu_') or MCfile.startswith('W2JetsToLNu_') or MCfile.startswith('W3JetsToLNu_') or MCfile.startswith('W4JetsToLNu_') or MCfile.startswith('WJetsToLNu_'):
                sumwEventsFile_grouped['DY'] += sumwEventsFile[flavor][channel][MCfile]
                sumw2EventsFile_grouped['DY'] += sumw2EventsFile[flavor][channel][MCfile]
                continue
            if MCfile.startswith('TTTo') or MCfile.startswith('ST_'):
                sumwEventsFile_grouped['TT'] += sumwEventsFile[flavor][channel][MCfile]
                sumw2EventsFile_grouped['TT'] += sumw2EventsFile[flavor][channel][MCfile]
                continue

        S1 = sumwEventsFile_grouped['DY']
        S2 = sumwEventsFile_grouped['TT']
        S_tot = S1 + S2

        if S_tot == 0:
            print(f"Warning: no events for {flavor}, {channel}")
            continue

        sigma_S1 = np.sqrt(sumw2EventsFile_grouped['DY'])
        sigma_S2 = np.sqrt(sumw2EventsFile_grouped['TT'])

        f1 = S1 / S_tot
        f2 = S2 / S_tot

        sigma_f1 = np.sqrt((S2 / S_tot**2 * sigma_S1) ** 2 + (S1 / S_tot**2 * sigma_S2) ** 2)
        sigma_f2 = sigma_f1  # Since f2 = 1 - f1

        #print(f"DY = {f1} ± {sigma_f1}")
        #print(f"ttbar = {f2} ± {sigma_f2}")

        out_dict[flavor][channel]['DY'] = np.array([f1, sigma_f1]).tolist()
        out_dict[flavor][channel]['ttbar'] = np.array([f2, sigma_f2]).tolist()

with open(output_path_results_file, 'w') as outfile:
    yaml.dump(out_dict, outfile)


def print_stat_unc_table(data):
    print("{:<10} {:<10} {:<10} {:<10} {:<10}".format("Flavor", "Channel", "Process", "Weight", "Stat. Unc."))
    print("-" * 60)

    for flavor, categories in data.items():
        for channel, processes in categories.items():
            for process, values in processes.items():
                weight = round(values[0], 3)
                stat_unc = round(values[1], 3)
                print(f"{flavor:<10} {channel:<10} {process:<10} {weight:<10} {stat_unc:<10}")


print(f'Summary: \n')
print_stat_unc_table(out_dict)
