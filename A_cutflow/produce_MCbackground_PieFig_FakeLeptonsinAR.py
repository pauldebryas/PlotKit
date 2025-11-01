import os
import yaml
import matplotlib.pyplot as plt

from common.helpers import load_ntuples, dict_files
from common.regions.regions import compute_region_mask

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2017'
#tag = tag used for anatuple production
tag = 'FinalProd'
#region_name = region where MC background is studied
region_name = 'SignalRegion'
#flavor = "Tau", "Muon" or "Electron"
flavor = "Tau"
#----------------------------------------------------------------------------------------------------------------

flavors = {
    "Tau": f"{region_name}_FtauAR_NotTrueLeptons",
    "Muon": f"{region_name}_FmuAR_NotTrueLeptons",
    "Electron": f"{region_name}_FeAR_NotTrueLeptons"
}

if flavor == "Tau":
    channels = ['tee','ttm', 'tte','tmm','tem']
    regionFake = 'FtauAR'
if flavor == "Muon":
    channels = ['tmm','ttm','tem']
    regionFake = 'FmuAR'
if flavor == "Electron":
    channels = ['tee','tte','tem']
    regionFake = 'FeAR'

output_path_fig = os.path.join(os.getenv("RUN_PATH"),f'A_cutflow/figures/{tag}/{period}/')
if not os.path.exists(output_path_fig):
    os.makedirs(output_path_fig) 

def save_sumwEvents(channels, tag, period, region_name, regionFake):
    my_dict = {}
    for channel in channels:
        my_dict[channel] = {}
        print(f'processing channel {channel}')
        path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/{channel}/anatuple/'
        dictfiles = dict_files(path_to_files)
        branches = load_ntuples(dictfiles)
        for FileName in branches.keys():
            if ~FileName.endswith(f'_{period}') & (FileName[:4] != 'HNL_') & (len(branches[FileName]) != 0):
                cut_region = compute_region_mask(branches[FileName], channel, 'MC', region_name)
                my_dict[channel][FileName] = sum(((branches[FileName]['genWeight'])[cut_region[f'{region_name}_{regionFake}_NotTrueLeptons']]).tolist())

    return my_dict

dict = save_sumwEvents(channels, tag, period, region_name, regionFake)

print(f'Plotting only Fake Lepton contribution in {regionFake}')
output_path_fig_file = output_path_fig + f'MCin{region_name}_FakeLepton_{regionFake}'

# Define a dictionary that maps each label to a specific color
label_color_map = {
    'WZ': '#007de9',
    'WminusHToTauTau': '#8567ad',
    'WplusHToTauTau': '#ff7a7a',
    'ZZ': '#f68e2f',
    'DYJetsToLL':'#1daab4',
    'TT': '#238d23',
    'WWW': '#e6ce06',
    'WJetsToLNu': '#4d2e07',
    'WW': '#7fff00',
    'ZHToTauTau': '#b2ffff',
    'TTWW': '#fed0ee',
    'TTZToLLNuNu_M-10': '#2a52be',
    'WWZ': '#e6ce06',
    'ttHToTauTau': '#ffe4e1', 
    'ST': '#fc0fc0', 
    'EWK': '#fc330f', 
    'others': '#717275',
}

MCsamples = dict

for channel in channels:
    print(f'... plotting channel {channel}')
    sumwtotal = sum([MCsamples[channel][MCfile] for MCfile in MCsamples[channel].keys()])
    MCsamples_grouped = {}
    MCsamples_grouped['DYJetsToLL'] = 0
    MCsamples_grouped['EWK'] = 0
    MCsamples_grouped['ST'] = 0
    MCsamples_grouped['TT'] = 0
    MCsamples_grouped['WJetsToLNu'] = 0
    MCsamples_grouped['others'] = 0
    for MCfile in MCsamples[channel].keys():
        if MCfile.startswith('DYJetsToLL_'):
            MCsamples_grouped['DYJetsToLL'] += MCsamples[channel][MCfile]
            continue
        if MCfile.startswith('EWK'):
            MCsamples_grouped['EWK'] += MCsamples[channel][MCfile]
            continue
        if MCfile.startswith('ST_'):
            MCsamples_grouped['ST'] += MCsamples[channel][MCfile]
            continue
        if MCfile.startswith('TTTo'):
            MCsamples_grouped['TT'] += MCsamples[channel][MCfile]
            continue
        if MCfile.endswith('JetsToLNu') | MCfile.startswith('WJetsToLNu_HT'):
            MCsamples_grouped['WJetsToLNu'] += MCsamples[channel][MCfile]
            continue
        if (100*MCsamples[channel][MCfile]/sumwtotal) <= 2.0:
            MCsamples_grouped['others'] += MCsamples[channel][MCfile]
            continue
        MCsamples_grouped[MCfile] = MCsamples[channel][MCfile]

    key_del = []
    for mykey in MCsamples_grouped.keys():
        if (100*MCsamples_grouped[mykey]/sumwtotal) <= 2.0:
            MCsamples_grouped['others'] += MCsamples_grouped[mykey]
            key_del.append(mykey)
            continue

    for kes in key_del:
        del MCsamples_grouped[kes]

    labels = MCsamples_grouped.keys()
    colors = [label_color_map[label] for label in labels]
    sizes = [MCsamples_grouped[MCfile] for MCfile in MCsamples_grouped.keys()] 
    fig, ax = plt.subplots()
    #plt.title(f'MC in Signal Region in {channel}')
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors)
    plt.savefig(output_path_fig_file + f'_{channel}.pdf')
