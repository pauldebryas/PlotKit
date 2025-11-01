import pickle
import os
import yaml
import matplotlib.pyplot as plt
import numpy as np

# parameters -----------------------------------------------------------------------------------------------------
period = '2016_HIPM'  
tag = 'FinalProd'
channels = ['tee','tmm','tte','ttm','tem']
mode = 'full'  # 'full' or 'afterChannelSelection'

# representative HNL mass points and colors
signal_masses = {
    'HNL_tau_M-50': 'brown',
    'HNL_tau_M-100': 'gray',
    'HNL_tau_M-300': 'darkred',
    'HNL_tau_M-500': 'darkslateblue',
    'HNL_tau_M-800': 'darkgreen',
}

# ROOT → matplotlib color map
color_map = {
    "kYellow": "gold",
    "kRed": "red",
    "kPink": "pink",
    "kMagenta": "magenta",
    "kViolet": "violet",
    "kBlue": "blue",
    "kAzure": "deepskyblue",
    "kCyan": "cyan",
    "kTeal": "teal",
    "kBlack": "black",
}

# main backgrounds by channel
main_backgrounds = {
    'tee': ["DY", "TT", "diBoson", "EWK"],
    'tmm': ["DY", "TT", "diBoson", "EWK"],
    'ttm': ["DY", "TT", "diBoson", "EWK", "WJetsToLNu", "ST"],
    'tte': ["DY", "TT", "diBoson", "EWK", "WJetsToLNu", "ST"],
    'tem': ["DY", "TT", "diBoson", "WJetsToLNu", "ST"],
}
#----------------------------------------------------------------------------------------------------------------

inputs_signal_file = '/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/common/config/all/inputs/inputs_AllSignal.yaml'
with open(inputs_signal_file, 'r') as f:
    inputs_signal = yaml.safe_load(f)

inputs_background_file = '/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/common/config/all/inputs/inputs_AllMCbackground.yaml'
with open(inputs_background_file, 'r') as f:
    inputs_background = yaml.safe_load(f)

figure_path = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/A_cutflow/figures/{tag}/{period}/cutflow_signal_vs_bkgTypes/'
if not os.path.exists(figure_path):
    os.makedirs(figure_path) 

print(f'load cutflow for period {period}')
cutflow_selection_signal = {}
cutflow_selection_background = {}

for channel in channels:
    cutflow_selection_signal[channel] = {}
    cutflow_selection_background[channel] = {}

    cutflow_path = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'cutflow_pkl')

    # load signal per sample
    cutflow_signal = {}
    for process in inputs_signal:
        if process['type'] == 'signal':
            for file in process['files']:
                sample = file.replace("_anatuple.root", "")
                cutflow_file = os.path.join(cutflow_path, file.replace("_anatuple.root", "_cutflow.pkl"))
                with open(cutflow_file, 'rb') as f:
                    cutflow_signal[sample] = pickle.load(f)

    # load background per group
    cutflow_background = {}
    for process in inputs_background:
        if process['type'] == 'background':
            group_name = process['name']
            cutflow_background[group_name] = {}
            for file in process['files']:
                sample = file.replace("_anatuple.root", "")
                cutflow_file = os.path.join(cutflow_path, file.replace("_anatuple.root", "_cutflow.pkl"))
                with open(cutflow_file, 'rb') as f:
                    cutflow_background[group_name][sample] = pickle.load(f)

    # steps (take from a known signal sample)
    list_of_keys = list(cutflow_signal['HNL_tau_M-20'].keys())
    cutflow_strings = [string for string in list_of_keys if string.startswith('sumw_')]
    cutflow_strings.remove('sumw_init')
    cutflow_strings.remove('sumw_corrections')

    # initialize storage
    cutflow_selection_signal[channel] = {mass: {step:0 for step in cutflow_strings} for mass in signal_masses.keys()}
    cutflow_selection_background[channel] = {}

    # fill signal
    for mass in signal_masses.keys():
        if mass in cutflow_signal:
            for step in cutflow_strings:
                cutflow_selection_signal[channel][mass][step] = cutflow_signal[mass][step][mass]

    # fill backgrounds grouped
    for group in inputs_background:
        gname = group['name']
        cutflow_selection_background[channel][gname] = {step:0 for step in cutflow_strings}
        for step in cutflow_strings:
            for sample in group['files']:
                sample_name = sample.replace("_anatuple.root", "")
                if sample_name in cutflow_background[gname]:
                    cutflow_selection_background[channel][gname][step] += cutflow_background[gname][sample_name][step][sample_name]

# titles -----------------------------------------------------------------------------------------------------
title_combined = {
    'tte':r'tte baseline cutflow',
    'ttm':r'ttm baseline cutflow',
    'tee':r'tee baseline cutflow',
    'tmm':r'tmm baseline cutflow',
    'tem':r'tem baseline cutflow',
}

# plotting -----------------------------------------------------------------------------------------------------
print('plotting combined cutflow with main backgrounds + other category')
for channel in channels:
    print(f'...{channel}')

    cut_full_name = list(next(iter(cutflow_selection_signal[channel].values())).keys())
    cut = [string.replace('sumw_', '') for string in cut_full_name]
    cut[0] = 'all'
    initial_steps = 'sumw_reweight'

    #only cutflow steps after channel_selection step
    if mode == 'afterChannelSelection':
        cut = ['channel_selection', 'tight_leptons_veto', 'HLT', 'Trigger_lepton_selection', 'Second_lepton_selection', 'Third_lepton_selection']
        cut_full_name = ['sumw_channel_selection', 'sumw_tight_leptons_veto', 'sumw_HLT', 'sumw_Trigger_lepton_selection', 'sumw_Second_lepton_selection', 'sumw_Third_lepton_selection']
        initial_steps = 'sumw_channel_selection'

    fig, ax = plt.subplots(1,1, figsize=(8, 5))

    signal_handles = []
    bkg_handles = []

    # signals
    for mass, color in signal_masses.items():
        if mass in cutflow_selection_signal[channel]:
            init_sig = np.array(cutflow_selection_signal[channel][mass][initial_steps])
            cutflow_channel_sig = [cutflow_selection_signal[channel][mass][step] for step in cut_full_name] 
            Eff_sig = cutflow_channel_sig/init_sig
            h = ax.errorbar(range(len(cut)), Eff_sig, fmt='o--', color=color, capsize=5, label=mass)
            signal_handles.append(h)

    # backgrounds (main ones + merge others)
    main_bkgs = main_backgrounds[channel]
    other_bkg = {step:0 for step in cut_full_name}

    for group in inputs_background:
        gname = group['name']
        color_name = group.get('color', 'black')
        color = color_map.get(color_name, 'black')

        if gname in main_bkgs:
            init_bkg = np.array(cutflow_selection_background[channel][gname][initial_steps])
            cutflow_channel_bkg = [cutflow_selection_background[channel][gname][step] for step in cut_full_name] 
            Eff_bkg = cutflow_channel_bkg/init_bkg
            h = ax.errorbar(range(len(cut)), Eff_bkg, fmt='s-', color=color, capsize=5, label=gname)
            bkg_handles.append(h)
        else:
            for step in cut_full_name:
                other_bkg[step] += cutflow_selection_background[channel][gname][step]

    # add merged "Other MC background"
    if sum(other_bkg.values()) > 0:
        init_bkg = np.array(other_bkg[initial_steps])
        cutflow_channel_bkg = [other_bkg[step] for step in cut_full_name] 
        Eff_bkg = cutflow_channel_bkg/init_bkg
        h = ax.errorbar(range(len(cut)), Eff_bkg, fmt='s-', color='lightgray', capsize=5, label='Other backgrounds')
        bkg_handles.append(h)

    ax.set_ylabel('Cut efficiency')
    if mode!= 'afterChannelSelection':
        ax.set_yscale('log')
    ax.set_xticks(range(len(cut)))
    ax.set_xticklabels(cut, rotation=45, ha='right')
    ax.grid(axis='x', color='black', linestyle=(0, (1, 10)), linewidth=0.5)

    # split legend into signal and background
    first_legend = ax.legend(handles=signal_handles, title="Signal", loc='upper right', fontsize=7)
    ax.add_artist(first_legend)
    ax.legend(handles=bkg_handles, title="MC Backgrounds", loc='lower left', fontsize=7, ncol=2)

    #ax.set_title(title_combined[channel] + f" ({period})")
    if mode == 'afterChannelSelection':
        plt.savefig(figure_path + f'BaselineCutflow_{channel}_{period}_afterChannelSelection.pdf', bbox_inches='tight')
    else:
        plt.savefig(figure_path + f'BaselineCutflow_{channel}_{period}.pdf', bbox_inches='tight')
