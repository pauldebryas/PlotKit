import pickle
import os
import yaml
import matplotlib.pyplot as plt
import numpy as np

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
periods = ['2016','2016_HIPM','2017','2018']
#tag = tag used for anatuple production
tag = 'FinalProd'
#channels
channels = ['tte', 'ttm','tee','tem','tmm']
#----------------------------------------------------------------------------------------------------------------

inputs_signal_file = '/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/common/config/all/inputs/inputs_AllSignal.yaml'
with open(inputs_signal_file, 'r') as f:
    inputs_signal = yaml.safe_load(f)

inputs_background_file = '/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/common/config/all/inputs/inputs_AllMCbackground.yaml'
with open(inputs_background_file, 'r') as f:
    inputs_background = yaml.safe_load(f)

figure_path = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/PlotKit/A_cutflow/figures/{tag}/cutflow/'
if not os.path.exists(figure_path):
    os.makedirs(figure_path) 

print('load cutflow')
cutflow_selection_signal = {}
cutflow_selection_background = {}
for period in periods:
    print(f'...{period}')
    cutflow_selection_signal[period] = {}
    cutflow_selection_background[period] = {}
    for channel in channels:
        cutflow_selection_signal[period][channel] = {}
        cutflow_selection_background[period][channel] = {}

        cutflow_path = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'cutflow_pkl')

        # load ntuples/pkl
        cutflow_signal = {}
        for process in inputs_signal:
            if process['type'] == 'signal':
                for file in process['files']:
                    sample = file.replace("_anatuple.root", "")
                    cutflow_file = os.path.join(cutflow_path, file.replace("_anatuple.root", "_cutflow.pkl"))
                    with open(cutflow_file, 'rb') as f:
                        cutflow_signal[sample] = pickle.load(f)
        
        # load ntuples/pkl
        cutflow_background = {}
        for process in inputs_background:
            if process['type'] == 'background':
                for file in process['files']:
                    sample = file.replace("_anatuple.root", "")
                    cutflow_file = os.path.join(cutflow_path, file.replace("_anatuple.root", "_cutflow.pkl"))
                    with open(cutflow_file, 'rb') as f:
                        cutflow_background[sample] = pickle.load(f)

        list_of_keys = list(cutflow_signal['HNL_tau_M-20'].keys())
        cutflow_strings = [string for string in list_of_keys if string.startswith('sumw_')]
        cutflow_strings.remove('sumw_init')
        cutflow_strings.remove('sumw_corrections')

        for step in cutflow_strings:
            cutflow_selection_signal[period][channel][step] = 0
            cutflow_selection_background[period][channel][step] = 0
            for HNLsamples in cutflow_signal.keys():
                cutflow_selection_signal[period][channel][step] += cutflow_signal[HNLsamples][step][HNLsamples]
            for MCbackground in cutflow_background.keys():  
                cutflow_selection_background[period][channel][step] += cutflow_background[MCbackground][step][MCbackground]

title_signal = {
    'tte':r'$\tau \tau e$ cutflow (signal)',
    'ttm':r'$\tau \tau \mu$ cutflow (signal)',
    'tee':r'$\tau e e$ cutflow (signal)',
    'tmm':r'$\tau \mu \mu$ cutflow (signal)',
    'tem':r'$\tau e \mu$ cutflow (signal)',
}

title_background = {
    'tte':r'$\tau \tau e$ cutflow (MC backgrounds)',
    'ttm':r'$\tau \tau \mu$ cutflow (MC backgrounds)',
    'tee':r'$\tau e e$ cutflow (MC backgrounds)',
    'tmm':r'$\tau \mu \mu$ cutflow (MC backgrounds)',
    'tem':r'$\tau e \mu$ cutflow (MC backgrounds)',
}

print('plotting cutflow')
for channel in channels:
    print(f'...{channel}')

    cut_full_name = list(cutflow_selection_signal[period][channel].keys())
    cut = [string.replace('sumw_', '') for string in cut_full_name]
    cut[0] = 'all'

     # signal
    fig, ax = plt.subplots(1,1, figsize=(8, 5))

    for period in periods:
        init = np.array(cutflow_selection_signal[period][channel]['sumw_reweight'])
        cutflow_channel = [cutflow_selection_signal[period][channel][step] for step in cut_full_name] 
        Eff =  cutflow_channel/init
        ax.errorbar(range(len(cut)), Eff, fmt='o--', capsize=5, label=period)
    
    ax.set_ylabel('Efficiency')
    ax.set_yscale('log')
    ax.set_xticks(range(len(cut)))
    ax.set_xticklabels(cut, rotation=45, ha='right')
    ax.grid(axis='x', color='black', linestyle=(0, (1, 10)), linewidth=0.5)
    ax.legend(prop={'size': 7})
    ax.set_title(title_signal[channel])
    plt.savefig(figure_path + f'signal_{channel}.pdf', bbox_inches='tight')

    # background
    fig, ax = plt.subplots(1,1, figsize=(8, 5))

    for period in periods:
        init = np.array(cutflow_selection_background[period][channel]['sumw_reweight'])
        cutflow_channel = [cutflow_selection_background[period][channel][step] for step in cut_full_name] 
        Eff =  cutflow_channel/init
        ax.errorbar(range(len(cut)), Eff, fmt='o--', capsize=5, label=period)
    
    ax.set_ylabel('Efficiency')
    ax.set_yscale('log')
    ax.set_xticks(range(len(cut)))
    ax.set_xticklabels(cut, rotation=45, ha='right')
    ax.grid(axis='x', color='black', linestyle=(0, (1, 10)), linewidth=0.5)
    ax.legend(prop={'size': 7})
    ax.set_title(title_background[channel])
    plt.savefig(figure_path + f'MCbackground_{channel}.pdf', bbox_inches='tight')
