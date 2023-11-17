import numpy as np
import awkward as ak
from hep_ml import reweight
import os
import yaml
import json
from sklearn.model_selection import train_test_split
from helpers import load_input_files, load_ntuples, draw_distributions, compute_ratio_witherr
from hist_makers.regions.GBReweighter_regions import compute_region_mask

#computing weights of BDT for background estimation
#------------------------------------------------------------------------------------
#parameters
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'AllHNLSamples'
#channel = channel of the weights 
channel = 'tmm'
#columns = inputs parameters of the BDT
columns = ['Tau_pt', 'Tau_Jet_pt', 'Tau_eta', 'Tau_charge', 'Tau_decayMode']
#hist_settings = parameters for control plots
hist_settings = {'Tau_pt': {'bins': 100, 
                            'density': True, 
                            'alpha': 0.7},
                'Tau_Jet_pt': {'bins': 100, 
                               'density': True, 
                               'alpha': 0.7},
                'Tau_eta': {'bins': 100, 
                            'density': True, 
                            'alpha': 0.7},
                'Tau_charge': {'bins': 2, 
                               'density': True, 
                               'alpha': 0.7},
                'Tau_decayMode': {'bins': 11, 
                                  'density': True, 
                                  'alpha': 0.7}
                 }
# parameters for BDT
nEstimators=50
LearningRate=0.1
MaxDepth=3
MinSamplesLeaf=1000
GBargs= {'subsample': 0.4}

#------------------------------------------------------------------------------------

inputs_cfg_file = f'{os.getenv("RUN_PATH")}/config/hnl/hnl_{channel}/inputs.yaml'
input_dir = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag, channel, 'anatuple')

with open(inputs_cfg_file, 'r') as f:
    inputs_cfg = yaml.safe_load(f)
input_files = load_input_files(inputs_cfg, input_dir)
if input_files['data'] == None:
    raise(f'missing data samples to compute weights')

savefig_dir = f'{os.getenv("RUN_PATH")}/GBReweighter/figures/{channel}/'
if not os.path.exists(savefig_dir):
    os.makedirs(savefig_dir)

saveresults_dir = f'{os.getenv("RUN_PATH")}/GBReweighter/results/'
if not os.path.exists(saveresults_dir):
    os.makedirs(saveresults_dir)

branches = load_ntuples(input_files['data'])
mask_region_data = compute_region_mask(branches, channel)

#------------------------------------------------------------------------------------
# print number of events in the different regions (check if enought stat)

print('N events in regions:')
print(f"N events in D_f {len(branches[mask_region_data['SidebandZveto_fail']])}")
print(f"N events in C_f {len(branches[mask_region_data['ControlRegion_fail']])}")
print(f"N events in S_f {len(branches[mask_region_data['SignalRegion_fail']])}")
print(f"N events in D_p {len(branches[mask_region_data['SidebandZveto_pass']])}")
print(f"N events in C_p {len(branches[mask_region_data['ControlRegion_pass']])}")
print('')
#------------------------------------------------------------------------------------

original = branches[mask_region_data['SidebandZveto_fail']]
target = branches[mask_region_data['SidebandZveto_pass']]

original = ak.to_dataframe(original[columns])
target = ak.to_dataframe(target[columns])

original_weights = np.ones(len(original))

# divide original samples into training ant test parts
original_train, original_test = train_test_split(original)
# divide target samples into training ant test parts
target_train, target_test = train_test_split(target)

original_weights_train = np.ones(len(original_train))
original_weights_test = np.ones(len(original_test))

#------------------------------------------------------------------------------------
print('Test part for target distribution (without reweighting):')

draw_distributions(original_test, target_test, original_weights_test, 'DistributionsTest', columns, hist_settings, savefig_dir)
print('')
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#print('Bins-based reweighting in n dimensions:')

#bins_reweighter = reweight.BinsReweighter(n_bins=20, n_neighs=1.)
#bins_reweighter.fit(original_train, target_train)

#bins_weights_test = bins_reweighter.predict_weights(original_test)
# validate reweighting rule on the test part comparing 1d projections
#draw_distributions(original_test, target_test, bins_weights_test, 'BinsReweighting', columns, hist_settings, savefig_dir)
#print('')
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
print('Gradient Boosted Reweighter:')

reweighter = reweight.GBReweighter(n_estimators=nEstimators, learning_rate=LearningRate, max_depth=MaxDepth, min_samples_leaf=MinSamplesLeaf, gb_args=GBargs)
reweighter.fit(original_train, target_train)

gb_weights_test = reweighter.predict_weights(original_test)
# validate reweighting rule on the test part comparing 1d projections
draw_distributions(original_test, target_test, gb_weights_test, 'GBReweighting', columns, hist_settings, savefig_dir)
print('')
#------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------
#Gradient Boosted Reweighter over all Z region

reweighter = reweight.GBReweighter(n_estimators=nEstimators, learning_rate=LearningRate, max_depth=MaxDepth, min_samples_leaf=MinSamplesLeaf, gb_args=GBargs)
reweighter.fit(original, target)

C_f = branches[mask_region_data['ControlRegion_fail']]

S_f = branches[mask_region_data['SidebandZveto_fail']]
S_p = branches[mask_region_data['SidebandZveto_pass']]
sf, sf_low, sf_up = compute_ratio_witherr(len(S_p),len(S_f))

gb_weights = reweighter.predict_weights(ak.to_dataframe(C_f[columns]))
event_nb = np.array(C_f['event'])

# saving
dict = {}
dict['event'] = event_nb.tolist()
dict['gb_weights'] = gb_weights.tolist()
dict['sf'] = sf
dict['sf_low'] = sf_low
dict['sf_up'] = sf_up

with open(f'{os.getenv("RUN_PATH")}/GBReweighter/results/gb_weights_ControlRegion_{channel}.json', 'w') as fp:
    json.dump(dict, fp)
#------------------------------------------------------------------------------------