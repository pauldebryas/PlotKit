import os
import pandas as pd
from torch import load
import pickle
import json

from C_DNNvar.helpers import process_channels, get_dnn_score_dict_torch_OddEven, DNN_flexible
from C_DNNvar.data_extractor_class import flatten_2D_list
from common.helpers import get_hnl_masses

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for anatuple production
tag = 'AddJETcorr'
#folderDNNName = name of the folder in which DNN model is stored
folderDNNName= 'AllYears_AllChannels'
#modelName = name of the DNN model you want to use
modelName = 'DNN3'
#----------------------------------------------------------------------------------------------------------------

features =  ['event', 'genWeight', 
            'charge_1', 'charge_2', 'charge_3', 
            'pt_1', 'pt_2', 'pt_3', 'pt_MET', 
            'eta_1', 'eta_2', 'eta_3',
            'mass_1', 'mass_2', 'mass_3',
            'phi_1', 'phi_2', 'phi_3', 'phi_MET', 
            'deltaphi_12', 'deltaphi_13', 'deltaphi_23', 
            'deltaphi_1MET', 'deltaphi_2MET', 'deltaphi_3MET',
            ['deltaphi_1(23)', 'deltaphi_2(13)', 'deltaphi_3(12)', 
            'deltaphi_MET(12)', 'deltaphi_MET(13)', 'deltaphi_MET(23)',
            'deltaphi_1(2MET)', 'deltaphi_1(3MET)', 'deltaphi_2(1MET)', 'deltaphi_2(3MET)', 'deltaphi_3(1MET)', 'deltaphi_3(2MET)'],
            'deltaeta_12', 'deltaeta_13', 'deltaeta_23', 
            ['deltaeta_1(23)', 'deltaeta_2(13)', 'deltaeta_3(12)'],
            'deltaR_12', 'deltaR_13', 'deltaR_23', 
            ['deltaR_1(23)', 'deltaR_2(13)', 'deltaR_3(12)'],
            'pt_123',
            'mt_12', 'mt_13', 'mt_23', 
            'mt_1MET', 'mt_2MET', 'mt_3MET',
            ['mt_1(23)', 'mt_2(13)', 'mt_3(12)',
            'mt_MET(12)', 'mt_MET(13)', 'mt_MET(23)',
            'mt_1(2MET)', 'mt_1(3MET)', 'mt_2(1MET)', 'mt_2(3MET)', 'mt_3(1MET)', 'mt_3(2MET)'],
            'mass_12', 'mass_13', 'mass_23',
            'mass_123',
            'Mt_tot',
            ['HNL_CM_angle_with_MET_1', 'HNL_CM_angle_with_MET_2'], 
            ['W_CM_angle_to_plane_1', 'W_CM_angle_to_plane_2'], ['W_CM_angle_to_plane_with_MET_1', 'W_CM_angle_to_plane_with_MET_2'],
            ['HNL_CM_mass_1', 'HNL_CM_mass_2'], 
            ['HNL_CM_mass_with_MET_1', 'HNL_CM_mass_with_MET_2'], 
            ['W_CM_angle_12','W_CM_angle_13', 'W_CM_angle_23', 'W_CM_angle_1MET', 'W_CM_angle_2MET', 'W_CM_angle_3MET']]

features.extend(['signal_label', 'channel', 'event_type', 'mass_hyp'])
flat_features = flatten_2D_list(features)

anatuple_path =  os.path.join('/eos/user/p/pdebryas/HNL/anatuple', period, tag)
channels = ['tee', 'tem', 'tmm', 'tte', 'ttm']

model_class = {}
scaler = {}
for parity in ['even', 'odd']:
    model_info_dir=os.path.join(os.getenv("RUN_PATH"), 'C_DNNvar', 'saved_models', f'{folderDNNName}_{parity}')

    #find the info of the model matching modelName
    model_info_list=pd.read_pickle(os.path.join(model_info_dir, 'model_info_list'))
    model_info = model_info_list[model_info_list["save_name"] == "/" + modelName]

    #load scaler
    scaler_path_orig = model_info['scaler_path'].iloc[0]
    scaler_path = os.path.join(model_info_dir,'scalers', scaler_path_orig.split('/')[-1])
    with open(scaler_path, 'rb') as f:
        scaler[parity] = pickle.load(f)

    #load selection of input variables
    input_variables = model_info['input_variables'].iloc[0]

    model_class[parity] = DNN_flexible(input_variables, model_info['hidden_layers'].iloc[0])
    model_class[parity].load_state_dict(load(os.path.join(model_info_dir, modelName + '.pt')))

All_channel_dict = process_channels(channels, flat_features, anatuple_path)

Mass_Hyp_vec = get_hnl_masses(period)

output_json_file = f'/eos/user/p/pdebryas/HNL/DNNscore/{period}/{folderDNNName}/'
os.makedirs(output_json_file, exist_ok=True)
for Mass_Hyp in Mass_Hyp_vec:
    print(f'...saving scores for mass {Mass_Hyp}')
    simple_dnndict = get_dnn_score_dict_torch_OddEven(All_channel_dict, model_class, input_variables, Mass_Hyp, scaler=scaler)
    output_filename = f'DNNscore_M{Mass_Hyp}.json'
    # Save the output dictionary as a JSON file
    with open(output_json_file + output_filename, "w") as json_file:
        json.dump(simple_dnndict, json_file, indent=4)