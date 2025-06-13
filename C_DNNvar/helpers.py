from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import torch

from C_DNNvar.data_extractor_class import Data_extractor_v5

def process_channels(channels, flat_features, path):
    """
    This function processes a list of channels and returns a dictionary. 
    For each channel, it creates a sub-dictionary that contains n dictionaries
    where n is the number of files in path. Each of these dictionaries contains 
    40 features as keys and corresponding numpy arrays as values.
    
    The resulting dictionary structure looks like:
    
    channels_data = {
        'channel1': {
            'filename1': {feature1: numpy array, feature2: numpy array, ..., feature40: numpy array},
            'filename2': {feature1: numpy array, feature2: numpy array, ..., feature40: numpy array},
            ...
        },
        'channel2': {
            'filename1': {feature1: numpy array, feature2: numpy array, ..., feature40: numpy array},
            'filename2': {feature1: numpy array, feature2: numpy array, ..., feature40: numpy array},
            ...
        },
        ...
    }
    
    Args:
    channels (list): List of channels to be processed.
    flat_features (list): List of features to be extracted.
    path (str): Path to the anatuple directory.
    
    Returns:
    dict: A dictionary with the described structure.
    """

    # produce data var needed by extractor
    values = []
    for i in range(len(flat_features)):
        values.append([])
    data = dict(zip(flat_features, values))

    All_channel_dict = {}
    pbar = tqdm(channels)

    for channel in pbar:
        pbar.set_description(f"Processing {channel}")
        All_channel_dict[channel] = {}

        # Initialize data extractor
        extractor = Data_extractor_v5(channel)

        # Extract data
        path_anatuple = os.path.join(path, channel, 'anatuple')
        data_dict4 = extractor(path_anatuple, data=data, with_mass_hyp = False)

        filelist = os.listdir(path_anatuple)
        for file in filelist:
            filename = file.replace('.root','')
            mask_file = (np.array(data_dict4['event_type']) == filename)
            print(f'{np.sum(mask_file)} events in {filename}')
            file_data = {}
            for feature in flat_features:
                # Skip unwanted features
                if feature == 'signal_labels' or feature == 'event_type':
                    continue

                # Populate file_data for the feature
                file_data[feature] = np.array(data_dict4[feature])[mask_file].flatten()

            # Add the background and signal data to the all channel dictionary
            All_channel_dict[channel][filename.replace('_anatuple','')] = file_data

    return All_channel_dict

def get_dnn_score_dict_torch_OddEven(data_dict_dnn, model_class, vars_list, masshyp, scaler=None):

    model_even =model_class['even']
    model_even.eval()

    model_odd =model_class['odd']
    model_odd.eval()

    output_dict = {}
    vars_list_copy= vars_list.copy()
    vars_list_copy.remove('signal_label')
    vars_list_copy.remove('weightNorm')

    for channel in  tqdm(data_dict_dnn.keys(), desc='channel', disable=True):
        output_dict[channel] = {}

        for filename in  tqdm(data_dict_dnn[channel].keys(), desc='filename', disable=True):
            output_dict[channel][filename] = {}
        
            data = pd.DataFrame.from_dict(data_dict_dnn[channel][filename])
            data['mass_hyp'] = masshyp

            even_df = data[data['event'] % 2 == 0]
            odd_df = data[data['event'] % 2 == 1]

            even_data_all_concat  = even_df[vars_list_copy]
            even_data_all_concat  = even_data_all_concat.to_numpy()

            odd_data_all_concat  = odd_df[vars_list_copy]
            odd_data_all_concat  = odd_data_all_concat.to_numpy()

            if len(even_data_all_concat) == 0:
                output_even_event  = []
                output_even_scores = []
            else:
                if scaler['even'] is not None:
                    even_data_all_concat=scaler['even'].transform(even_data_all_concat)

                even_data_all_concat_tensor = torch.tensor(even_data_all_concat).float()

                with torch.no_grad():
                    output=model_odd(even_data_all_concat_tensor)
                scores=output.numpy()
                
                output_even_event  = even_df['event'].to_numpy()
                output_even_scores = np.array(scores[:len(even_df)].flatten())

            if len(odd_data_all_concat) == 0:
                output_odd_event  = []
                output_odd_scores = []
            else:
                if scaler['odd'] is not None:
                    odd_data_all_concat=scaler['odd'].transform(odd_data_all_concat)

                odd_data_all_concat_tensor = torch.tensor(odd_data_all_concat).float()
                
                with torch.no_grad():
                    output=model_even(odd_data_all_concat_tensor)
                scores=output.numpy()
                
                output_odd_event  = odd_df['event'].to_numpy()
                output_odd_scores = np.array(scores[:len(odd_df)].flatten())

            output_dict[channel][filename]['event'] = np.concatenate([output_even_event, output_odd_event]).tolist()
            output_dict[channel][filename]['scores']= np.concatenate([output_even_scores, output_odd_scores]).tolist()

    return output_dict

class DNN_flexible(torch.nn.Module):
    def __init__(self, input_vars, hidden_layer_sizes):
        super(DNN_flexible, self).__init__()

        layer_sizes = [len(input_vars)-2] + hidden_layer_sizes + [1]

        self.layers = torch.nn.ModuleList()
        for i in range(len(layer_sizes) - 1):
            self.layers.append(torch.nn.Linear(layer_sizes[i], layer_sizes[i + 1]))
        
        self.relu = torch.nn.ReLU()
        self.sigmoid = torch.nn.Sigmoid()

    def forward(self, x):
        for i in range(len(self.layers) - 1):
            x = self.relu(self.layers[i](x))
        x = self.sigmoid(self.layers[-1](x))
        return x
    

def get_dnn_score_dict_torch_simple(data_dict_dnn, model_class, vars_list, masshyp, scaler=None):
    """
    Given a dictionary containing data for various channels and a deep learning model, this method computes the model's
    scores for each event. The method also modifies the original dictionary to 
    include these scores.

    Parameters:
    data_dict_dnn (dict): A dictionary where each key is a channel and its corresponding value is a nested dictionary 
                          with keys filename, each associated with a DataFrame containing features for 
                          each instance in the corresponding category. 
    model_class (str): trained model score for score calculation.
    vars_list (list): A list of feature names that the model uses for prediction. 
                      It should include 'signal_label' and 'weightNorm' but they will be removed inside the function.
    masshyp (float): The mass hypothesis under consideration.
    scaler (object, optional): An instance of a preprocessing scaler if the data needs to be scaled. 
                               The default is None, indicating that no scaling is required.

    Returns:
    dict: The modified dictionary where each key DataFrame now includes a new 'scores' column 
          containing the model's scores for each event.
    """

    model_even =model_class['even']
    model_even.eval()

    model_odd =model_class['odd']
    model_odd.eval()

    output_dict = {}
    vars_list_copy= vars_list.copy()
    vars_list_copy.remove('signal_label')
    vars_list_copy.remove('weightNorm')

    for channel in  tqdm(data_dict_dnn.keys(), desc='channel', disable=True):
        output_dict[channel] = {}

        for filename in  tqdm(data_dict_dnn[channel].keys(), desc='filename', disable=True):
            output_dict[channel][filename] = {}
        
            data = pd.DataFrame.from_dict(data_dict_dnn[channel][filename])
            data['mass_hyp'] = masshyp

            even_df = data[data['event'] % 2 == 0]
            odd_df = data[data['event'] % 2 == 1]

            even_data_all_concat  = even_df[vars_list_copy]
            even_data_all_concat  = even_data_all_concat.to_numpy()

            odd_data_all_concat  = odd_df[vars_list_copy]
            odd_data_all_concat  = odd_data_all_concat.to_numpy()

            if len(even_data_all_concat) == 0:
                output_even_event  = []
                output_even_scores = []
            else:
                if scaler['even'] is not None:
                    even_data_all_concat=scaler['even'].transform(even_data_all_concat)

                even_data_all_concat_tensor = torch.tensor(even_data_all_concat).float()

                with torch.no_grad():
                    output=model_odd(even_data_all_concat_tensor)
                scores=output.numpy()
                
                output_even_event  = (even_df['event'].to_numpy()).tolist()
                output_even_scores = np.array(scores[:len(even_df)].flatten()).tolist()

            if len(odd_data_all_concat) == 0:
                output_odd_event  = []
                output_odd_scores = []
            else:
                if scaler['odd'] is not None:
                    odd_data_all_concat=scaler['odd'].transform(odd_data_all_concat)

                odd_data_all_concat_tensor = torch.tensor(odd_data_all_concat).float()
                
                with torch.no_grad():
                    output=model_even(odd_data_all_concat_tensor)
                scores=output.numpy()
                
                output_odd_event  = (odd_df['event'].to_numpy()).tolist()
                output_odd_scores = np.array(scores[:len(odd_df)].flatten()).tolist()

            output_dict[channel][filename]['event'] = np.concatenate([output_even_event, output_odd_event])
            output_dict[channel][filename]['scores']= np.concatenate([output_even_scores, output_odd_scores])

    return output_dict