from uproot import open
import os
from fnmatch import filter
from copy import deepcopy
import vector
from vector import arr
import numpy as np
from functools import reduce
from operator import iconcat

from numpy import array, zeros_like, sqrt, cos, sin, arctan, exp, pi, unique, empty, concatenate, ones
from numpy.random import choice

# Global variables
output_vars_v5 = ['event', 'genWeight', 
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
                  ['HNL_CM_angle_with_MET_1', 'HNL_CM_angle_with_MET_2', 'HNL_CM_angle_with_MET_3'], 
                  ['W_CM_angle_to_plane_1', 'W_CM_angle_to_plane_2', 'W_CM_angle_to_plane_3'], ['W_CM_angle_to_plane_with_MET_1', 'W_CM_angle_to_plane_with_MET_2', 'W_CM_angle_to_plane_with_MET_3'],
                  ['HNL_CM_mass_1', 'HNL_CM_mass_2', 'HNL_CM_mass_3'], 
				  ['HNL_CM_mass_with_MET_1', 'HNL_CM_mass_with_MET_2', 'HNL_CM_mass_with_MET_3'], 
                  ['W_CM_angle_12','W_CM_angle_13', 'W_CM_angle_23', 'W_CM_angle_1MET', 'W_CM_angle_2MET', 'W_CM_angle_3MET'],
                  #'n_tauh',
                  ['px_1', 'py_1', 'pz_1', 'E_1', 'px_2', 'py_2', 'pz_2', 'E_2', 'px_3', 'py_3', 'pz_3', 'E_3'],
                  ['moth_mass_12', 'moth_mass_13', 'moth_mass_23', 'moth_pt_12', 'moth_pt_13', 'moth_pt_23', 'moth_eta_12', 'moth_eta_13', 'moth_eta_23', 'moth_phi_12', 'moth_phi_13', 'moth_phi_23', 'moth_px_12', 'moth_px_13', 'moth_px_23', 'moth_py_12', 'moth_py_13', 'moth_py_23', 'moth_pz_12', 'moth_pz_13', 'moth_pz_23', 'moth_E_12', 'moth_E_13', 'moth_E_23'],
                  'E_tot']


#===================================================================================================

class Data_extractor():
    """
    A Data_extractor extracts data from a folder of root files containing the anatuples.
    It takes a channel as argument : channel = "tee" "tem" "tmm" "tte" or "ttm"
    When called, it returns the variables of interest for the DNN training
    """
    def __init__(self, channel, raw_vars_general, raw_vars_lepton1, raw_vars_lepton2, raw_vars_lepton3, output_vars, functions, input_vars):
        """
        -channel : flavour of the 3 prompt leptons present in the decay. channel = "tee" "tem" "tmm" "tte" or "ttm"
        -raw_vars_general : names of variables in the root files that will be loaded and which are present only once, and not for each lepton
        -raw_vars_lepton(1,2,3) : end of names of variables in the root files that will be loaded and which are defined for a specific lepton.
                                The naming convention for such variables is L_X where L = Electron(1,2), Muon(1,2), Tau(1,2). Only specify
                                _X, since L will be deduced from the channel
        -output_vars : names of variable of interest that will be created by the data extractor
        -functions : functions that will be used to compute the output_vars (one function for each output_vars in the right order). If the 
                     corresponding output variable is already present as raw variable, put None as a function.
        -input_vars : list of lists of variables that are passed to the functions to compute the output_vars. If the variable in question 
                      is specific to one lepton, then "(1,2,3)_X" will be converted to lepton(1,2,3)_X. 
                      For example, in tee channel "3_mass"->"Electron2_mass"

        """
        self.channel = channel
        if self.channel == "tee":
            self.n_taus = 1
            self.lepton1 = "Tau"
            self.lepton2 = "Electron1"
            self.lepton3 = "Electron2"
        elif self.channel == "tem":
            self.n_taus = 1
            self.lepton1 = "Tau"
            self.lepton2 = "Electron"
            self.lepton3 = "Muon"
        elif self.channel == "tmm":
            self.n_taus = 1
            self.lepton1 = "Tau"
            self.lepton2 = "Muon1"
            self.lepton3 = "Muon2"
        elif self.channel == "tte":
            self.n_taus = 2
            self.lepton1 = "Tau1"
            self.lepton2 = "Tau2"
            self.lepton3 = "Electron"
        elif self.channel == "ttm":
            self.n_taus = 2
            self.lepton1 = "Tau1"
            self.lepton2 = "Tau2"
            self.lepton3 = "Muon"
        else:
            raise ValueError("The channel name \""+channel+"\" is not valid")
        self.raw_vars = raw_vars_general
        for var in raw_vars_lepton1:
            self.raw_vars.append(self.lepton1+var)
        for var in raw_vars_lepton2:
            self.raw_vars.append(self.lepton2+var)
        for var in raw_vars_lepton3:
            self.raw_vars.append(self.lepton3+var)
        
        self.input_vars = replace_prefix_in_list(input_vars, to_replace=['1','2','3'], replace_by=[self.lepton1, self.lepton2, self.lepton3])

        self.functions = functions
        self.output_vars = output_vars
        self.flat_output_vars = flatten_2D_list(output_vars)


    def __call__(self, path, signal_prefix = ['HNL'], real_data_prefix = ['EGamma', 'SingleMuon', 'Tau'], data = None, file_list = None, with_mass_hyp = True):
        """
        Arguments :
            -path : the path to the root files
            -signal_prefix : beginning of names of the files containing the signal (here "HNL"). It can be a string or a list of strings
            -real_data_prefix : beginning of filenames that correspond to real data, and that will be ignored
            -data : dictionnary to which the extracted data will be appended (if None, the dictionary will be created)
            -file_list : list of root files from which data will be extracted (if None, all root files present in path will be used).
            -with_mass_hyp : if True, the data will contain , the HNL mass hypothesis in GeV for the signal events, and a random choice 
                             among the different hypothesis for background events
        Output : 
            -data : dictionary containing the event indices, the variables of interest, the label of the event, and the type of event.
                    By default, data will contain the entries "signal_label" (1 for signal, 0 for background), "channel" and "event_type" (name of the 
                    file in which the events were taken)
        """
        total_keys = deepcopy(self.flat_output_vars)
        total_keys.extend(['signal_label', 'channel', 'event_type'])
        total_keys.append('mass_hyp')
        value_list = []
        for i in range(len(self.flat_output_vars)):
            value_list.append(empty((0,)))
        data = dict(zip(self.flat_output_vars, value_list))

        data['mass_hyp'] = []
        data['signal_label'] = []
        data['channel'] = []
        data['event_type'] = []


        if set(list(data.keys())) != set(total_keys):
            raise KeyError("The data keys don't match the names of the variable created by the data extractor : ", list(data.keys()), total_keys)

        if file_list == None:
            file_list = filter(os.listdir(path), '*.root')

        # Create a list of all considered HNL mass hypothesis
        if type(signal_prefix) != list:
                signal_prefix = [signal_prefix]

        mass_hyps = []
        if with_mass_hyp:
            for filename in file_list:
                for prefix in signal_prefix:
                    if filename[:len(prefix)] == prefix:
                        mass_hyps.append(isolate_int(filename, separators=['-', '_'])[0])
            mass_hyps = unique(array(mass_hyps))
        weightsum1=0
        weightsum2=0
        numsum2=0
        
        for filename in file_list:
            RealData = False
            for prefix in real_data_prefix:
                if filename[:len(prefix)] == prefix:
                    RealData = True

            with open(os.path.join(path,filename)) as file:
                if len(file.keys()) == 0:
                    print(f"The ROOT file {filename} is empty (no keys).")
                else:
                    anatuple = open(os.path.join(path,filename))['Events;1'].arrays(self.raw_vars, library='np')

                    n = len(anatuple[list(anatuple.keys())[0]])

                    if n==0:
                        continue

                    anatuple['channel'] = [self.channel]*n


                    # Creation of the data
                    for i, var in enumerate(self.output_vars):
                        if self.functions[i] == None:
                            data[var] = concatenate((data[var], anatuple[self.input_vars[i][0]]))
                        else:
                            outputs = self.functions[i](*call_dict_with_list(anatuple, self.input_vars[i]))
                            if type(var) == list:
                                for j,v in enumerate(var):
                                    data[v] = concatenate((data[v], outputs[j]))
                            else:
                                data[var] = concatenate((data[var], outputs))

                    label = 0
                    mass = ones((n,))
                    for prefix in signal_prefix:
                        if filename[:len(prefix)] == prefix:
                            label = 1
                            if with_mass_hyp:
                                mass *= isolate_int(filename,separators=['-', '_'])[0]
                    if label == 0 and with_mass_hyp:
                        mass = choice(mass_hyps, n)
                    
                    # Add mass hypothesis
                    if with_mass_hyp:
                        if 'mass_hyp' in data.keys():
                            data['mass_hyp'] = concatenate((data['mass_hyp'], mass))
                        else:
                            data['mass_hyp'] = mass
                    else:
                        if 'mass_hyp' in data.keys():
                            data['mass_hyp'] = concatenate((data['mass_hyp'], np.full((n,), None)))
                        else:
                            data['mass_hyp'] = np.full((n,), None)



                    # Add signal label (by default)
                    if 'signal_label' in data.keys():
                        data['signal_label'] = concatenate((data['signal_label'], ones((n,))*label))
                    else:
                        data['signal_label'] = ones((n,))*label

                    # Add channel (by default)
                    if 'channel' in data.keys():
                        data['channel'].extend([self.channel]*n)
                    else:
                        data['channel'] = [self.channel]*n

                    # Add event type (by default)
                    if 'event_type' in data.keys():
                        data['event_type'].extend([filename.replace('.root','')]*n)
                    else:
                        data['event_type'] = [filename.replace('.root','')]*n
 
        return data


'''
        for filename in file_list:
            RealData = False
            for prefix in real_data_prefix:
                if filename[:len(prefix)] == prefix:
                    RealData = True
           
            # # Raw data loading
            # limit_charge = 3
            # limit_tau_jet = 5
            # limit_em_iso = 0.15

            # cut = ''
            # if self.channel == 'tte':
            #     cut = '(abs(Tau1_charge + Tau2_charge + Electron_charge) < {}) & (Tau1_idDeepTau2018v2p5VSjet >= {}) & (Tau2_idDeepTau2018v2p5VSjet >= {}) & (Electron_pfRelIso03_all < {})'.format(limit_charge, limit_tau_jet, limit_tau_jet, limit_em_iso)

            # if self.channel == 'tee':
            #     cut = '(abs(Tau_charge + Electron1_charge + Electron2_charge) < {}) & (Tau_idDeepTau2018v2p5VSjet >= {}) & (Electron1_pfRelIso03_all < {}) & (Electron2_pfRelIso03_all < {})'.format(limit_charge, limit_tau_jet, limit_em_iso, limit_em_iso)

            # if self.channel == 'tem':
            #     cut = '(abs(Tau_charge + Electron_charge + Muon_charge) < {}) & (Tau_idDeepTau2018v2p5VSjet >= {}) & (Electron_pfRelIso03_all < {}) & (Muon_pfRelIso03_all < {})'.format(limit_charge, limit_tau_jet, limit_em_iso, limit_em_iso)

            # if self.channel == 'tmm':
            #     cut = '(abs(Tau_charge + Muon1_charge + Muon2_charge) < {}) & (Tau_idDeepTau2018v2p5VSjet >= {}) & (Muon1_pfRelIso03_all < {}) & (Muon2_pfRelIso03_all < {})'.format(limit_charge, limit_tau_jet, limit_em_iso, limit_em_iso)

            # if self.channel == 'ttm':
            #     cut = '(abs(Tau1_charge + Tau2_charge + Muon_charge) < {}) & (Tau1_idDeepTau2018v2p5VSjet >= {}) & (Tau2_idDeepTau2018v2p5VSjet >= {}) & (Muon_pfRelIso03_all < {})'.format(limit_charge, limit_tau_jet, limit_tau_jet, limit_em_iso)            
            
            # anatuple_before_cut = open(os.path.join(path,filename))['Events;1'].arrays(self.raw_vars, library='np') # type: ignore
            # weightsum_before_cut = anatuple_before_cut['genWeight'].sum()
            # weightsum1 += weightsum_before_cut
            # # print('weightsum before cut : ', weightsum_before_cut)
            
            anatuple = open(os.path.join(path,filename))['Events;1'].arrays(self.raw_vars, library='np') #open(os.path.join(path,filename))['Events;1'].arrays(self.raw_vars, cut=cut, library='np') # type: ignore
            #weightsum_after_cut = anatuple['genWeight'].sum()
            #weightsum2 += weightsum_after_cut
            #numsum2 += len(anatuple['genWeight'])

            n = len(anatuple[list(anatuple.keys())[0]])

            if n==0:
                continue

            anatuple['channel'] = [self.channel]*n


            # Creation of the data
            for i, var in enumerate(self.output_vars):
                if self.functions[i] == None:
                    data[var] = concatenate((data[var], anatuple[self.input_vars[i][0]]))
                else:
                    outputs = self.functions[i](*call_dict_with_list(anatuple, self.input_vars[i]))
                    if type(var) == list:
                        for j,v in enumerate(var):
                            data[v] = concatenate((data[v], outputs[j]))
                    else:
                        data[var] = concatenate((data[var], outputs))

            label = 0
            mass = ones((n,))
            for prefix in signal_prefix:
                if filename[:len(prefix)] == prefix:
                    label = 1
                    if with_mass_hyp:
                        mass *= isolate_int(filename,separators=['-', '_'])[0]
            if label == 0 and with_mass_hyp:
                mass = choice(mass_hyps, n)
            
            # Add mass hypothesis
            if with_mass_hyp:
                if 'mass_hyp' in data.keys():
                    data['mass_hyp'] = concatenate((data['mass_hyp'], mass))
                else:
                    data['mass_hyp'] = mass
            else:
                if 'mass_hyp' in data.keys():
                    data['mass_hyp'] = concatenate((data['mass_hyp'], np.full((n,), None)))
                else:
                    data['mass_hyp'] = np.full((n,), None)



            # Add signal label (by default)
            if 'signal_label' in data.keys():
                data['signal_label'] = concatenate((data['signal_label'], ones((n,))*label))
            else:
                data['signal_label'] = ones((n,))*label

            # Add channel (by default)
            if 'channel' in data.keys():
                data['channel'].extend([self.channel]*n)
            else:
                data['channel'] = [self.channel]*n

            # Add event type (by default)
            if 'event_type' in data.keys():
                data['event_type'].extend([filename.replace('.root','')]*n)
            else:
                data['event_type'] = [filename.replace('.root','')]*n
        
        # print('weightsum before cut : ', weightsum1)
        # print('weightsum after cut : ', weightsum2)
        # print('numsum after cut : ', numsum2)
        # weightsum= data['genWeight'].sum()
        # print("weightsum = ", weightsum)

        return data
'''

class Data_extractor_v5(Data_extractor):
    def __init__(self, channel):
        output_vars = deepcopy(output_vars_v5)
        functions =[None, None,                 # event, genWeight
                    None, None, None,           # charges
                    None, None, None, None,     # pts
                    None, None, None,           # etas
                    None, None, None,           # masses
                    None, None, None, None,         # phis
                    deltaphi, deltaphi, deltaphi, 
                    deltaphi, deltaphi, deltaphi,
                    deltaphi3,
                    deltaeta, deltaeta, deltaeta, 
                    deltaeta3,
                    deltaR, deltaR, deltaR, 
                    deltaR3,
                    sum_pt, 
                    transverse_mass, transverse_mass, transverse_mass, 
                    transverse_mass, transverse_mass, transverse_mass,
                    transverse_mass3,
                    invariant_mass, invariant_mass, invariant_mass,
                    invariant_mass,
                    total_transverse_mass, 
                    HNL_CM_angles_with_MET, 
                    W_CM_angles_to_plane, W_CM_angles_to_plane_with_MET,
			        HNL_CM_masses,
                    HNL_CM_masses_with_MET, 
                    W_CM_angles,
                    #count_tauh,
                    p4calc,
                    motherpair_vals,
                    Energy_tot]
        raw_vars_general = ['event', 'genWeight', 'MET_pt', 'MET_phi']
        lepton_specific = ['_eta', '_mass', '_phi', '_pt', '_charge'] #, '_genPartFlav']
        raw_vars_lepton1 = lepton_specific
        raw_vars_lepton2 = lepton_specific
        raw_vars_lepton3 = lepton_specific
        input_vars = [['event'], ['genWeight'], 
			        ['1_charge'], ['2_charge'], ['3_charge'], 
			        ['1_pt'], ['2_pt'], ['3_pt'], ['MET_pt'],
			        ['1_eta'], ['2_eta'], ['3_eta'], 
			        ['1_mass'], ['2_mass'], ['3_mass'],
                    ['1_phi'], ['2_phi'], ['3_phi'], ['MET_phi'],
			        ['1_phi', '2_phi'], ['1_phi', '3_phi'], ['2_phi', '3_phi'], 
			        ['1_phi', 'MET_phi'], ['2_phi', 'MET_phi'], ['3_phi', 'MET_phi'], 
			        ['1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        ['1_eta', '2_eta'], ['1_eta', '3_eta'], ['2_eta', '3_eta'], 
			        ['1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        ['1_eta', '2_eta', '1_phi', '2_phi'], ['1_eta', '3_eta', '1_phi', '3_phi'], ['2_eta', '3_eta', '2_phi', '3_phi'], 
			        ['1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        [['1_pt', '2_pt', '3_pt'],['1_phi', '2_phi', '3_phi'],['1_eta', '2_eta', '3_eta'], ['1_mass', '2_mass', '3_mass']], 
			        ['1_pt', '2_pt', '1_phi', '2_phi'], ['1_pt', '3_pt', '1_phi', '3_phi'], ['2_pt', '3_pt', '2_phi', '3_phi'],
			        ['1_pt', 'MET_pt', '1_phi', 'MET_phi'], ['2_pt', 'MET_pt', '2_phi', 'MET_phi'], ['3_pt', 'MET_pt', '3_phi', 'MET_phi'],
			        ['1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        [['1_pt', '2_pt'],['1_phi', '2_phi'],['1_eta', '2_eta'], ['1_mass', '2_mass']], [['1_pt', '3_pt'],['1_phi', '3_phi'],['1_eta', '3_eta'], ['1_mass', '3_mass']], [['2_pt', '3_pt'],['2_phi', '3_phi'],['2_eta', '3_eta'], ['2_mass', '3_mass']], 	
                    [['1_pt', '2_pt', '3_pt'],['1_phi', '2_phi', '3_phi'],['1_eta', '2_eta', '3_eta'], ['1_mass', '2_mass', '3_mass']], 
			        ['1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi'],
			        [ '1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        [ '1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'], ['1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        [ '1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'], 
			        [ '1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        ['1_pt', '2_pt', '3_pt', 'MET_pt', '1_phi', '2_phi', '3_phi', 'MET_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
			        #['channel', '1_genPartFlav', '2_genPartFlav', '3_genPartFlav'],
                    ['1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
                    ['1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass'],
                    ['1_pt', '2_pt', '3_pt', '1_phi', '2_phi', '3_phi', '1_eta', '2_eta', '3_eta', '1_mass', '2_mass', '3_mass', 'MET_pt']
                    ]
        super().__init__(channel, raw_vars_general=raw_vars_general, raw_vars_lepton1=raw_vars_lepton1, raw_vars_lepton2=raw_vars_lepton2, 
                         raw_vars_lepton3=raw_vars_lepton3, output_vars=output_vars, functions=functions, input_vars=input_vars)

def flatten_2D_list(multi_dim_list):
    new_list = []
    for ele in multi_dim_list:
        if type(ele) is list:
            new_list.append(ele)
        else:
            new_list.append([ele])
    return reduce(iconcat, new_list, [])

def count_tauh(channel, genPartFlavs_1, genPartFlavs_2, genPartFlavs_3):
    """
    Input : 
        -channel : string of three characters corresponding to the three prompt leptons in the decay
        -genPartFlavs : 3 (1 for each lepton) arguments describing the flavour of genParticle
    Output :
        -number of hadronic taus present in the event (either 0, 1 or 2) 
    """
    # if len(args) == 1:
    #     if len(args[0]) != 4:
    #         raise TypeError("Wrong number of arguments")
    #     channel = args[0][0][0]
    #     genPartFlavs = args[0][1:]
    # elif len(args) == 4:
    #     channel = args[0][0]
    #     genPartFlavs = args[1:]
    # else:
    #     raise TypeError("Wrong number of arguments")
    channel = channel[0]
    is_list = False
    genPartFlavs = [genPartFlavs_1, genPartFlavs_2, genPartFlavs_3]
    if type(genPartFlavs[0]) == list:
        is_list = True
        for lepton_flav in genPartFlavs:
            lepton_flav = np.array(lepton_flav)
    n_tauh = np.zeros_like(genPartFlavs[0]).astype('int64')
    for i, lepton_flav in enumerate(genPartFlavs):
        if channel[i] == 't':
            n_tauh += (lepton_flav==5).astype('int64')
    
    if is_list:
        n_tauh = n_tauh.tolist()
    
    return n_tauh 

def replace_prefix_in_list(list_, to_replace, replace_by):
    """
    Input :
        -list_ : python list of strings, potentially multidimensional
        -to_replace : list of characters or substrings that will be replaced in each element of the list
        -replace_by : list of characters or substrings that will replace the "to_replace" elements
    Output :
        -list with the same structure as the input list, with the replaced characters
    """
    if type(list_) != list:
        for i,s in enumerate(to_replace):
            if list_[:len(s)] == s:
                list_ = list_.replace(list_[:len(s)],replace_by[i])
        return list_
    else:
        sublist = []
        for el in list_:
            sublist.append(replace_prefix_in_list(el, to_replace, replace_by))
        return sublist
    
def isolate_int(string, separators):
    if type(separators) != list:
       separators = [separators]
    ints = []

    for i in range(1,len(separators)):
       string = string.replace(separators[i], separators[0])

    for z in string.split(separators[0]):
       if z.isdigit():
          ints.append(int(z))

    return ints

def call_dict_with_list(dictionary, list_):
    """
    Input :
        -python dictionary
        -python list (potentially multidimensional) of entries
    Output :
        -list with the same structure as the input list, but with the keys replaced by the values of the dictionary at the corresponding keys 
    """
    if type(list_) != list:
        return dictionary[list_]
    else:
        sublist = []
        for el in list_:
            sublist.append(call_dict_with_list(dictionary, el))
        return sublist

# kinematic functions ===================================================================================================

def eta_to_theta(eta):
    return 2 * arctan(exp(-eta))

def vec_3D(norm, phi, theta, hep_form=True):
    np_output = array(
        [norm * sin(theta) * cos(phi), norm * sin(theta) * sin(phi), norm * cos(theta)]
    )
    if hep_form:
        output = arr({"x": np_output[0], "y": np_output[1], "z": np_output[2]})
        return output
    return array

def deltaphi(phi1, phi2):
    """
    Arguments:
        -phi1 : azimuthal angle of the first particle
        -phi2 : azimuthal angle of the second particle
    """

    return abs((((phi1 - phi2) + pi) % (2 * pi)) - pi)

def deltaphi3(
    pt_1,
    pt_2,
    pt_3,
    pt_MET,
    phi_1,
    phi_2,
    phi_3,
    phi_MET,
    eta_1,
    eta_2,
    eta_3,
    mass_1,
    mass_2,
    mass_3,
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the leptons
        -pt_MET : norm of missing transverse momentum
        -phi_1,2,3 : azimuthal angle of the leptons
        -phi_MET : azimuthal angle of the missing transverse momentum
        -eta_1,2,3 : pseudorapidity of the leptons
        -eta_1,2,3 : mass of the leptons
    Output:
        -All combinations of deltaphi between 1 object and the sum of two others
    """
    n = len(pt_1)
    eta_MET = []
    mass_MET = []
    if type(pt_1) == list:
        for i in range(n):
            eta_MET.append(0)
            mass_MET.append(0)
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})
    vectors = [vector_1, vector_2, vector_3, vector_MET]

    groups = [
        [0, 1, 2],
        [1, 0, 2],
        [2, 0, 1],
        [3, 0, 1],
        [3, 0, 2],
        [3, 1, 2],
        [0, 1, 3],
        [0, 2, 3],
        [1, 0, 3],
        [1, 2, 3],
        [2, 0, 3],
        [2, 1, 3],
    ]
    deltaphis = []

    for comb in groups:
        deltaphis.append(
            deltaphi(vectors[comb[0]].phi, (vectors[comb[1]] + vectors[comb[2]]).phi)
        )

    return deltaphis

def deltaeta(eta1, eta2):
    """
    Arguments:
        -eta1 : pseudorapidity of the first particle
        -eta2 : pseudorapidity of the second particle
    """

    return abs(eta1 - eta2)

def deltaeta3(
    pt_1, pt_2, pt_3, phi_1, phi_2, phi_3, eta_1, eta_2, eta_3, mass_1, mass_2, mass_3
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the leptons
        -phi_1,2,3 : azimuthal angle of the leptons
        -eta_1,2,3 : pseudorapidity of the leptons
        -eta_1,2,3 : mass of the leptons
    Output:
        -All combinations of deltaeta between 1 lepton and the sum of two others
    """
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})
    vectors = [vector_1, vector_2, vector_3]

    groups = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
    deltaetas = []

    for comb in groups:
        deltaetas.append(
            deltaeta(vectors[comb[0]].eta, (vectors[comb[1]] + vectors[comb[2]]).eta)
        )

    return deltaetas

def deltaR(eta1, eta2, phi1, phi2):
    """
    Arguments:
        -eta1 : pseudorapidity of the first particle
        -eta2 : pseudorapidity of the second particle
        -phi1 : azimuthal angle of the first particle
        -phi2 : azimuthal angle of the second particle
    """

    return sqrt(deltaeta(eta1, eta2) ** 2 + deltaphi(phi1, phi2) ** 2)

def deltaR3(
    pt_1, pt_2, pt_3, phi_1, phi_2, phi_3, eta_1, eta_2, eta_3, mass_1, mass_2, mass_3
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the leptons
        -phi_1,2,3 : azimuthal angle of the leptons
        -eta_1,2,3 : pseudorapidity of the leptons
        -eta_1,2,3 : mass of the leptons
    Output:
        -All combinations of deltaR between 1 lepton and the sum of two others
    """
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})
    vectors = [vector_1, vector_2, vector_3]

    groups = [[0, 1, 2], [1, 0, 2], [2, 0, 1]]
    deltaRs = []

    for comb in groups:
        vec_sum = vectors[comb[1]] + vectors[comb[2]]
        deltaRs.append(
            deltaR(vectors[comb[0]].eta, vec_sum.eta, vectors[comb[0]].phi, vec_sum.phi)
        )

    return deltaRs

def sum_pt(pts, phis, etas, masses):
    """
    Aguments :
        -pts : transverse momentum of the particles
        -phis : azimuthal angles of the particles
        -etas : pseudorapidity of the particles
        -masses : masses of the particles
    All arguments have 2 coordinates :
        -the first component corresponds to the type of particle (muon, tau, MET)
        -The second coordinate corresponds to the event.
    Output :
        -Sum of the transverse momentum of the three leptons and the MET
    """

    p_tot = arr({"pt": pts[0], "phi": phis[0], "eta": etas[0], "M": masses[0]})
    for i in range(1, len(masses)):
        p_tot += arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
    return p_tot.pt

def transverse_mass(pt_1, pt_2, phi_1, phi_2):
    """
    Arguments :
        -pt_1 : transverse momentum of the first particle
        -pt_2 : transverse momentum of the second particle
        -phi_1 : azimuthal angle of the first particle
        -phi_2 : azimuthal angle of the second particle
    """
    
    result = 2.0 * pt_1 * pt_2 * (1.0 - np.cos(phi_1 - phi_2))

    # mask= result <0
    # mask_result=result[mask]
    # if len(mask_result)>0:
    #     print("mask_result", mask_result)
    
    # if result <0:
    #     print("pt1, pt2, phi1, phi2", pt_1, pt_2, phi_1, phi_2)
    # result[result < 0] = 0.0
    return np.sqrt(result)

def transverse_mass3(
    pt_1,
    pt_2,
    pt_3,
    pt_MET,
    phi_1,
    phi_2,
    phi_3,
    phi_MET,
    eta_1,
    eta_2,
    eta_3,
    mass_1,
    mass_2,
    mass_3,
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the leptons
        -pt_MET : norm of missing transverse momentum
        -phi_1,2,3 : azimuthal angle of the leptons
        -phi_MET : azimuthal angle of the missing transverse momentum
        -eta_1,2,3 : pseudorapidity of the leptons
        -eta_1,2,3 : mass of the leptons
    Output:
        -All combinations of transverse masses between 1 object and the sum of two others
    """
    n = len(pt_1)
    eta_MET = []
    mass_MET = []
    if type(pt_1) == list:
        for i in range(n):
            eta_MET.append(0)
            mass_MET.append(0)
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})
    vectors = [vector_1, vector_2, vector_3, vector_MET]

    groups = [
        [0, 1, 2],
        [1, 0, 2],
        [2, 0, 1],
        [3, 0, 1],
        [3, 0, 2],
        [3, 1, 2],
        [0, 1, 3],
        [0, 2, 3],
        [1, 0, 3],
        [1, 2, 3],
        [2, 0, 3],
        [2, 1, 3],
    ]
    transverse_masses = []

    for comb in groups:
        vec_sum = vectors[comb[0]] + vectors[comb[1]]
        transverse_masses.append(
            transverse_mass(
                vectors[comb[0]].pt, vec_sum.pt, vectors[comb[0]].phi, vec_sum.phi
            )
        )

    return transverse_masses

def total_transverse_mass(pt_1, pt_2, pt_3, pt_miss, phi_1, phi_2, phi_3, phi_miss):
    """
    Arguments :
        -pt_1 : transverse momentum of the first particle
        -pt_2 : transverse momentum of the second particle
        -pt_3 : transverse momentum of the third particle
        -pt_miss : missing transverse momentum
        -phi_1 : azimuthal angle of the first particle
        -phi_2 : azimuthal angle of the second particle
        -phi_3 : azimuthal angle of the third particle
        -phi_miss : azimuthal angle of missing particles
    """

    return sqrt(
        transverse_mass(pt_1, pt_2, phi_1, phi_2) ** 2
        + transverse_mass(pt_1, pt_3, phi_1, phi_3) ** 2
        + transverse_mass(pt_2, pt_3, phi_2, phi_3) ** 2
        + transverse_mass(pt_1, pt_miss, phi_1, phi_miss) ** 2
        + transverse_mass(pt_2, pt_miss, phi_2, phi_miss) ** 2
        + transverse_mass(pt_3, pt_miss, phi_3, phi_miss) ** 2
    )

def invariant_mass(pts, phis, etas, masses):
    """
    Aguments :
        -pts : transverse momentum of the particles
        -phis : azimuthal angles of the particles
        -etas : pseudorapidity of the particles
        -masses : masses of the particles
    All arguments have 2 coordinates :
        -the first component corresponds to the type of particle (muon, tau, MET)
        -The second coordinate corresponds to the event.
    Output : invariant mass of the sum of the 4-vectors of all objects passed to the function
    """
    p_tot = arr({"pt": pts[0], "phi": phis[0], "eta": etas[0], "M": masses[0]})
    for i in range(1, len(masses)):
        p_tot += arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
    return p_tot.mass

def HNL_CM_angles_with_MET(
    pt_1, pt_2, pt_3, pt_MET,
    phi_1, phi_2, phi_3, phi_MET,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
):
    """
    Arguments :
        -pt_1,2,3,MET : transverse momentum of the three leptons and the missing momentum
        -phi_1,2,3,MET : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[angle1, angle2, angle3] : angles between the 3 possible lepton pairs in their rest frame.
    """
    
    n = len(pt_1)

    if type(pt_1) == list:
        eta_MET = [0] * n
        mass_MET = [0] * n
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)

    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})

    # Lepton properties in dictionaries for easier looping
    pts = {1: pt_1, 2: pt_2, 3: pt_3}
    phis = {1: phi_1, 2: phi_2, 3: phi_3}
    etas = {1: eta_1, 2: eta_2, 3: eta_3}
    masses = {1: mass_1, 2: mass_2, 3: mass_3}

    # All possible pairs
    pair_candidate = [[1, 2], [1, 3], [2, 3]]
    angles = []

    for pair in pair_candidate:
        i, j = pair
        vector_i = arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
        vector_j = arr({"pt": pts[j], "phi": phis[j], "eta": etas[j], "M": masses[j]})

        vector_tot = vector_i + vector_j + vector_MET
        vector_i = vector_i.boostCM_of_p4(vector_tot)
        vector_j = vector_j.boostCM_of_p4(vector_tot)

        angle = vector_i.deltaangle(vector_j)
        angles.append(angle)

    return angles

def W_CM_angles_to_plane(
    pt_1, pt_2, pt_3,
    phi_1, phi_2, phi_3,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
):
    """
    Arguments :
        -pt_1,2,3 : transverse momentum of the three leptons
        -phi_1,2,3 : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[angle1, angle2, angle3] : angles between 1 lepton and the plane formed by the 2 other leptons in the rest frame 
                                    of the 3 leptons. There are 3 possible lepton pairs -> 3 angles in the output.
    """
    
    # Lepton properties in dictionaries for easier looping
    pts = {1: pt_1, 2: pt_2, 3: pt_3}
    phis = {1: phi_1, 2: phi_2, 3: phi_3}
    etas = {1: eta_1, 2: eta_2, 3: eta_3}
    masses = {1: mass_1, 2: mass_2, 3: mass_3}

    # All possible pairs
    pair_candidate = [[1, 2], [1, 3], [2, 3]]
    angles = []

    for pair in pair_candidate:
        i, j = pair
        k = next(value for value in [1, 2, 3] if value != i and value != j)

        vector_first = arr({"pt": pts[k], "phi": phis[k], "eta": etas[k], "M": masses[k]})
        vector_i = arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
        vector_j = arr({"pt": pts[j], "phi": phis[j], "eta": etas[j], "M": masses[j]})
        
        vector_tot = vector_i + vector_j + vector_first
        vector_i = vector_i.boostCM_of_p4(vector_tot)
        vector_j = vector_j.boostCM_of_p4(vector_tot)
        vector_first = vector_first.boostCM_of_p4(vector_tot)
        
        normal = vector_i.to_Vector3D().cross(vector_j.to_Vector3D())
        angle = vector_first.deltaangle(normal)
        
        angles.append(abs(pi / 2 - angle))
    return angles

def W_CM_angles_to_plane_with_MET(
    pt_1, pt_2, pt_3, pt_MET,
    phi_1, phi_2, phi_3, phi_MET,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
):
    """
    Arguments :
        -pt_1,2,3,MET : transverse momentum of the three leptons and the missing momentum
        -phi_1,2,3,MET : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[angle1, angle2, angle3] : angles between 1 lepton and the plane formed by the 2 other leptons in the rest frame 
                                    of the 3 leptons. There are 3 possible lepton pairs -> 3 angles in the output.
    """
    
    n = len(pt_1)
    if type(pt_1) == list:
        eta_MET = [0] * n
        mass_MET = [0] * n
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)
    
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})

    # Lepton properties in dictionaries for easier looping
    pts = {1: pt_1, 2: pt_2, 3: pt_3}
    phis = {1: phi_1, 2: phi_2, 3: phi_3}
    etas = {1: eta_1, 2: eta_2, 3: eta_3}
    masses = {1: mass_1, 2: mass_2, 3: mass_3}

    # All possible pairs
    pair_candidate = [[1, 2], [1, 3], [2, 3]]
    angles = []

    for pair in pair_candidate:
        i, j = pair
        k = next(value for value in [1, 2, 3] if value != i and value != j)

        vector_first = arr({"pt": pts[k], "phi": phis[k], "eta": etas[k], "M": masses[k]})
        vector_i = arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
        vector_j = arr({"pt": pts[j], "phi": phis[j], "eta": etas[j], "M": masses[j]})
        
        vector_tot = vector_i + vector_j + vector_first + vector_MET
        vector_i = vector_i.boostCM_of_p4(vector_tot)
        vector_j = vector_j.boostCM_of_p4(vector_tot)
        vector_first = vector_first.boostCM_of_p4(vector_tot)
        
        normal = vector_i.to_Vector3D().cross(vector_j.to_Vector3D())
        angle = vector_first.deltaangle(normal)
        
        angles.append(abs(pi / 2 - angle))
    return angles

def W_CM_angles(
    pt_1,
    pt_2,
    pt_3,
    pt_MET,
    phi_1,
    phi_2,
    phi_3,
    phi_MET,
    eta_1,
    eta_2,
    eta_3,
    mass_1,
    mass_2,
    mass_3,
):
    """
    Arguments :
        -pt_1,2,3,MET : transverse momentum of the three leptons and the missing momentum
        -phi_1,2,3,MET : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -Angle between the momenta of 2 objects in the center of mass frame of the first W boson (without considering missing momentum),
         for all 6 combination of 2 objects
    """

    n = len(pt_1)
    eta_MET = []
    mass_MET = []
    if type(pt_1) == list:
        for i in range(n):
            eta_MET.append(0)
            mass_MET.append(0)
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})

    vectors = [vector_1, vector_2, vector_3, vector_MET]

    vector_tot = vector_1 + vector_2 + vector_3

    for i in range(len(vectors)):
        vectors[i] = vectors[i].boostCM_of_p4(vector_tot)

    pairs = [[0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]]
    angles = []

    for pair in pairs:
        angle = vectors[pair[0]].deltaangle(vectors[pair[1]])
        angles.append(angle)

    return angles

def W_CM_angles_with_MET(
    pt_1,
    pt_2,
    pt_3,
    pt_MET,
    phi_1,
    phi_2,
    phi_3,
    phi_MET,
    eta_1,
    eta_2,
    eta_3,
    mass_1,
    mass_2,
    mass_3,
):
    """
    Arguments :
        -pt_1,2,3,MET : transverse momentum of the three leptons and the missing momentum
        -phi_1,2,3,MET : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -Angle between the momenta of 2 objects in the center of mass frame of the first W boson (with missing momentum),
         for all 6 combination of 2 objects
    """

    n = len(pt_1)
    eta_MET = []
    mass_MET = []
    if type(pt_1) == list:
        for i in range(n):
            eta_MET.append(0)
            mass_MET.append(0)
    else:
        eta_MET = zeros_like(pt_1)
        mass_MET = zeros_like(pt_1)
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})
    vector_1 = arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "M": mass_1})
    vector_2 = arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "M": mass_2})
    vector_3 = arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "M": mass_3})

    vectors = [vector_1, vector_2, vector_3, vector_MET]

    vector_tot = vector_1 + vector_2 + vector_3 + vector_MET

    for i in range(len(vectors)):
        vectors[i] = vectors[i].boostCM_of_p4(vector_tot)

    pairs = [[0, 1], [0, 2], [1, 2], [0, 3], [1, 3], [2, 3]]
    angles = []

    for pair in pairs:
        angle = vectors[pair[0]].deltaangle(vectors[pair[1]])
        angles.append(angle)

    return angles

def HNL_CM_masses(
    pt_1, pt_2, pt_3,
    phi_1, phi_2, phi_3,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
):
    """
    Arguments :
        -pt_1,2,3 : transverse momentum of the three leptons
        -phi_1,2,3 : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[HNL_mass1, HNL_mass2, HNL_mass3] : invariant mass of the sum of 2 leptons. 
                    There are 3 possible lepton pairs -> 3 masses in the output.
    """

    # Lepton properties in dictionaries for easier looping
    pts = {1: pt_1, 2: pt_2, 3: pt_3}
    phis = {1: phi_1, 2: phi_2, 3: phi_3}
    etas = {1: eta_1, 2: eta_2, 3: eta_3}
    masses = {1: mass_1, 2: mass_2, 3: mass_3}

    # All possible pairs
    pair_candidate = [[1, 2], [1, 3], [2, 3]]
    HNL_masses = []

    for pair in pair_candidate:
        i, j = pair

        vector_i = arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
        vector_j = arr({"pt": pts[j], "phi": phis[j], "eta": etas[j], "M": masses[j]})
        
        vector_tot = vector_i + vector_j
        HNL_masses.append(vector_tot.mass)
    return HNL_masses

def HNL_CM_masses_with_MET(
    pt_1, pt_2, pt_3,
    pt_MET,
    phi_1, phi_2, phi_3,
    phi_MET,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
):
    """
    Arguments :
        -pt_1,2,3,MET : transverse momentum of the three leptons and the missing momentum
        -phi_1,2,3,MET : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[HNL_mass1, HNL_mass2, HNL_mass3] : invariant mass of the sum of 2 leptons. 
                    There are 3 possible lepton pairs -> 3 masses in the output.
    """

    # MET properties
    n = len(pt_1)
    eta_MET = [0] * n
    mass_MET = [0] * n
    vector_MET = arr({"pt": pt_MET, "phi": phi_MET, "eta": eta_MET, "M": mass_MET})

    # Lepton properties in dictionaries for easier looping
    pts = {1: pt_1, 2: pt_2, 3: pt_3}
    phis = {1: phi_1, 2: phi_2, 3: phi_3}
    etas = {1: eta_1, 2: eta_2, 3: eta_3}
    masses = {1: mass_1, 2: mass_2, 3: mass_3}

    # All possible pairs
    pair_candidate = [[1, 2], [1, 3], [2, 3]]
    HNL_masses = []

    for pair in pair_candidate:
        i, j = pair

        vector_i = arr({"pt": pts[i], "phi": phis[i], "eta": etas[i], "M": masses[i]})
        vector_j = arr({"pt": pts[j], "phi": phis[j], "eta": etas[j], "M": masses[j]})
        
        vector_tot = vector_i + vector_j + vector_MET
        HNL_masses.append(vector_tot.mass)
    return HNL_masses

def p4calc(
    pt_1, pt_2, pt_3,
    phi_1, phi_2, phi_3,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3
):
    """
    Arguments :
        -pt_1,2,3 : transverse momentum of the three leptons
        -phi_1,2,3 : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output :
        -[px,py,pz,energy] : 4-momentum of each lepton. Since there are three,
        there are three lists
    """

    # Create vector arrays for each particle
    particle1 = vector.arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "mass": mass_1})
    particle2 = vector.arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "mass": mass_2})
    particle3 = vector.arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "mass": mass_3})

    # Extract 4-momentum components for each particle
    particles = [particle1, particle2, particle3]
    components = ['px', 'py', 'pz', 'energy']
    p4_list = []

    for particle in particles:
        for component in components:
            p4_list.append(getattr(particle, component))
    # print(p4_list)
    # assert(5==3)
    
    p4_list_flat = np.array(p4_list).flatten().tolist()
    return p4_list
    # return ([px_1, py_1, pz_1, energy_1, px_2, py_2, pz_2, energy_2, px_3, py_3, pz_3, energy_3])

def motherpair_vals(
    pt_1, pt_2, pt_3,
    phi_1, phi_2, phi_3,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the three leptons
        -phi_1,2,3 : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output:
        -add_feat_array : mass, pt, eta, phi, px, py, pz, energy of the mother particles
    """

    # Create vector arrays for each particle
    particle1 = vector.arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "mass": mass_1})
    particle2 = vector.arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "mass": mass_2})
    particle3 = vector.arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "mass": mass_3})

    # Create mother particles by summing up the particle vectors
    p4_mother12 = particle1 + particle2
    p4_mother23 = particle2 + particle3
    p4_mother13 = particle1 + particle3

    pairs = ['12', '13', '23']
    motherpairs = [p4_mother12, p4_mother13, p4_mother23]
    features_toadd = ['mass', 'pt', 'eta', 'phi', 'px', 'py', 'pz', 'energy']

    # Initialize the feature list
    features_list = []

    # Populate the feature list
    for feature in features_toadd:
        for i, pair in enumerate(pairs):
            features_list.append(getattr(motherpairs[i], feature))

    # features_list_flat = np.array(features_list).flatten().tolist()
    # print(features_list)
    return features_list

def Energychecker(input, output):
    input_np=input
    output=output
    inputvars2=['1_eta', '1_mass', '1_phi', '1_pt', '2_eta', '2_mass', '2_phi', '2_pt', '3_eta', '3_mass', '3_phi', '3_pt', 'MET_phi', 'MET_pt']
    inputdict={var: val for var, val in zip(inputvars2, input_np)}

    particle1 = vector.obj(phi=inputdict['1_phi'], pt=inputdict['1_pt'], eta=inputdict['1_eta'], mass=inputdict['1_mass'])
    particle2 = vector.obj(phi=inputdict['2_phi'], pt=inputdict['2_pt'], eta=inputdict['2_eta'], mass=inputdict['2_mass'])
    particle3 = vector.obj(phi=inputdict['3_phi'], pt=inputdict['3_pt'], eta=inputdict['3_eta'], mass=inputdict['3_mass'])

    Etot=particle1.E+particle2.E+particle3.E+inputdict['MET_pt']
    print("etot calculated vs predicted", Etot, output)

def Energy_tot(
    pt_1, pt_2, pt_3,
    phi_1, phi_2, phi_3,
    eta_1, eta_2, eta_3,
    mass_1, mass_2, mass_3,
    MET_pt
    
):
    """
    Arguments:
        -pt_1,2,3 : transverse momentum of the three leptons
        -phi_1,2,3 : azimuthal angle of the three leptons
        -eta_1,2,3 : pseudorapidity of the three leptons
        -mass_1,2,3 : mass of the three leptons
    Output:
        -Energy_tot : total energy of the three leptons
    """

    # Create vector arrays for each particle
    # particle1 = vector.arr({"pt": pt_1, "phi": phi_1, "eta": eta_1, "mass": mass_1})
    # particle2 = vector.arr({"pt": pt_2, "phi": phi_2, "eta": eta_2, "mass": mass_2})
    # particle3 = vector.arr({"pt": pt_3, "phi": phi_3, "eta": eta_3, "mass": mass_3})
    particle1=vector.array({"pt": pt_1, "phi": phi_1, "eta": eta_1, "mass": mass_1})
    particle2=vector.array({"pt": pt_2, "phi": phi_2, "eta": eta_2, "mass": mass_2})
    particle3=vector.array({"pt": pt_3, "phi": phi_3, "eta": eta_3, "mass": mass_3})

    # Create mother particles by summing up the particle vectors
    E_1=particle1.energy
    E_2=particle2.energy
    E_3=particle3.energy

    # Calculate the total energy
    Energy_tot = E_1+E_2+E_3+MET_pt

    return Energy_tot

