import numpy as np
from matplotlib import pyplot as plt
import uproot
import os
from hep_ml.metrics_utils import ks_2samp_weighted
from statsmodels.stats.proportion import proportion_confint

def compute_ratio_witherr(N_num, N_den):
    ''' Return the ratio between N_num and N_den with err computed using Clopper-Pearson interval based on Beta distribution with alpha=0.32
    '''
    conf_int = proportion_confint(count=N_num, nobs=N_den, alpha=0.32, method='beta')
    ratio = N_num / N_den
    ratio_low = ratio - conf_int[0]
    ratio_up = conf_int[1] - ratio
    return ratio, ratio_low, ratio_up

def load_input_files(inputs_cfg, input_dir):
    input_files = {}
    for input in inputs_cfg:
        if 'files' not in input.keys():
            raise f'Missing files dict for data in {input.keys()}'
        files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
        for file in files_list:
            if os.path.isfile(file) == False:
                print('WARNING: ' + file + ' is missing')
                files_list.remove(file)
        input_files[input['name']] = files_list
    return input_files

def load_ntuples(files):
    list = []
    for file_path in files:
        # Load the dataframes
        DataUproot = uproot.open(file_path)
        list.append(DataUproot['Event;1'].arrays())
    return np.concatenate(list)

def draw_distributions(original, target, new_original_weights, name_fig, columns, hist_settings, savefig_dir):
    if not os.path.exists(f'{savefig_dir}/{name_fig}/'):
        os.makedirs(f'{savefig_dir}/{name_fig}/')
    
    for id, column in enumerate(columns, 1):
        plt.figure()
        xlim = np.percentile(np.hstack([target[column]]), [0.01, 99.99])
        plt.hist(original[column], weights=new_original_weights, range=xlim, **hist_settings[column])
        plt.hist(target[column], range=xlim, **hist_settings[column])
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(original[column], target[column], weights1=new_original_weights, weights2=np.ones(len(target), dtype=float)))
        plt.savefig(f'{savefig_dir}/{name_fig}/{name_fig}_{column}.pdf')
        
def check_ks_of_expression(expression, original_test, target_test, original_weights_test, bins_weights_test, gb_weights_test):
    col_original = original_test.eval(expression, engine='python')
    col_target = target_test.eval(expression, engine='python')
    w_target = np.ones(len(col_target), dtype='float')
    print('No reweight   KS:', ks_2samp_weighted(col_original, col_target, weights1=original_weights_test, weights2=w_target))
    print('Bins reweight KS:', ks_2samp_weighted(col_original, col_target, weights1=bins_weights_test, weights2=w_target))
    print('GB Reweight   KS:', ks_2samp_weighted(col_original, col_target, weights1=gb_weights_test, weights2=w_target))