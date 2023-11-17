import numpy as np
import os
import json
from matplotlib import pyplot as plt

# comparison figure
save_fig = f'{os.getenv("RUN_PATH")}/GBReweighter/figures/weights_distribution.pdf'

weights_dict = {}
for channel in ['tee','tmm']:
    #load the weights
    with open(f'{os.getenv("RUN_PATH")}/GBReweighter/results/gb_weights_ControlRegion_{channel}.json', 'r') as file:
        weights = json.load(file)
    
    weights_dict[channel] = np.array(weights['gb_weights'])


plt.figure()
for channel in ['tee','tmm']:
    plt.hist(weights_dict[channel], label=f'{channel} weights',bins = 20, alpha= 0.5, density=True) 
plt.legend()
plt.savefig(save_fig)