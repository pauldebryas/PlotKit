import uproot
import os
from tqdm import tqdm
import numpy as np
from hist_makers.regions.FakeRate_regions import compute_region_mask
import yaml

def dict_files(directory):
    files = {}
    for entry in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, entry)):
            files[entry.replace('_anatuple.root','')] = os.path.join(directory, entry)
    return files

def have_common_elements(list1, list2):
    return bool(set(list1) & set(list2))

def load_ntuples(samples, TreeName = 'Events;1'):
  processes = samples.keys()
  branches = {}
  for p in processes:
      file_path = samples[p]
      # Load the dataframes
      with uproot.open(file_path) as DataUproot:
        branches[p] = DataUproot['Events;1'].arrays()
  return branches

channels = ['tee','tmm','tem','tte','ttm']

EventsNb = {}
for channel in channels:
  EventsNb[channel] = {}
  print(f'processing channel {channel}')
  path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/2018/WithCorrections/{channel}/anatuple/'
  dictfiles = dict_files(path_to_files)
  for FileName in dictfiles.keys():
    if FileName[-5:] == '_2018':
        branches = load_ntuples({FileName:dictfiles[FileName]})
        cut_region = compute_region_mask(branches[FileName], channel, 'data', False)
        list = []
        for XXX in ['FFF','FFP','FPF','FPP','PFF','PFP','PPF']:
            list.append(((branches[FileName]['event'])[cut_region[f'AppRegion{XXX}']]).tolist())
        EventsNb[channel] = np.concatenate(list).tolist()

for channel1 in channels:
    for channel2 in channels:
        if channel1 == channel2:
            continue
        else:
            print(f'... comparing {channel1} and {channel2}')
            if have_common_elements(EventsNb[channel1], EventsNb[channel2]):
                print('true')


'''

N_events = {}
sumw_events = {}
for channel in channels:
    #print(f'For {channel}')
    N_events[channel] = {}
    sumw_events[channel] = {}
    for tag in tags:
        N_events[channel][tag] = {}
        sumw_events[channel][tag] = {}
        for HNLmass in inputs[channel][tag].keys():
            cut_region = compute_region_mask(inputs[channel][tag][HNLmass], channel, 'MC', IsSideband = False)
            cut_region['InvertedSignalRegion'] = ~cut_region['SignalRegion']
            N_events[channel][tag][HNLmass] = len(np.array(inputs[channel][tag][HNLmass]['genWeight']).flatten()[cut_region['SignalRegion']])
            sumw_events[channel][tag][HNLmass] = np.sum(np.array(inputs[channel][tag][HNLmass]['genWeight']).flatten()[cut_region['SignalRegion']])
print(N_events)    
print(sumw_events)

print(f'   HNL300:')
print(f'     Nevents:')    
print(f'      - Non Orth: {N_events[channel]["NotOrthoInAR"]["300"]}')
print(f'      - Orth: {N_events[channel]["AddOrthogonality"]["300"]}')
print(f'     Sumw:')    
print(f'      - Non Orth: {sumw_events[channel]["NotOrthoInAR"]["300"]}')
print(f'      - Orth: {sumw_events[channel]["AddOrthogonality"]["300"]}')
'''

'''
HNLMasses_x = [int(HNLmass.replace('HNL','')) for HNLmass in inputs['tee']["NotOrthoInAR"].keys()]

ratio_tte = []
ratio_tee = []
for HNLmass in HNLMasses_x:
    ratio_tte.append(N_events['tte']["AddOrthogonality"][f'HNL{str(HNLmass)}']/N_events['tte']["NotOrthoInAR"][f'HNL{str(HNLmass)}'] )
    ratio_tee.append(N_events['tee']["AddOrthogonality"][f'HNL{str(HNLmass)}']/N_events['tee']["NotOrthoInAR"][f'HNL{str(HNLmass)}'] )

# Plot the ratio
plt.figure(figsize=(10, 6))
plt.plot(HNLMasses_x, 100*np.array(ratio_tte), marker='o', linestyle='-', color='b', label='Ratio tte')
plt.plot(HNLMasses_x, 100*np.array(ratio_tee), marker='o', linestyle='-', color='r', label='Ratio tee')
plt.xlabel('HNL masses')
plt.ylabel('Ratio 100*(newSel/oldSel)')
plt.legend()
plt.grid(True)
plt.savefig('ratio_Anatuple.pdf')
#plt.show()
'''