import os
import yaml
import numpy as np
import matplotlib.pyplot as plt

from common.helpers import load_ntuples, dict_files
from common.regions.regions import compute_region_mask

# parameters -----------------------------------------------------------------------------------------------------
period = '2017'
tag = 'FinalProd'
flavor = 'Muon'  #"Muon" or "Electron" or "Tau"
region_name = 'DYRegionFR' # "DYRegionFR" or "ttbarRegionFR"
# ----------------------------------------------------------------------------------------------------------------

if flavor == "Muon":
    if region_name == "DYRegionFR":
        channels = ["Zmu"]
    if region_name == "ttbarRegionFR":
        channels = ["llmu"]
if flavor == "Electron":
    if region_name == "DYRegionFR":
        channels = ["Ze"]
    if region_name == "ttbarRegionFR":
        channels = ["lle"]
if flavor == "Tau":
    if region_name == "DYRegionFR":
        channels = ["tll"]
    if region_name == "ttbarRegionFR":
        channels = ["tem"]

output_path_fig = os.path.join(os.getenv("RUN_PATH", "."), f'A_cutflow/figures/{tag}/{period}/FakeLeptonOrigin/')
os.makedirs(output_path_fig, exist_ok=True)

output_path_fig_file = os.path.join(output_path_fig, f'MCin{region_name}_FakeLepton_{flavor}_origin')

# color map for labels
label_color_map = {
    'lightquark': '#007de9',
    'heavyquark': '#8567ad',
    'gluon': '#ff7a7a',
    'others': '#717275',
}

my_dict = {}
for channel in channels:
    my_dict[channel] = {}
    print(f'processing channel {channel}')
    if flavor in ['Electron', 'Muon']:
        path_to_files = f'/eos/user/p/pdebryas/HNL_LLFF/anatuple/{period}/{tag}/{channel}/anatuple/'
    else:
        path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/{channel}/anatuple/'
    dictfiles = dict_files(path_to_files)
    branches = load_ntuples(dictfiles)

    for FileName in branches.keys():
        # IMPORTANT FIXES:
        # - use Python 'not' and 'and' for string checks instead of bitwise ~ and &
        # - skip files that end with _{period} or start with HNL_ or have no branches
        if (not FileName.endswith(f'_{period}')) and (FileName[:4] != 'HNL_') and (len(branches[FileName]) != 0):
            #print(f'  -> file: {FileName}')
            my_dict[channel][FileName] = {}

            # compute region mask (returns a dict of boolean masks per region)
            cut_region = compute_region_mask(branches[FileName], channel, 'MC', region_name)
            # the mask selecting the "NotTrueLeptons" entries in the chosen region
            mask_key = f'{region_name}_PassLooseNotTightWP_{flavor}IsNotPromptLepton'
            if mask_key not in cut_region:
                # defensive: skip if region mask not present
                print(f'     WARNING: mask {mask_key} not in compute_region_mask result for {FileName}. Skipping.')
                continue
            base_mask = cut_region[mask_key]

            #------
            ValidJet_mask = np.concatenate(np.abs((branches[FileName][f'{flavor}_MatchingJetValid']).tolist())) == 1
            base_mask = base_mask & ValidJet_mask
            #------

            # build masks in the full array, restricted by base_mask only once
            full_flavour_array = np.concatenate(np.abs((branches[FileName][f'{flavor}_MatchingJetPartonFlavour']).tolist()))
            flavour_in_base = full_flavour_array[base_mask]

            genw_base = np.array((branches[FileName]['genWeight'][base_mask]).tolist())

            # now build masks on the restricted flavour array
            cut_LLdquark = np.array(flavour_in_base == 1)
            cut_LLuquark = np.array(flavour_in_base == 2)
            cut_LLsquark = np.array(flavour_in_base == 3)
            cut_LLcquark = np.array(flavour_in_base == 4)
            cut_LLbquark = np.array(flavour_in_base == 5)
            cut_LLtquark = np.array(flavour_in_base == 6)
            cut_LLgluon =  np.array(flavour_in_base == 21)
            cut_LLothers = ~(cut_LLdquark | cut_LLuquark | cut_LLsquark | cut_LLcquark | cut_LLbquark | cut_LLtquark | cut_LLgluon)

            my_dict[channel][FileName]['sumw'] = np.sum(genw_base)
            my_dict[channel][FileName]['sumw_lightquark'] = np.sum(genw_base[cut_LLdquark | cut_LLuquark | cut_LLsquark])
            my_dict[channel][FileName]['sumw_heavyquark'] = np.sum(genw_base[cut_LLcquark | cut_LLbquark | cut_LLtquark])
            my_dict[channel][FileName]['sumw_gluon']     =  np.sum(genw_base[cut_LLgluon])
            my_dict[channel][FileName]['sumw_others']     = np.sum(genw_base[cut_LLothers])
#print(my_dict)

# Helper to safe-initialize sum dicts
def zero_group_dict():
    return {
        'sumw': 0.0,
        'sumw_lightquark': 0.0,
        'sumw_heavyquark': 0.0,
        'sumw_gluon': 0.0,
        'sumw_others': 0.0,
    }

def plot_pie_from_dict(d, title, outname):
    # d must have the 4 keys below
    labels = ['lightquark', 'heavyquark', 'gluon', 'others']
    sizes = [d.get('sumw_lightquark', 0.0),
             d.get('sumw_heavyquark', 0.0),
             d.get('sumw_gluon', 0.0),
             d.get('sumw_others', 0.0)]

    # if everything is zero, skip plotting
    if sum(sizes) == 0:
        print(f'   -> nothing to plot for {title} (all zeros).')
        return
    if any(s < 0 for s in sizes):
        print(f'   -> nothing to plot for {title} (negative values detected).')
        return
    colors = [label_color_map[l] for l in labels]
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close(fig)
    print(f'   -> saved {outname}')

# Now build and plot per channel
for channel in channels:
    print(f'... plotting channel {channel}')
    if not my_dict.get(channel):
        print(f'   WARNING: no files collected for channel {channel}.')
        continue

    # Initialize accumulators
    my_All = zero_group_dict()
    my_DY = zero_group_dict()
    # my_ST = zero_group_dict()
    my_TT = zero_group_dict()
    # my_Wjets = zero_group_dict()

    # Accumulate per-file into the right category
    for MCfile, values in my_dict[channel].items():
        # skip files that didn't get processed fully
        if 'sumw' not in values:
            continue
        # All (everything)
        for k in my_All.keys():
            my_All[k] += values.get(k, 0.0)

        # Identify groups by filename patterns
        if MCfile.startswith('DYJetsToLL_'):
            for k in my_DY.keys():
                my_DY[k] += values.get(k, 0.0)
            continue

        # if MCfile.endswith('JetsToLNu') or MCfile.startswith('WJetsToLNu_HT'):
        #     for k in my_Wjets.keys():
        #         my_Wjets[k] += values.get(k, 0.0)
        #     continue

        # if MCfile.startswith('ST_'):
        #     for k in my_ST.keys():
        #         my_ST[k] += values.get(k, 0.0)
        #     continue

        if MCfile.startswith('TTTo') or MCfile.startswith('TTJets') or MCfile.startswith('TT_'):
            for k in my_TT.keys():
                my_TT[k] += values.get(k, 0.0)
            continue

        # If none of the above groups matched, they are still included in All (already added).

    # Produce pies: All, DY, TT, ST, Wjets (only if they have content)
    print(my_All)
    plot_pie_from_dict(my_All, f'{channel} - All MC (origin of fake {flavor})', output_path_fig_file + f'_{channel}_All.pdf')
    if (my_DY['sumw'] > 0) and (region_name == "DYRegionFR"):
        plot_pie_from_dict(my_DY, f'{channel} - DY (origin of fake {flavor})', output_path_fig_file + f'_{channel}_DY.pdf')
    if (my_TT['sumw'] > 0) and (region_name == "ttbarRegionFR"):
        plot_pie_from_dict(my_TT, f'{channel} - TT (origin of fake {flavor})', output_path_fig_file + f'_{channel}_TT.pdf')
    # if my_ST['sumw'] > 0:
    #     plot_pie_from_dict(my_ST, f'{channel} - ST (origin of fake {flavor})', output_path_fig_file + f'_{channel}_ST.pdf')
    # if my_Wjets['sumw'] > 0:
    #     plot_pie_from_dict(my_Wjets, f'{channel} - Wjets (origin of fake {flavor})', output_path_fig_file + f'_{channel}_Wjets.pdf')
