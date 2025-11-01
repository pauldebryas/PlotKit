import yaml
import matplotlib.pyplot as plt

from common.helpers import load_ntuples, dict_files
from common.regions.regions import compute_region_mask

def save_sumwEvents(channels, tag, period, output_file, region_name):
    my_dict_PassTightWP = {}
    my_dict_PassTightWP_TrueLeptons = {}
    for channel in channels:
        my_dict_PassTightWP[channel] = {}
        my_dict_PassTightWP_TrueLeptons[channel] = {}
        print(f'processing channel {channel}')
        path_to_files = f'/eos/user/p/pdebryas/HNL/anatuple/{period}/{tag}/{channel}/anatuple/'
        dictfiles = dict_files(path_to_files)
        branches = load_ntuples(dictfiles)
        for FileName in branches.keys():
            if ~FileName.endswith(f'_{period}') & (FileName[:4] != 'HNL_') & (len(branches[FileName]) != 0):
                cut_region = compute_region_mask(branches[FileName], channel, 'MC', region_name)
                my_dict_PassTightWP[channel][FileName] = sum(((branches[FileName]['genWeight'])[cut_region[f'{region_name}_PassTightWP']]).tolist())
                my_dict_PassTightWP_TrueLeptons[channel][FileName] = sum(((branches[FileName]['genWeight'])[cut_region[f'{region_name}_PassTightWP_TrueLeptons']]).tolist())

    my_dict = {}
    my_dict['PassTightWP'] = my_dict_PassTightWP
    my_dict['PassTightWP_TrueLeptons'] = my_dict_PassTightWP_TrueLeptons

    with open(output_file, 'w') as outfile:
        yaml.dump(my_dict, outfile)
    return my_dict

def plot_pieFig(channels, output_file, output_path_fig, onlyTL = False):

    # Define a dictionary that maps each label to a specific color
    label_color_map = {
        'WZ': '#007de9',
        'WminusHToTauTau': '#8567ad',
        'WplusHToTauTau': '#ff7a7a',
        'ZZ': '#f68e2f',
        'DYJetsToLL':'#1daab4',
        'TT': '#238d23',
        'WWW': '#e6ce06',
        'WJetsToLNu': '#4d2e07',
        'WW': '#7fff00',
        'ZHToTauTau': '#b2ffff',
        'TTWW': '#fed0ee',
        'TTZToLLNuNu_M-10': '#2a52be',
        'WWZ': '#e6ce06',
        'ttHToTauTau': '#ffe4e1', 
        'ST': '#fc0fc0', 
        'EWK': '#fc330f', 
        'others': '#717275',
    }

    # Read YAML file
    with open(output_file, 'r') as stream:
        SumwFile = yaml.safe_load(stream)

    if onlyTL:
        MCsamples = SumwFile['PassTightWP_TrueLeptons']
    else:
        MCsamples = SumwFile['PassTightWP']

    for channel in channels:
        print(f'... plotting channel {channel}')
        sumwtotal = sum([MCsamples[channel][MCfile] for MCfile in MCsamples[channel].keys()])
        MCsamples_grouped = {}
        MCsamples_grouped['DYJetsToLL'] = 0
        MCsamples_grouped['EWK'] = 0
        MCsamples_grouped['ST'] = 0
        MCsamples_grouped['TT'] = 0
        MCsamples_grouped['WJetsToLNu'] = 0
        MCsamples_grouped['others'] = 0
        for MCfile in MCsamples[channel].keys():
            if MCfile.startswith('DYJetsToLL_'):
                MCsamples_grouped['DYJetsToLL'] += MCsamples[channel][MCfile]
                continue
            if MCfile.startswith('EWK'):
                MCsamples_grouped['EWK'] += MCsamples[channel][MCfile]
                continue
            if MCfile.startswith('ST_'):
                MCsamples_grouped['ST'] += MCsamples[channel][MCfile]
                continue
            if MCfile.startswith('TTTo'):
                MCsamples_grouped['TT'] += MCsamples[channel][MCfile]
                continue
            if MCfile.endswith('JetsToLNu') | MCfile.startswith('WJetsToLNu_HT'):
                MCsamples_grouped['WJetsToLNu'] += MCsamples[channel][MCfile]
                continue
            if (100*MCsamples[channel][MCfile]/sumwtotal) <= 2.4:
                MCsamples_grouped['others'] += MCsamples[channel][MCfile]
                continue
            MCsamples_grouped[MCfile] = MCsamples[channel][MCfile]

        key_del = []
        for mykey in MCsamples_grouped.keys():
            if (100*MCsamples_grouped[mykey]/sumwtotal) <= 2.4:
                MCsamples_grouped['others'] += MCsamples_grouped[mykey]
                key_del.append(mykey)
                continue

        for kes in key_del:
            del MCsamples_grouped[kes]
 
        labels = MCsamples_grouped.keys()
        colors = [label_color_map[label] for label in labels]
        sizes = [MCsamples_grouped[MCfile] for MCfile in MCsamples_grouped.keys()] 
        fig, ax = plt.subplots()
        #plt.title(f'MC in Signal Region in {channel}')
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors)
        plt.savefig(output_path_fig + f'_{channel}.pdf')
    return

def save_info_XsecUnc(channels, input_file, output_file):

    # Read YAML file
    with open(input_file, 'r') as stream:
        SumwFile = yaml.safe_load(stream)
    MCsamples = SumwFile['PassTightWP_TrueLeptons']

    XsecFile = '/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/AnatupleProd/config/crossSections13TeV.yaml'
    with open(XsecFile, 'r') as f:
        Xsec = yaml.safe_load(f)

    my_dict = {}
    for channel in channels:
        my_dict[channel] = {}
        print(f'------ channel {channel} ------')
        sumwtotal = sum([MCsamples[channel][MCfile] for MCfile in MCsamples[channel].keys()])
        impactMCsamples_sorted = dict(sorted(MCsamples[channel].items(), key=lambda item: item[1], reverse=True))
        integrated_sum = 0
        for MC in impactMCsamples_sorted.keys():
            my_dict[channel][MC] = {}
            integrated_sum += (100*impactMCsamples_sorted[MC])/sumwtotal
            if MC in Xsec.keys():
                if Xsec[MC]['unc'] == None:
                    RelXsec_unc = None
                    #print(f'!!!!! No unc found for {MC} !!!!! ')
                else:
                    RelXsec_unc = Xsec[MC]['unc']/Xsec[MC]['crossSec']
                    #print(f'rel unc: {RelXsec_unc}')
            else:
                RelXsec_unc = None
                #print(f'!!!!! No key found for {MC} !!!!!')
            my_dict[channel][MC]['prop'] = impactMCsamples_sorted[MC]/sumwtotal
            my_dict[channel][MC]['xsec'] = Xsec[MC]['crossSec']
            my_dict[channel][MC]['relunc'] = RelXsec_unc
            if RelXsec_unc == None:
                print(f'{MC} --> prop: {round(100*impactMCsamples_sorted[MC]/sumwtotal,2)}%, rel unc: None  --- integrated: {round(integrated_sum,2)}%')
            else:
                print(f'{MC} --> prop: {round(100*impactMCsamples_sorted[MC]/sumwtotal,2)}%, rel unc: {round(RelXsec_unc*100,2)}%  --- integrated: {round(integrated_sum,2)}%')
            if integrated_sum > 99:
                break
        print('')

    with open(output_file, 'w') as outfile:
        yaml.dump(my_dict, outfile)

    return my_dict
