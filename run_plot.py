from Plotter import Plotter, load
import os
import argparse

parser = argparse.ArgumentParser(description='Running code to produce plots.')
parser.add_argument('--channel', required=True, type=str, help="Channel to plot.")
parser.add_argument('--tag', required=True, type=str, help="Name of the file where anatuple are stored.")
args = parser.parse_args()

channel = args.channel
if channel not in ['ttm', 'tmm', 'tte', 'tee', 'tem']:
    raise "channel not valid"

input_folder = f'/eos/user/p/pdebryas/HNL/anatuple/nanoV10/{args.tag}/{channel}/anatuple'
#common config
page_cfg = 'config/shared/cms_stacked.yaml'
page_cfg_custom = 'config/shared/2018.yaml'
hist_maker_file = f'hist_makers/hnl_maker.py'
#channel config
hist_cfg = f'config/hnl/hnl_{channel}/histograms.yaml'
inputs_cfg = f'config/hnl/hnl_{channel}/inputs.yaml'
output_dir = f'output/{channel}/'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if channel == 'ttm':
    custom_title ="cat_text=#tau_{h}#tau_{h}#mu"
    var_to_plot = ['Tau1_pt', 'Tau2_pt', 'Muon_pt', 'pt_sum_ttm']

if channel == 'tmm':
    custom_title ="cat_text=#tau_{h}#mu#mu"
    var_to_plot = ['Tau_pt', 'Muon1_pt', 'Muon2_pt', 'pt_sum_tmm']

if channel == 'tte':
    custom_title ="cat_text=#tau_{h}#tau_{h} e"
    var_to_plot = ['Tau1_pt', 'Tau2_pt', 'Electron_pt', 'pt_sum_tte']

if channel == 'tee':
    custom_title ="cat_text=#tau_{h} e e"
    var_to_plot = ['Tau_pt', 'Electron1_pt', 'Electron2_pt', 'pt_sum_tee']

if channel == 'tem':
    custom_title ="cat_text=#tau_{h} e #mu"
    var_to_plot = ['Tau_pt', 'Electron_pt', 'Muon_pt', 'pt_sum_tem']

custom = None if custom_title is None else dict(item.split('=') for item in custom_title.split(','))
plotter = Plotter(page_cfg=page_cfg, page_cfg_custom=page_cfg_custom, hist_cfg=hist_cfg, inputs_cfg=inputs_cfg)
hist_maker = load(hist_maker_file)

for var in var_to_plot:
    print("---------------------- Plotting "+var + "----------------------")
    hists = hist_maker.make_histograms(input_folder, hist_name=var, hist_cfg=plotter.hist_cfg, inputs_cfg=plotter.inputs_cfg, channel = channel)
    plotter.plot(var, hists, output_dir+f'{var}.pdf', custom=custom)
