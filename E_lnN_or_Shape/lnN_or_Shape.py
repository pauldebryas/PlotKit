import ROOT
import os
import yaml

from common.helpers import get_hnl_masses 

# can be took from ssh://cern/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/config.ini
tag = 'AddJETcorr'
period = '2016_HIPM'
channels = ['tee','tmm','tem','tte','ttm']
plot_fig = False

def var_from_channel(channel):
    config_file_channel = f'{os.getenv("RUN_PATH")}/common/config/all/histograms/histograms_{channel}.yaml'
    with open(config_file_channel, 'r') as f:
        config_channel = yaml.safe_load(f)
    vars_raw =  list(config_channel.keys())

    HNLMassRange = get_hnl_masses(period) 
    vars = []
    for var in vars_raw:
        for HNLMass in HNLMassRange:
            var_name = f'{var}_HNLMass{str(HNLMass)}'
            vars.append(var_name)
    return vars

def plot_fig_ratio(h1, h2, corr, ud, channel, var, fit_param = None):
    print(f'... plotting corr for {corr} {ud}')

    # Create a canvas for the ratio plot
    c1 = ROOT.TCanvas("c1", "Title", 800, 600)

    # Set different colors for h1 and h2
    h2.SetLineColor(ROOT.kBlue)
    h1.SetLineColor(ROOT.kBlack)

    h2_copy = h2.Clone("hist_ratio_up")
    for i in range(1, h1.GetNbinsX() + 1):  # Loop over bins (assuming 1D histogram)
        h2_copy.SetBinError(i, 0.00001)

    # Create a TRatioPlot object
    rp = ROOT.TRatioPlot(h2_copy, h1)

    # Customize the canvas and ratio plot
    c1.SetTicks(0, 1)

    # Draw the canvas and ratio plot
    c1.Draw()
    rp.Draw()

    # Add label for the x-axis using LaTeX
    config_file_channel = f'{os.getenv("RUN_PATH")}/common/config/all/histograms/histograms_{channel}.yaml'
    with open(config_file_channel, 'r') as f:
        config_channel = yaml.safe_load(f)
    config = config_channel
    var_raw = var.split('_HNLMass')[0]
    HNLMass = var.split('_HNLMass')[1]
    x_title = config[var_raw]['x_title']
    if var_raw == 'DNNscore':
        x_title = x_title + f'{HNLMass} GeV'
    h1.GetXaxis().SetTitle(x_title)  # Set x-axis label to eta (eta) using LaTeX

    if fit_param != None:
       # Fit a constant function to the ratio histogram
       fit_func = ROOT.TF1("fit_func", str(fit_param), h1.GetXaxis().GetXmin(), h1.GetXaxis().GetXmax())
       fit_func.SetLineColor(ROOT.kRed)
       rp.GetLowerPad().cd()  # Switch to the lower pad
       fit_func.Draw("same")

    # Add legend
    rp.GetUpperPad().cd()
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # Define legend position (x1, y1, x2, y2)
    legend.AddEntry(h1, process, "l")
    legend.AddEntry(h2, f"{corr} {ud}", "l")
    if fit_param != None:
        legend.AddEntry(fit_func, f"Fit Parameter: {fit_param:.{4}f}", "l")  # Add fit parameter to legend
    legend.Draw()

    # Save the ratio plot canvas as a PNG image
    save_dir = f'{os.getenv("RUN_PATH")}/E_lnN_or_Shape/figures/{period}/{tag}/{channel}/{var}/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    c1.SaveAs(save_dir + f"{process}_{corr}_{ud}.pdf")
    return

def constant_function(x,par):
    return par[0]

def GetChi2(histogram):
    fit_func = ROOT.TF1("fit_func", constant_function, 0, 10, 1)
    fit_func.SetParameter(0, 1.0)
    histogram.Fit(fit_func, "qn")
    chi2 = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()
    p_value = ROOT.TMath.Prob(chi2, ndf)
    fit_param = fit_func.GetParameter(0)
    fit_param_error = fit_func.GetParError(0)
    return chi2,p_value,fit_param,fit_param_error

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    channels = ['tee','tmm','tte','ttm','tem']
    for channel in channels:
        print(f'--------------------------------------------------- {channel} ---------------------------------------------------')
        vars = var_from_channel(channel)
        for var in vars:
            print(f'For var {var}:')
            # Open the ROOT file
            file_path = f'{os.getenv("RUN_PATH")}/D_RootHist/results/{period}/{tag}/{channel}/FakeRate/SignalRegion/TH1_{var}_HNL.root'
            file = ROOT.TFile.Open(file_path)

            # Check if file is opened successfully
            if not file or file.IsZombie():
                print("Error: Unable to open file")
                exit(1)

            # Get list of keys in the "signal_region/" folder
            keys = file.Get("signal_region").GetListOfKeys()

            # Loop through keys to find corrections
            corr_list = {}
            for key in keys:
                name = key.GetName()
                split_name = name.split('_')
                process = split_name[0]
                if process not in corr_list.keys():
                    corr_list[process] = []

                if (split_name[-1].endswith('Up')) or (split_name[-1].endswith('Down')):
                    if split_name[-1].endswith('Up'):
                        name = name.replace(process+'_', "")[:-2]
                    if split_name[-1].endswith('Down'):
                        name = name.replace(process+'_', "")[:-4]
                    if name not in corr_list[process]:
                        corr_list[process].append(name)

            # Check for each corrections if it's lnN or Shape unc.
            lnN_or_Shape = {}
            for process in corr_list.keys():
                if process in ['data']:
                    pass
                else:
                    print(f'- Process: {process}')
                    lnN_or_Shape[process] = {}
                    for corr in corr_list[process]:
                        lnN_or_Shape[process][corr] = {}
                        #print(f'... checking lnN or shape for {process}_{corr} ...')

                        # Get histograms h1 and h2 from the file and check if valid
                        h1 = file.Get(f"signal_region/{process}")
                        h2_Up = file.Get(f"signal_region/{process}_{corr}Up")
                        h2_Down = file.Get(f"signal_region/{process}_{corr}Down")
                        if not(h1 and h2_Up and h2_Down):
                            print(f"Error: Unable to retrieve histograms nom/Up/Down from the file.")
                            exit()  # Exit the script if histograms are not retrieved successfully

                        # compute ratio plots
                        h2_ratio_Up = h2_Up.Clone("hist_ratio_up")
                        h2_ratio_Down = h2_Down.Clone("hist_ratio_Down")
                        h2_ratio_Up.Divide(h1)
                        h2_ratio_Down.Divide(h1)

                        # Loop over bins to compute error (no classic prop error. since events are basically the same)
                        for i in range(1, h1.GetNbinsX() + 1):  # Loop over bins (assuming 1D histogram)
                            bin_value_h1 = h1.GetBinContent(i)  # Get bin content
                            bin_error_h2_Up = h2_Up.GetBinError(i)  # Get bin error
                            bin_error_h2_Down = h2_Down.GetBinError(i)  # Get bin error
                            if bin_value_h1 != 0:
                                h2_ratio_Up.SetBinError(i, bin_error_h2_Up/bin_value_h1) 
                                h2_ratio_Down.SetBinError(i, bin_error_h2_Down/bin_value_h1)
                            else:
                                h2_ratio_Up.SetBinError(i, 0) 
                                h2_ratio_Down.SetBinError(i, 0)

                        chi2_Up, p_value_Up, fit_param_Up, fit_param_error_Up = GetChi2(h2_ratio_Up)
                        chi2_Down, p_value_Down, fit_param_Down, fit_param_error_Down = GetChi2(h2_ratio_Down)

                        if (p_value_Up == 1.0) & (p_value_Down == 1.0):
                            lnN_or_Shape[process][corr]['type'] = 'lnN'
                            if min(chi2_Up, chi2_Down) == chi2_Up:
                                lnN_or_Shape[process][corr]['fit_param'] = fit_param_Down
                                lnN_or_Shape[process][corr]['fit_param_err'] = fit_param_error_Down
                                lnN_or_Shape[process][corr]['chi2'] = chi2_Down
                                lnN_or_Shape[process][corr]['p_value'] = p_value_Down
                            if min(chi2_Up, chi2_Down) == chi2_Down:
                                lnN_or_Shape[process][corr]['fit_param'] = fit_param_Up
                                lnN_or_Shape[process][corr]['fit_param_err'] = fit_param_error_Up
                                lnN_or_Shape[process][corr]['chi2'] = chi2_Up
                                lnN_or_Shape[process][corr]['p_value'] = p_value_Up
                        else:
                            lnN_or_Shape[process][corr]['type'] = 'Shape'
                            if min(chi2_Up, chi2_Down) == chi2_Up:
                                lnN_or_Shape[process][corr]['fit_param'] = fit_param_Down
                                lnN_or_Shape[process][corr]['fit_param_err'] = fit_param_error_Down
                                lnN_or_Shape[process][corr]['chi2'] = chi2_Down
                                lnN_or_Shape[process][corr]['p_value'] = p_value_Down
                            if min(chi2_Up, chi2_Down) == chi2_Down:
                                lnN_or_Shape[process][corr]['fit_param'] = fit_param_Up
                                lnN_or_Shape[process][corr]['fit_param_err'] = fit_param_error_Up
                                lnN_or_Shape[process][corr]['chi2'] = chi2_Up
                                lnN_or_Shape[process][corr]['p_value'] = p_value_Up
                        
                        if plot_fig:
                            if (not process.startswith('HNL')) & (lnN_or_Shape[process][corr]['p_value'] < 0.97):
                                if lnN_or_Shape[process][corr]['type'] == 'lnN':
                                    plot_fig_ratio(h1, h2_Up, corr, 'Up', channel, var, fit_param_Up)
                                    plot_fig_ratio(h1, h2_Down, corr, 'Down', channel, var, fit_param_Down)
                                else:
                                    plot_fig_ratio(h1, h2_Up, corr, 'Up', channel, var)
                                    plot_fig_ratio(h1, h2_Down, corr, 'Down', channel, var)                        
                            
            # Write the dictionary to the YAML file
            
            save_dir = f'{os.getenv("RUN_PATH")}/E_lnN_or_Shape/results/{period}/{tag}/{channel}/'
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

            file_path = save_dir + f'{var}.yaml'
            with open(file_path, 'w') as yaml_file:
                yaml.dump(lnN_or_Shape, yaml_file, default_flow_style=False)
            print('')
        
        print('')
        print('')
