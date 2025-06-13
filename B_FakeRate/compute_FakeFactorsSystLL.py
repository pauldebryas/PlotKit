import os
import yaml
import numpy as np
import correctionlib.convert
import hist
import ROOT

from common.helpers import load_ntuples
from common.regions.regions import compute_region_mask
from D_RootHist.hist_makers.helpers import ToRootHist_val
from B_FakeRate.helpers import compute_ptcorr, apply_FR_methodLL, ratio_func, AdaptBin

#computing systematic uncertainty on Fake Factors

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for LL FR computation
tag  = 'LightLeptFFV2'
#bins = either int and bin size such that statistic is the same for all bins or 'config' and pt and eta bins from config file 
bins = 'config'
#----------------------------------------------------------------------------------------------------------------
Leptons = ['Electron','Muon']

output_results = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/')
if not os.path.exists(output_results):
    os.makedirs(output_results) 

# Load FR
FR = {}
FR['Electron_ttbar'] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactorsElectron/ttbarRegionFR/P_fake_TL_Electron.json'))
FR['Electron_DY'] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactorsElectron/DYRegionFR/P_fake_TL_Electron.json'))
FR['Muon_ttbar'] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactorsMuon/ttbarRegionFR/P_fake_TL_Muon.json'))
FR['Muon_DY'] =  correctionlib.CorrectionSet.from_file(os.path.join(os.getenv("RUN_PATH"), f'B_FakeRate/results/{tag}/{period}/FakeFactorsMuon/DYRegionFR/P_fake_TL_Muon.json'))

# load Corr factor from YAML file
with open(os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/results/{tag}/{period}/correctionFactorsLL.yml'), 'r') as file:
    CorrFactor = yaml.safe_load(file)

Results = {}
for Lepton in Leptons:
    Results[Lepton] = {}
    for RegionUnc in ['ttbarRegionValidation', 'DYRegionValidation']:
        RegionName = RegionUnc.replace('RegionValidation', '')
        Results[Lepton][RegionName] = {}
        print(f'Computing syst error on Fake Factors for {Lepton} evaluated in {RegionUnc}:')

        W_DY_ttbar = {}
        W_DY_ttbar['DY'] = [0, 0]
        W_DY_ttbar['ttbar']= [0, 0]
        if Lepton == 'Electron':
            output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactors{Lepton}/')
            if RegionUnc == 'ttbarRegionValidation':
                channel = 'lle'
                W_DY_ttbar['ttbar'][0] = 1.
            if RegionUnc == 'DYRegionValidation':
                channel = 'Ze'
                W_DY_ttbar['DY'][0] = 1.
        if Lepton == 'Muon':
            output_figures = os.path.join(os.getenv("RUN_PATH"),f'B_FakeRate/figures/{tag}/{period}/FakeFactors{Lepton}/')
            if RegionUnc == 'ttbarRegionValidation':
                channel = 'llmu'
                W_DY_ttbar['ttbar'][0] = 1.
            if RegionUnc == 'DYRegionValidation':
                channel = 'Zmu'
                W_DY_ttbar['DY'][0] = 1.

        if channel == None:
            print('Unidentified Lepton name: can be Electron or Muon')
            raise

        if not os.path.exists(output_figures):
            os.makedirs(output_figures)  

        if bins == 'config':
            nbins = []
            with open(f'{os.getenv("RUN_PATH")}/common/config/all/bin_FR.yaml', 'r') as f:
                binsconfig = yaml.safe_load(f)
            nbins.append(binsconfig[Lepton]['pt'])
            nbins.append(binsconfig[Lepton]['abseta'])
        else:
            nbins = bins
        
        # load input files
        inputs = {}
        print(f'Load inputs for channel {channel}')
        inputs_cfg_file = f'{os.getenv("RUN_PATH")}/common/config/all/inputs/inputs_TrueLepton.yaml'
        input_dir = os.path.join('/eos/user/p/pdebryas/HNL_LLFF/anatuple', period, tag, channel, 'anatuple')
        if Lepton == 'Electron':
            inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_tte/inputs_data.yaml' 
        elif Lepton == 'Muon':
            inputs_cfg_file_data = f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_ttm/inputs_data.yaml' 

        with open(inputs_cfg_file, 'r') as f:
            inputs_cfg = yaml.safe_load(f)

        with open(inputs_cfg_file_data, 'r') as f:
            inputs_cfg_data = yaml.safe_load(f)

        input_files = {}
        exclude_list = [] #['WJetsToLNu_HT-100To200_anatuple.root', 'WJetsToLNu_HT-1200To2500_anatuple.root']
        exclude_list_files = [os.path.join(input_dir, elem) for elem in exclude_list] 

        for input in inputs_cfg:
            if 'files' not in input.keys():
                raise f'Missing files dict for {input}'
            files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
            files_list_new = list(files_list)
            for file in files_list:
                if (os.path.isfile(file) == False) or (file in exclude_list_files):
                    print('WARNING: ' + file + ' is missing')
                    files_list_new.remove(file)
            input_files[input['name']] = files_list_new

        for input in inputs_cfg_data:
            if input['name'] == 'data':
                if 'files' not in input.keys():
                    raise f'Missing files dict for {input}'
                files_list = [os.path.join(input_dir, elem) for elem in input['files']] 
                files_list_new = list(files_list)
                for file in files_list:
                    if os.path.isfile(file) == False:
                        print('WARNING: ' + file + ' is missing')
                        files_list_new.remove(file)
                input_files[input['name']] = files_list_new

        inputs[channel] = load_ntuples(input_files)

        TrueLepton_pt = []
        TrueLepton_weight = []
        FakeLepton_pt = []
        FakeLepton_weight = []
        Data_pt = []
        Data_weights = []

        cut_region = compute_region_mask(inputs[channel]['TrueLepton'], channel, 'MC', RegionUnc)
        pt_corr =  compute_ptcorr(inputs[channel]['TrueLepton'], Lepton, CorrFactor)
        TrueLepton_pt = np.concatenate( (TrueLepton_pt, pt_corr[cut_region[f'{RegionUnc}_PassTightWP_{Lepton}IsPromptLepton']]) )
        TrueLepton_weight = np.concatenate( (TrueLepton_weight, np.array(inputs[channel]['TrueLepton']['genWeight']).flatten()[cut_region[f'{RegionUnc}_PassTightWP_{Lepton}IsPromptLepton']]) )

        cut_region = compute_region_mask(inputs[channel]['data'], channel, 'data', RegionUnc)
        pt_corr =  compute_ptcorr(inputs[channel]['data'], Lepton, CorrFactor)
        Data_pt = np.concatenate( (Data_pt, pt_corr[cut_region[f'{RegionUnc}_PassTightWP']]) )
        Data_weights = np.concatenate( (Data_weights, np.array(inputs[channel]['data']['genWeight']).flatten()[cut_region[f'{RegionUnc}_PassTightWP']]) )

        fakes_pt, weights = apply_FR_methodLL(FR, W_DY_ttbar, inputs[channel], channel, Lepton, RegionUnc, CorrFactor)
        FakeLepton_pt = np.concatenate( (FakeLepton_pt, fakes_pt) )
        FakeLepton_weight = np.concatenate( (FakeLepton_weight, weights) )

        if isinstance(nbins, int):
            pt_bins = AdaptBin(FakeLepton_pt, nbins).tolist()
            pt_bins = [round(i) for i in pt_bins]
            pt_bins[-1] += 1 
            pt_bins[0] -= 1 
        elif len(nbins) == 2:
            pt_bins = np.array(nbins[0])
            pt_bins[0]  = 1
            
        hists = {}
        hists['PromptLeptons'] = hist.Hist.new.Variable(pt_bins, name='x', overflow=True).Weight()
        hists['PromptLeptons'].fill(x=TrueLepton_pt, weight=TrueLepton_weight)
        hists['MissIDPred'] = hist.Hist.new.Variable(pt_bins, name='x', overflow=True).Weight()
        hists['MissIDPred'].fill(x=FakeLepton_pt, weight=FakeLepton_weight)
        hists['data'] = hist.Hist.new.Variable(pt_bins, name='x', overflow=True).Weight()
        hists['data'].fill(x=Data_pt, weight=Data_weights)
        hists['ExpBackground'] = hist.Hist.new.Variable(pt_bins, name='x', overflow=True).Weight()
        hists['ExpBackground'].fill(x=np.concatenate( (Data_pt,TrueLepton_pt)), weight=np.concatenate( (Data_weights,(-1)*TrueLepton_weight)))

        # Combine hist1 and hist2
        #stacked_hist = hists['TrueLepton'] + hists['FakeBackground']

        #print(f'binning: {pt_bins}')
        #print(f'PromptLeptons: {hists["PromptLeptons"].values()}')
        #print(f'FakeBackground: {hists["FakeBackground"].values()}')
        #print('')

        Results[Lepton][RegionName]['ratio'] = abs(ratio_func(hists['ExpBackground'], hists['MissIDPred'])).tolist()
        Results[Lepton][RegionName]['bin']   = pt_bins.tolist()

        # Setup ROOT style
        ROOT.gStyle.SetOptStat(0)

        # Create canvas
        c1 = ROOT.TCanvas("c1", "A ratio example")

        # Convert histograms to ROOT TH1 objects
        h1 = ToRootHist_val(hists['MissIDPred'])
        h2 = ToRootHist_val(hists['ExpBackground'])

        # Set line colors and fill colors for stack components
        h1.SetLineColor(ROOT.kBlue-10)
        h1.SetFillColor(ROOT.kBlue-10)
        h2.SetLineColor(ROOT.kBlack)
        h2.SetMarkerStyle(20)

        # Create the ratio plot with the total expected background and the data
        rp = ROOT.TRatioPlot(h2, h1)

        # Customize the canvas
        c1.SetTicks(0, 1)
        c1.SetLogx()
        
        max_val = max([max(hists['ExpBackground'].values()), max(hists['MissIDPred'].values())])

        h1.SetMinimum(0)
        h1.SetMaximum(20 + max_val)

        h2.SetMinimum(0)
        h2.SetMaximum(20 + max_val)
        h2.GetYaxis().SetTitle('Events / bin')

        rp.GetXaxis().SetTitle(f'{Lepton} pt [GeV]')

        # Draw the ratio plot
        rp.Draw()

        # Get the upper pad to draw the stack histograms and the data histogram
        upper_pad = rp.GetUpperPad()
        upper_pad.cd()

        # Draw the stacked histograms
        h1.Draw("hist same")

        # Draw the data histogram on the same canvas
        h2.Draw("same E")

        # Update the upper pad to refresh the axis
        upper_pad.Update()

        # Add a legend
        legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)
        legend.AddEntry(h1, "Miss id. rate prediction", "f")
        legend.AddEntry(h2, "Data", "lep")
        legend.Draw()

        # Save the plot
        c1.SaveAs(output_figures + f"Syst_Unc_{Lepton}FFs_{RegionUnc}.pdf")

        '''
        # Create a stack for the background components
        hs = ROOT.THStack("hs", "")

        # Add components to the stack
        hs.Add(h1_fake)
        hs.Add(h1_true)

        # Create a histogram for the total expected background by adding the components
        h1 = h1_fake.Clone("h1")
        h1.Add(h1_true)

        # Create the ratio plot with the total expected background and the data
        rp = ROOT.TRatioPlot(h2, h1)

        # Customize the canvas
        c1.SetTicks(0, 1)
        c1.SetLogx()
        h1.GetXaxis().SetTitle(f'{Lepton} pt corr [GeV]')
        h1.GetYaxis().SetTitle('Events / bin')
        h1.SetMinimum(0)
        h1.SetMaximum(20 + max([max(hists['ExpBackground'].values()), max(hists['data'].values())]))

        # Draw the ratio plot
        rp.Draw()

        # Get the upper pad to draw the stack histograms and the data histogram
        upper_pad = rp.GetUpperPad()
        upper_pad.cd()

        # Set the same y-axis range for the stack
        hs.SetMinimum(0)
        hs.SetMaximum(20 + max([max(hists['ExpBackground'].values()), max(hists['data'].values())]))

        # Draw the stacked histograms
        hs.Draw("hist same")

        # Set the same y-axis range for the stack
        h2.SetMinimum(0)
        h2.SetMaximum(20 + max([max(hists['ExpBackground'].values()), max(hists['data'].values())]))

        # Draw the data histogram on the same canvas
        h2.Draw("same E")

        # Update the upper pad to refresh the axis
        upper_pad.Update()

        # Add a legend
        legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)
        legend.AddEntry(h1_fake, "Fake Background", "f")
        legend.AddEntry(h1_true, "True Lepton (MC)", "f")
        legend.AddEntry(h2, "Data", "lep")
        legend.Draw()

        # Save the plot
        c1.SaveAs(output_figures + f"Syst_Unc_{Lepton}FFs_{RegionUnc}.pdf")
        '''

# Save dictionary to a YAML file
with open(output_results + 'Syst_Unc_LLFF.yaml', 'w') as file:
    yaml.dump(Results, file, default_flow_style=False)


