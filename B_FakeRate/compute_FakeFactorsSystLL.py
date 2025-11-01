import os
import yaml
import numpy as np
import correctionlib.convert
import hist
import ROOT
import math

from common.helpers import load_ntuples
from common.regions.regions import compute_region_mask
from D_RootHist.hist_makers.helpers import ToRootHist_val
from B_FakeRate.helpers import compute_ptcorr, apply_FR_methodLL, ratio_func, AdaptBin

#computing systematic uncertainty on Fake Factors

#parameters -----------------------------------------------------------------------------------------------------
#period = period of data taking
period = '2018'
#tag = tag used for LL FR computation
tag  = 'FinalProd'
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

        Results[Lepton][RegionName]['ratio'] = abs(ratio_func(hists['ExpBackground'], hists['MissIDPred'], hists['PromptLeptons'], hists['data'])).tolist()
        Results[Lepton][RegionName]['bin']   = pt_bins.tolist()

        # --- convert coffea hist -> ROOT TH1 and enable Sumw2 ---
        h_data   = ToRootHist_val(hists['data']);       h_data.Sumw2()
        h_missid = ToRootHist_val(hists['MissIDPred']); h_missid.Sumw2()
        h_prompt = ToRootHist_val(hists['PromptLeptons']); h_prompt.Sumw2()

        # --- style ---
        h_missid.SetLineColor(ROOT.kBlue-10)
        h_missid.SetFillColor(ROOT.kBlue-10)

        h_prompt.SetLineColor(ROOT.kGreen+2)
        h_prompt.SetFillColor(ROOT.kGreen+2)

        h_data.SetLineColor(ROOT.kBlack)
        h_data.SetMarkerStyle(20)

        # --- build stack and total MC ---
        stack = ROOT.THStack("stack","")
        stack.Add(h_missid)
        stack.Add(h_prompt)

        h_mc_total = h_missid.Clone("h_mc_total")
        h_mc_total.Add(h_prompt)

        # --- combined MC stat uncertainty (top pad band) ---
        h_unc = h_mc_total.Clone("h_unc")
        n_bins = h_unc.GetNbinsX()
        for ib in range(1, n_bins+1):
            err = math.sqrt(h_missid.GetBinError(ib)**2 + h_prompt.GetBinError(ib)**2)
            h_unc.SetBinError(ib, err)
        h_unc.SetFillColor(ROOT.kGray+2)
        h_unc.SetFillStyle(3354)   # hatched
        h_unc.SetLineWidth(0)

        # --- axis ranges and titles on underlying histograms ---
        h_mc_total.SetMinimum(0)
        h_data.SetMinimum(0)
        max_val = max(h_mc_total.GetMaximum(), h_data.GetMaximum())
        h_mc_total.SetMaximum(20 + max_val)
        h_data.SetMaximum(20 + max_val)

        # --- canvas and ratio plot (canvas must be created before TRatioPlot) ---
        ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "A ratio example")
        c1.SetTicks(0, 1)
        c1.SetLogx()

        rp = ROOT.TRatioPlot(h_data, h_mc_total)
        rp.Draw()

        # configure lower (ratio) axis labels / range
        lower_graph = rp.GetLowerRefGraph()
        lower_graph.GetXaxis().SetTitle("p_{T}^{corr} [GeV]")
        lower_graph.GetYaxis().SetTitle("ratio")
        lower_graph.SetMinimum(0.5)
        lower_graph.SetMaximum(1.5)

        # --- upper pad: draw stack, MC uncertainty band, data ---
        upper_pad = rp.GetUpperPad()
        upper_pad.cd()
        stack.Draw("hist same")
        h_unc.Draw("E2 same")   # hashed uncertainty band on top of stack
        h_data.Draw("E same")
        #upper_pad.GetYaxis().SetTitle("Events / bin")
        upper_pad.Update()

        # --- lower pad: draw relative uncertainty band around 1 ---
        lower_pad = rp.GetLowerPad()
        lower_pad.cd()

        h_rel_unc = h_mc_total.Clone("h_rel_unc")
        for ib in range(1, h_rel_unc.GetNbinsX()+1):
            mc = h_mc_total.GetBinContent(ib)
            err = math.sqrt(h_missid.GetBinError(ib)**2 + h_prompt.GetBinError(ib)**2)
            if mc > 0:
                h_rel_unc.SetBinContent(ib, 1.0)
                h_rel_unc.SetBinError(ib, err / mc)
            else:
                h_rel_unc.SetBinContent(ib, 1.0)
                h_rel_unc.SetBinError(ib, 0.0)

        h_rel_unc.SetFillColor(ROOT.kGray+2)
        h_rel_unc.SetFillStyle(3354)
        h_rel_unc.SetLineWidth(0)
        h_rel_unc.Draw("E2 same")

        # redraw ratio graph on top
        lower_graph.Draw("same p")
        lower_graph.GetXaxis().SetTitle("p_{T}^{corr} [GeV]")
        lower_pad.Update()

        # --- legend + save ---
        upper_pad.cd()
        legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.85)
        legend.SetBorderSize(0)
        legend.AddEntry(h_data, "Data", "lep")
        legend.AddEntry(h_prompt, "Prompt leptons", "f")
        legend.AddEntry(h_missid, "Miss id. prediction", "f")
        legend.AddEntry(h_unc, "MC stat. unc.", "f")
        legend.Draw()

        c1.SaveAs(output_figures + f"Syst_Unc_{Lepton}FFs_{RegionUnc}.pdf")


# Save dictionary to a YAML file
with open(output_results + 'Syst_Unc_LLFF.yaml', 'w') as file:
    yaml.dump(Results, file, default_flow_style=False)
