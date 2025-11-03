import ROOT
import os
import yaml
import numpy as np
from ROOT import TCanvas, TGraph, TLegend
import CMS_lumi, tdrstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)

tag = 'FinalProd'
period = '2018' # 'All' if all year combined, else '2018', ...

chTolatex = {
    'tee':'#tau e e',
    'tee_ss':'#tau e^{#pm} e^{#pm} (SS)',
    'tee_os':'#tau e^{#pm} e^{#mp} (OS)',
    'tmm':'#tau #mu #mu',
    'tmm_ss':'#tau #mu^{#pm} #mu^{#pm} (SS)',
    'tmm_os':'#tau #mu^{#pm} #mu^{#mp} (OS)'
}

# PLOT upper limits
def PlotExpLimits(V2_lim, outputfile, period, channel, V2_lim_ss, V2_lim_os, V2_lim_combined):
    """
    Plot expected |V_{Nτ}|^2 limits vs m_HNL with 68%/95% bands from V2_lim,
    plus median-only curves for SS/OS (dotted) and combined (solid).
    - V2_lim*, for each HNL mass, should map to a 5-element iterable:
      [ -2σ, -1σ, median, +1σ, +2σ ]
    """
    # --- CMS style -----------------------------------------------------------
    CMS_lumi.cmsText = "CMS Preliminary"
    CMS_lumi.cmsTextSize = 0.65
    CMS_lumi.outOfFrame = True
    tdrstyle.setTDRStyle()

    # --- Prep ranges & graphs -----------------------------------------------
    HNL_mass_range = list(V2_lim.keys())
    N = len(HNL_mass_range)

    yellow = TGraph(2*N)   # 95% band
    green  = TGraph(2*N)   # 68% band
    median = TGraph(N)     # inclusive median from V2_lim

    # Additional medians
    g_med_ss       = TGraph(N)  # SS (dotted)
    g_med_os       = TGraph(N)  # OS (dotted)
    g_med_combined = TGraph(N)  # combined (solid)

    up2s = []
    for i in range(N):
        m = HNL_mass_range[i]
        l_incl = V2_lim[m]
        l_ss   = V2_lim_ss[m]
        l_os   = V2_lim_os[m]
        l_comb = V2_lim_combined[m]

        up2s.append(l_incl[4])

        # bands from inclusive (V2_lim)
        yellow.SetPoint(    i,    m, l_incl[4])   # +2σ
        green.SetPoint(     i,    m, l_incl[3])   # +1σ
        median.SetPoint(    i,    m, l_incl[2])   # median
        green.SetPoint(  2*N-1-i, m, l_incl[1])   # -1σ
        yellow.SetPoint( 2*N-1-i, m, l_incl[0])   # -2σ

        # medians for SS / OS / combined
        g_med_ss.SetPoint(i, m, l_ss[2])
        g_med_os.SetPoint(i, m, l_os[2])
        g_med_combined.SetPoint(i, m, l_comb[2])

    # --- Canvas & frame ------------------------------------------------------
    W, H = 800, 600
    T, B, L, R = 0.08*H, 0.12*H, 0.12*W, 0.04*W
    c = TCanvas("c", "c", 100, 100, W, H)
    c.SetFillColor(0); c.SetBorderMode(0)
    c.SetFrameFillStyle(0); c.SetFrameBorderMode(0)
    c.SetLeftMargin(L/W); c.SetRightMargin(R/W)
    c.SetTopMargin(T/H);  c.SetBottomMargin(B/H)
    c.SetTickx(0); c.SetTicky(0)
    c.SetGrid()
    c.cd()

    frame = c.DrawFrame(1.4, 0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.0)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("|V_{N#tau}|^{2}")
    frame.GetXaxis().SetTitle("m_{HNL} [GeV]")
    frame.SetMinimum(1e-3)
    frame.SetMaximum(1)
    frame.GetXaxis().SetLimits(min(HNL_mass_range), max(HNL_mass_range))

    c.SetLogy(); c.SetLogx()

    # --- Draw bands ----------------------------------------------------------
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')

    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')

    # --- Draw medians --------------------------------------------------------
    # Inclusive (from V2_lim): dashed black (kept from original)
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(1)  # solid
    median.Draw('Lsame')

    # SS: dotted, blue
    g_med_ss.SetLineColor(ROOT.kBlue + 1)
    g_med_ss.SetLineWidth(2)
    g_med_ss.SetLineStyle(2)  # dotted
    g_med_ss.Draw('Lsame')

    # OS: dotted, red
    g_med_os.SetLineColor(ROOT.kRed + 1)
    g_med_os.SetLineWidth(2)
    g_med_os.SetLineStyle(2)  # dotted
    g_med_os.Draw('Lsame')

    # Combined: solid, black (thicker)
    g_med_combined.SetLineColor(1)
    g_med_combined.SetLineWidth(2)
    g_med_combined.SetLineStyle(2)  # dotted
    g_med_combined.Draw('Lsame')

    # --- Luminosity stamp ----------------------------------------------------
    if period in ['2018','2017','2016','2016_HIPM']:
        CMS_lumi.CMS_lumi(c, 13, 11)
    else:
        CMS_lumi.CMS_lumi(c, 15, 11)

    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    # --- Legend --------------------------------------------------------------
    x1, y1 = 0.60, 0.30
    x2, y2 = x1 + 0.28, y1 + 0.22
    legend = TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)

    # Period text (match original behavior)
    run_label = f"Run2 {period}" if period in ['2018','2017','2016','2016_HIPM'] else "Run2"

    # Labels (use dict safely in case a key is missing)
    ch_label_incl = chTolatex.get(channel, channel)
    ch_label_ss   = chTolatex.get(f"{channel}_ss", f"{ch_label_incl} (SS)")
    ch_label_os   = chTolatex.get(f"{channel}_os", f"{ch_label_incl} (OS)")

    legend.AddEntry(g_med_combined, f"{ch_label_incl} (SS+OS comb.)", 'L')
    legend.AddEntry(g_med_ss,       f"{ch_label_ss}", 'L')
    legend.AddEntry(g_med_os,       f"{ch_label_os}", 'L')
    legend.AddEntry(median,         f"{ch_label_incl} (baseline)", 'L')
    legend.AddEntry(green,          "#pm 1 std. deviation", 'f')
    legend.AddEntry(yellow,         "#pm 2 std. deviation", 'f')
    legend.Draw()

    # --- Output --------------------------------------------------------------
    print(" ")
    c.SaveAs(outputfile)
    c.Close()

    # # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
 
    # # CMS style
    # CMS_lumi.cmsText = "CMS Preliminary"
    # #CMS_lumi.extraText = "Preliminary"
    # CMS_lumi.cmsTextSize = 0.65
    # CMS_lumi.outOfFrame = True
    # tdrstyle.setTDRStyle()

    # HNL_mass_range = list(V2_lim.keys())
    # N = len(HNL_mass_range)
    # yellow = TGraph(2*N)    # yellow band
    # green = TGraph(2*N)     # green band
    # median = TGraph(N)      # median line
    
    # up2s = [ ]
    # for i in range(N):
    #     limit = V2_lim[HNL_mass_range[i]]
    #     up2s.append(limit[4])
    #     yellow.SetPoint(    i,    HNL_mass_range[i], limit[4] ) # + 2 sigma
    #     green.SetPoint(     i,    HNL_mass_range[i], limit[3] ) # + 1 sigma
    #     median.SetPoint(    i,    HNL_mass_range[i], limit[2] ) # median
    #     green.SetPoint(  2*N-1-i, HNL_mass_range[i], limit[1] ) # - 1 sigma
    #     yellow.SetPoint( 2*N-1-i, HNL_mass_range[i], limit[0] ) # - 2 sigma
    # W = 800
    # H  = 600
    # T = 0.08*H
    # B = 0.12*H
    # L = 0.12*W
    # R = 0.04*W
    # c = TCanvas("c","c",100,100,W,H)
    # #c.SetLogy()
    # c.SetFillColor(0)
    # c.SetBorderMode(0)
    # c.SetFrameFillStyle(0)
    # c.SetFrameBorderMode(0)
    # c.SetLeftMargin( L/W )
    # c.SetRightMargin( R/W )
    # c.SetTopMargin( T/H )
    # c.SetBottomMargin( B/H )
    # c.SetTickx(0)
    # c.SetTicky(0)
    # c.SetGrid()
    # c.cd()
    # frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    # frame.GetYaxis().CenterTitle()
    # frame.GetYaxis().SetTitleSize(0.05)
    # frame.GetXaxis().SetTitleSize(0.05)
    # frame.GetXaxis().SetLabelSize(0.04)
    # frame.GetYaxis().SetLabelSize(0.04)
    # frame.GetYaxis().SetTitleOffset(1.0)
    # frame.GetXaxis().SetNdivisions(508)
    # frame.GetYaxis().CenterTitle(True)
    # frame.GetYaxis().SetTitle("|V_{N#tau}|^{2}")
    # frame.GetXaxis().SetTitle("m_{HNL} [GeV]")
    # frame.SetMinimum(1e-3)
    # #frame.SetMaximum(max(up2s)*1.05)
    # frame.SetMaximum(1)
    # frame.GetXaxis().SetLimits(min(HNL_mass_range),max(HNL_mass_range))
 
    # c.SetLogy()
    # c.SetLogx()

    # yellow.SetFillColor(ROOT.kOrange)
    # yellow.SetLineColor(ROOT.kOrange)
    # yellow.SetFillStyle(1001)
    # yellow.Draw('F')

    # green.SetFillColor(ROOT.kGreen+1)
    # green.SetLineColor(ROOT.kGreen+1)
    # green.SetFillStyle(1001)
    # green.Draw('Fsame')
 
    # median.SetLineColor(1)
    # median.SetLineWidth(2)
    # median.SetLineStyle(2)
    # median.Draw('Lsame')
 
    # if period in ['2018','2017','2016', '2016_HIPM']:
    #     CMS_lumi.CMS_lumi(c,13,11)
    # else:
    #     CMS_lumi.CMS_lumi(c,15,11)

    # ROOT.gPad.SetTicks(1,1)
    # frame.Draw('sameaxis')

    # x1 = 0.6
    # x2 = x1 + 0.24
    # y1 = 0.30
    # y2 = y1 + 0.16
    # legend = TLegend(x1,y1,x2,y2)
    # legend.SetFillStyle(0)
    # legend.SetBorderSize(0)
    # legend.SetTextSize(0.041)
    # legend.SetTextFont(42)

    # if period in ['2018','2017','2016', '2016_HIPM']:
    #     legend.AddEntry(median, f"Run2 {period} expected for {chTolatex[channel]}",'L')
    # else:
    #     legend.AddEntry(median, f"Run2 expected for {chTolatex[channel]}",'L')
    
    # legend.AddEntry(green, "#pm 1 std. deviation",'f')
    # legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
    # legend.Draw()

    # print(" ")
    # c.SaveAs(outputfile)
    # c.Close()

def main():
    channels = ['tee','tmm']

    if period in ['2018','2017','2016', '2016_HIPM']:
        BDV_file_name = f'BDV_inputs_{period}.yaml'
    elif period == 'All':
        BDV_file_name = 'BDV_inputs.yaml'
    else:
        print('period parmeter unknown')

    with open(f"F_Exp_limit/results/{tag}_OGchannels/{BDV_file_name}", 'r') as yaml_file:
        BDV = yaml.load(yaml_file, Loader=yaml.FullLoader)

    with open(f"F_Exp_limit/results/{tag}/{BDV_file_name}", 'r') as yaml_file:
        BDV_ssos = yaml.load(yaml_file, Loader=yaml.FullLoader)

    outputdir = f'{os.getenv("RUN_PATH")}/F_Exp_limit/figures/{tag}/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    CombinePath = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}_OGchannels/'

    CombinePath_SS_OS = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}/'

    for channel in channels:
        if period in ['2018','2017','2016', '2016_HIPM']:
            fig_name =  f'tll_study_{channel}_{period}.pdf'
            input_file = f'intPoints_{channel}_{period}.yaml'
            input_file_ss = f'intPoints_{channel}_ss_{period}.yaml'
            input_file_os = f'intPoints_{channel}_os_{period}.yaml'
            input_file_combined = f'intPoints_all_channels_{channel}_{period}.yaml'
        else:
            fig_name =  f'tll_study_{channel}.pdf'
            input_file = f'intPoints_{channel}.yaml'
            input_file_ss = f'intPoints_{channel}_ss.yaml'
            input_file_os = f'intPoints_{channel}_os.yaml'
            input_file_combined = f'intPoints_all_channels_{channel}.yaml'

        outputfile = outputdir + fig_name
        V2_lim = {}
        V2_lim_ss = {}
        V2_lim_os = {}
        V2_lim_combined = {}

        for HNL_mass in list(BDV[channel].keys()):
            #print(f'HNL_mass: {HNL_mass}')
            file_path = os.path.join(CombinePath, input_file)
            with open(file_path, 'r') as yaml_file:
                intPoints = yaml.load(yaml_file, Loader=yaml.FullLoader)

            file_path_ss = os.path.join(CombinePath_SS_OS, input_file_ss)
            with open(file_path_ss, 'r') as yaml_file:
                intPoints_ss = yaml.load(yaml_file, Loader=yaml.FullLoader)

            file_path_os = os.path.join(CombinePath_SS_OS, input_file_os)
            with open(file_path_os, 'r') as yaml_file:
                intPoints_os = yaml.load(yaml_file, Loader=yaml.FullLoader)

            file_path_combined = os.path.join(CombinePath_SS_OS, input_file_combined)
            with open(file_path_combined, 'r') as yaml_file:
                intPoints_combined = yaml.load(yaml_file, Loader=yaml.FullLoader)

            #print(f'BDV: {BDV[channel][HNL_mass]}')
            V2_lim[HNL_mass] = np.array(intPoints[BDV[channel][HNL_mass]][HNL_mass])
            V2_lim_ss[HNL_mass] = np.array(intPoints_ss[BDV_ssos[f'{channel}_ss'][HNL_mass]][HNL_mass])
            V2_lim_os[HNL_mass] = np.array(intPoints_os[BDV_ssos[f'{channel}_os'][HNL_mass]][HNL_mass])
            V2_lim_combined[HNL_mass] = np.array(intPoints_combined[HNL_mass])
            # if len(V2_lim[HNL_mass]) == 5:
            #     print(f'HNL_mass: {HNL_mass}')
            #     print(f'V2_lim: {V2_lim[HNL_mass]}')
            # else:
            #     print(f'HNL_mass: {HNL_mass}')
            #     print(f'BDV: {BDV[channel][HNL_mass]}')
            #     print(f'!!!!!! V2_lim: {V2_lim[HNL_mass]} !!!!!!')
            #     print(" ")
            #     V2_lim[HNL_mass] = np.array([V2_lim[HNL_mass], V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass]])
        PlotExpLimits(V2_lim, outputfile, period, channel, V2_lim_ss, V2_lim_os, V2_lim_combined)
 
if __name__ == '__main__':
    main()
