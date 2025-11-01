import ROOT
import os
import yaml
import numpy as np
from ROOT import TCanvas, TGraph, TLegend
import CMS_lumi, tdrstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)

tag = 'FinalProd'
period = 'All' # 'All' if all year combined, else '2018', ...

chTolatex = {
    'ttm':'#tau #tau #mu',
    'tem':'#tau e #mu',
    'tee':'#tau e e',
    'tmm':'#tau #mu #mu',
    'tte':'#tau #tau e'
}

# PLOT upper limits
def PlotExpLimits(V2_lim, outputfile, period, channel):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
 
    # CMS style
    CMS_lumi.cmsText = "CMS Preliminary"
    #CMS_lumi.extraText = "Preliminary"
    CMS_lumi.cmsTextSize = 0.65
    CMS_lumi.outOfFrame = True
    tdrstyle.setTDRStyle()

    HNL_mass_range = list(V2_lim.keys())
    N = len(HNL_mass_range)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    
    up2s = [ ]
    for i in range(N):
        limit = V2_lim[HNL_mass_range[i]]
        up2s.append(limit[4])
        yellow.SetPoint(    i,    HNL_mass_range[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    HNL_mass_range[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    HNL_mass_range[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, HNL_mass_range[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, HNL_mass_range[i], limit[0] ) # - 2 sigma
    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    #c.SetLogy()
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
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
    #frame.SetMaximum(max(up2s)*1.05)
    frame.SetMaximum(1)
    frame.GetXaxis().SetLimits(min(HNL_mass_range),max(HNL_mass_range))
 
    c.SetLogy()
    c.SetLogx()

    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')

    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')
 
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')
 
    if period in ['2018','2017','2016', '2016_HIPM']:
        CMS_lumi.CMS_lumi(c,13,11)
    else:
        CMS_lumi.CMS_lumi(c,15,11)

    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.6
    x2 = x1 + 0.24
    y1 = 0.30
    y2 = y1 + 0.16
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)

    if period in ['2018','2017','2016', '2016_HIPM']:
        legend.AddEntry(median, f"Run2 {period} expected for {chTolatex[channel]}",'L')
    else:
        legend.AddEntry(median, f"Run2 expected for {chTolatex[channel]}",'L')
    
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
    legend.Draw()

    print(" ")
    c.SaveAs(outputfile)
    c.Close()

def main():
    channels = ['tee','tmm','tte','ttm','tem']

    if period in ['2018','2017','2016', '2016_HIPM']:
        BDV_file_name = f'BDV_inputs_{period}.yaml'
    elif period == 'All':
        BDV_file_name = 'BDV_inputs.yaml'
    else:
        print('period parmeter unknown')

    with open(f"F_Exp_limit/results/{tag}/{BDV_file_name}", 'r') as yaml_file:
        BDV = yaml.load(yaml_file, Loader=yaml.FullLoader)

    outputdir = f'{os.getenv("RUN_PATH")}/F_Exp_limit/figures/{tag}/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    CombinePath = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}/'

    for channel in channels:
        if period in ['2018','2017','2016', '2016_HIPM']:
            fig_name =  f'ExpLimit_BDV_{channel}_{period}.pdf'
            input_file = f'intPoints_{channel}_{period}.yaml'
        else:
            fig_name =  f'ExpLimit_BDV_{channel}.pdf'
            input_file = f'intPoints_{channel}.yaml'

        outputfile = outputdir + fig_name
        V2_lim = {}
        for HNL_mass in list(BDV[channel].keys()):
            #print(f'HNL_mass: {HNL_mass}')
            file_path = os.path.join(CombinePath, input_file)
            with open(file_path, 'r') as yaml_file:
                intPoints = yaml.load(yaml_file, Loader=yaml.FullLoader)
            #print(f'BDV: {BDV[channel][HNL_mass]}')
            V2_lim[HNL_mass] = np.array(intPoints[BDV[channel][HNL_mass]][HNL_mass])
            if len(V2_lim[HNL_mass]) == 5:
                print(f'HNL_mass: {HNL_mass}')
                #print(f'V2_lim: {V2_lim[HNL_mass]}')
            else:
                print(f'HNL_mass: {HNL_mass}')
                print(f'BDV: {BDV[channel][HNL_mass]}')
                print(f'!!!!!! V2_lim: {V2_lim[HNL_mass]} !!!!!!')
                print(" ")
                V2_lim[HNL_mass] = np.array([V2_lim[HNL_mass], V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass]])
        PlotExpLimits(V2_lim, outputfile, period, channel)
 
if __name__ == '__main__':
    main()
