import ROOT
import os
import yaml
import numpy as np
from ROOT import TCanvas, TGraph, TLegend
import CMS_lumi, tdrstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)

tag = 'AddJETcorr'

# PLOT upper limits
def PlotExpLimits(V2_lim, outputfile):
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
        print(HNL_mass_range[i])
        limit = V2_lim[HNL_mass_range[i]]
        if len(limit) ==1:
           print(limit)
           new_lim = []
           new_lim.append(limit[0])
           new_lim.append(limit[0])
           new_lim.append(limit[0])
           new_lim.append(limit[0])
           new_lim.append(limit[0])
           limit = new_lim
        #print(limit)
        up2s.append(limit[4])
        yellow.SetPoint(    i,    HNL_mass_range[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    HNL_mass_range[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    HNL_mass_range[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, HNL_mass_range[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, HNL_mass_range[i], limit[0] ) # - 2 sigma

    mass = [20.16456139061248, 26.8352346902974,36.29431136096506, 49.61356550140213,58.614110947900635, 74.59242495694417, 83.86810310829574, 98.11748111433404,137.82937927094517, 175.65989929607008, 233.7480085567492, 338.99815299846193,491.5926903904981, 720.589843462258, 928.0560570265741]
    lim =  [0.0005019100441728138,0.0005122034362867965,0.0005420013832510946,0.0007123808769013648,0.0015507257967259716,0.023750341304095056,0.27220281774884153,0.10849934664589993,0.08759253331025571,0.10899607477957826,0.1308625474466631,0.22973062550925316,0.4744717187020259,1.0160725186361823,2.096582937391052]
    EXO_22_011 = TGraph(len(mass))
    for j in range(len(mass)):
        EXO_22_011.SetPoint(j,mass[j],lim[j])

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

    EXO_22_011.SetLineColor(ROOT.kRed)
    EXO_22_011.SetLineWidth(2)
    EXO_22_011.SetLineStyle(2)
    EXO_22_011.Draw('Lsame')
 
    CMS_lumi.CMS_lumi(c,15,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.45
    x2 = x1 + 0.24
    y1 = 0.30
    y2 = y1 + 0.16
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(median, "Run2 expected",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
    legend.AddEntry(EXO_22_011, "EXO_22_011 Expected limit",'L')
    legend.Draw()

    print(" ")
    c.SaveAs(outputfile)
    c.Close()

def main():

    outputdir = f'{os.getenv("RUN_PATH")}/F_Exp_limit/figures/{tag}/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    CombinePath = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}/'

    fig_name =  f'ExpLimit_BDV_all_channels.pdf'
    #fig_name =  f'ExpLimit_BDV_tem_tmm_tee.pdf'
    input_file_name = f'intPoints_all_channels.yaml'
    #input_file_name = 'intPoints_tem_tmm_tee.yaml'

    outputfile = outputdir + fig_name
    input_file_path = os.path.join(CombinePath, input_file_name)

    with open(input_file_path, 'r') as yaml_file:
        intPoints = yaml.load(yaml_file, Loader=yaml.FullLoader)

    V2_lim = {}
    for HNL_mass in list(intPoints.keys()):
        if HNL_mass in [30, 1000]:
            continue
        print(f'HNL_mass: {HNL_mass}')
        V2_lim[HNL_mass] = np.array(intPoints[HNL_mass])
    PlotExpLimits(V2_lim, outputfile)
 
if __name__ == '__main__':
    main()
