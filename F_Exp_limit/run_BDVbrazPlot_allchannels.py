import ROOT
import os
import yaml
import numpy as np
from ROOT import TCanvas, TGraph, TLegend
import CMS_lumi, tdrstyle
ROOT.gROOT.SetBatch(ROOT.kTRUE)

tag = 'FinalProd'
period = 'All' # 'All' if all year combined, else '2018', ...

# PLOT upper limits
def PlotExpLimits(V2_lim, outputfile, period):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
 
    # CMS style
    CMS_lumi.cmsText = "CMS simulation/Private work"
    CMS_lumi.extraText = "V_{e N}:V_{#mu N}:V_{#tau N} = 0:0:1"
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

    # Adding EXO-22-011 expected limit
    mass_EXO_22_011 = [20.16456139061248, 26.8352346902974,36.29431136096506, 49.61356550140213,58.614110947900635, 74.59242495694417, 83.86810310829574, 98.11748111433404,137.82937927094517, 175.65989929607008, 233.7480085567492, 338.99815299846193,491.5926903904981, 720.589843462258, 928.0560570265741]
    lim_EXO_22_011 =  [0.0005019100441728138,0.0005122034362867965,0.0005420013832510946,0.0007123808769013648,0.0015507257967259716,0.023750341304095056,0.27220281774884153,0.10849934664589993,0.08759253331025571,0.10899607477957826,0.1308625474466631,0.22973062550925316,0.4744717187020259,1.0160725186361823,2.096582937391052]
    EXO_22_011 = TGraph(len(mass_EXO_22_011))
    for j in range(len(mass_EXO_22_011)):
        EXO_22_011.SetPoint(j,mass_EXO_22_011[j],lim_EXO_22_011[j])

    # Adding ttl limit
    mass_ttl_limit = [
        40, 50, 60, 70, 75, 85, 100, 125, 150,
        200, 250, 300, 350, 400, 450, 500, 600, 
        800, 900, 1000
    ]

    lim_ttl_limit = [
        0.025217828737826516,   # 40
        0.0212220413989603,     # 50
        0.031130621442366026,   # 60
        0.055003546923927754,   # 70
        0.05523438425025445,    # 75
        0.1771370184425153,     # 85
        0.06701337685216789,    # 100
        0.0532670441765728,     # 125
        0.062417585704402825,   # 150
        0.08455804171252534,    # 200
        0.091012348555452,      # 250
        0.10532761100125833,    # 300
        0.13524945209733263,    # 350
        0.17918577506971164,    # 400
        0.1651330420969023,     # 450
        0.2925241133477115,     # 500
        0.37160267305315205,    # 600 (only one value)
        1.2116160049627793,     # 800
        2.043302761891477,      # 900
        3.3390557328470485      # 1000
    ]
    ttl_limit = TGraph(len(mass_ttl_limit))
    for j in range(len(mass_ttl_limit)):
        ttl_limit.SetPoint(j,mass_ttl_limit[j],lim_ttl_limit[j])

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
    #c.SetGrid()
    c.cd()
    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("|V_{#tau N}|^{2}")
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

    # Drawing EXO-22-011 expected limit
    EXO_22_011.SetLineColor(ROOT.kRed)
    EXO_22_011.SetLineWidth(2)
    EXO_22_011.SetLineStyle(1)
    EXO_22_011.Draw('Lsame')

    # Drawing ttl expected limit
    ttl_limit.SetLineColor(ROOT.kBlue)
    ttl_limit.SetLineWidth(2)
    ttl_limit.SetLineStyle(2)
    ttl_limit.Draw('Lsame')

    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.5
    x2 = x1 + 0.24
    y1 = 0.20
    y2 = y1 + 0.25
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)

    if period in ['2018','2017','2016', '2016_HIPM']:
        legend.AddEntry(median, f"Run2 {period} expected",'L')
    else:
        legend.AddEntry(median, "Asymptotic CL expected",'L')
    legend.AddEntry(ttl_limit, "#tau #tau l channels only",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
    legend.AddEntry(EXO_22_011, "EXO-22-011 CL expected",'L')
    
    legend.Draw()

    print(" ")
    c.SaveAs(outputfile)
    c.Close()

def main():

    outputdir = f'{os.getenv("RUN_PATH")}/F_Exp_limit/figures/{tag}/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    CombinePath = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}/'

    if period in ['2018','2017','2016', '2016_HIPM']:
        period_tag = '_' + period
    else:
        period_tag = '_all_years'

    fig_name =  f'ExpLimit_BDV_all_channels{period_tag}.pdf'
    input_file_name = f'intPoints_all_channels{period_tag}.yaml'

    outputfile = outputdir + fig_name
    input_file_path = os.path.join(CombinePath, input_file_name)

    with open(input_file_path, 'r') as yaml_file:
        intPoints = yaml.load(yaml_file, Loader=yaml.FullLoader)

    V2_lim = {}
    for HNL_mass in list(intPoints.keys()):
        V2_lim[HNL_mass] = np.array(intPoints[HNL_mass])
        if len(V2_lim[HNL_mass]) == 5:
            print(f'HNL_mass: {HNL_mass}')
            #print(f'V2_lim: {V2_lim[HNL_mass]}')
        else:
            print(f'HNL_mass: {HNL_mass}')
            print(f'!!!!!! V2_lim: {V2_lim[HNL_mass]} !!!!!!')
            print(" ")
            V2_lim[HNL_mass] = np.array([V2_lim[HNL_mass], V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass],V2_lim[HNL_mass]])
    
    PlotExpLimits(V2_lim, outputfile, period)
 
if __name__ == '__main__':
    main()
