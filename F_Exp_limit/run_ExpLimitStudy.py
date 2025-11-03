import ROOT
import os
import numpy as np
import yaml
from array import array
ROOT.gROOT.SetBatch(ROOT.kTRUE)

tag = 'FinalProd'
period = '2018' # 'All' if all year combined, else '2018', ...

# MAIN
def main():
    # channels = ['tee','tmm','tte','ttm','tem']
    channels = ['tee_ss','tee_os','tmm_ss','tmm_os']
    BDV = {} #BDV[HNL_mass][channel]

    for channel in channels:
        save_dir = f"F_Exp_limit/figures/{tag}/"
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        save_BDV = f"F_Exp_limit/results/{tag}/"
        if not os.path.exists(save_BDV):
            os.makedirs(save_BDV)

        if period in ['2018','2017','2016', '2016_HIPM']:
            figure_name = f"ExpLimitStudy_{channel}_{period}.pdf"
            BDV_file_name = f'BDV_inputs_{period}.yaml'
            input_file = f'intPoints_{channel}_{period}.yaml'
        elif period == 'All':
            figure_name = f"ExpLimitStudy_{channel}_all_years.pdf"
            BDV_file_name = 'BDV_inputs.yaml'
            input_file = f'intPoints_{channel}.yaml'
        else:
            print('period parmeter unknown')

        save_res_file = f'/afs/cern.ch/user/p/pdebryas/HNL_analysis/Analysis/LimitEstimation/CMSSW_11_3_4/src/HNLAnalysis/results/intPoints/{tag}/' + input_file
        # Read the existing YAML content from the file
        with open(save_res_file, 'r') as file:
            yaml_content = yaml.safe_load(file)

        var_to_compare = list(yaml_content.keys())

        HNL_mass_range = list(yaml_content[var_to_compare[0]].keys()) 

        # Canvas with log scale for y-axis
        canvas = ROOT.TCanvas("canvas", "Log Scale Y-Axis", 800,600)
        canvas.SetLogy()
        legend = ROOT.TLegend(0.5, 0.2, 0.8, 0.4)

        # Values for plot
        x_values = []
        y_values = {var: [] for var in var_to_compare}

        for HNL_mass in HNL_mass_range:
            # Check if all required variables are not None for this mass
            if all(yaml_content[var].get(HNL_mass) is not None for var in var_to_compare):
                x_values.append(HNL_mass)
                for var in var_to_compare:
                    y_values[var].append(yaml_content[var][HNL_mass])
        
        #here
        x_values = []
        y_values = {var: [] for var in var_to_compare}

        for HNL_mass in HNL_mass_range:
            valid = True
            values = {}

            for var in var_to_compare:
                val = yaml_content[var].get(HNL_mass)
                if val is None:
                    valid = False
                    break
                elif len(val) == 1:
                    values[var] = val[0]
                elif len(val) == 5:
                    values[var] = val[2]
                else:
                    print(f'Error for HNL mass {HNL_mass} and var {var}: Unexpected length {len(val)}')
                    values[var] = val[2] if len(val) > 2 else None
                    if values[var] is None:
                        valid = False
                        break

            if valid:
                x_values.append(HNL_mass)
                for var in var_to_compare:
                    y_values[var].append(values[var])

        '''
        x_values = HNL_mass_range
        y_values = {}
        for var in var_to_compare:
            y_values[var] = []
            for HNL_mass in HNL_mass_range:
                if yaml_content[var][HNL_mass] == None:
                    # remove mass value and do not add y value
                elif len(yaml_content[var][HNL_mass]) == 1:
                    y_values[var].append(yaml_content[var][HNL_mass][0])
                elif len(yaml_content[var][HNL_mass]) == 5:
                    y_values[var].append(yaml_content[var][HNL_mass][2])
                else:
                    print(f'error for HNL mass {HNL_mass} and var {var}')
                    y_values[var].append(yaml_content[var][HNL_mass][2])
        '''

        # Create a TGraph
        var = var_to_compare[0]
        graph1 = ROOT.TGraph(len(x_values), array('d', x_values), array('d', y_values[var]))
        graph1.SetLineStyle(1)
        graph1.SetLineColor(1)
        graph1.Draw("AC")  # A: Axis, P: Points (markers)
        legend.AddEntry(graph1, var, "l")

        # Create a TGraph for the second variable
        var = var_to_compare[1]
        graph2 = ROOT.TGraph(len(x_values), array('d', x_values), array('d', y_values[var]))
        graph2.SetLineStyle(1)
        graph2.SetLineColor(2)
        graph2.Draw("C SAME")  # C: Curve, SAME: Draw on the same canvas
        legend.AddEntry(graph2, var, "l")  # Add entry to the legend for graph2

        # Create a TGraph
        var = var_to_compare[2]
        graph3 = ROOT.TGraph(len(x_values), array('d', x_values), array('d', y_values[var]))
        graph3.SetLineStyle(1)
        graph3.SetLineColor(3)
        graph3.Draw("C SAME") 
        legend.AddEntry(graph3, var, "l")

        # Create a TGraph
        var = var_to_compare[3]
        graph4 = ROOT.TGraph(len(x_values), array('d', x_values), array('d', y_values[var]))
        graph4.SetLineStyle(1)
        graph4.SetLineColor(4)
        graph4.Draw("C SAME") 
        legend.AddEntry(graph4, var, "l")

        #Create a TGraph
        var = var_to_compare[4]
        graph5 = ROOT.TGraph(len(x_values), array('d', x_values), array('d', y_values[var]))
        graph5.SetLineStyle(1)
        graph5.SetLineColor(5)
        graph5.Draw("C SAME") 
        legend.AddEntry(graph5, var, "l")

        canvas.Update()

        # Customize axis labels and legend
        graph1.GetXaxis().SetTitle('m_{HNL} [GeV]') 
        graph1.GetYaxis().SetTitle("|V_{N#tau}|^{2}")
        graph1.GetYaxis().SetRangeUser(5e-3, 1.)
        graph1.GetXaxis().SetRangeUser(20, 700)
        legend.SetBorderSize(1)  # No border
        legend.SetFillColor(0)   # Transparent background
        legend.Draw()

        canvas.Show()
        # Show the canvas

        canvas.SaveAs(save_dir + figure_name)


        BDV[channel] = {}
        i = 0
        for mass in x_values:
            y_values_available = []
            var_available = []
            for var in y_values.keys():
                y_values_available.append(y_values[var][i]) 
                var_available.append(var)  
            BDV[channel][mass] = var_available[np.argmin(np.array(y_values_available))]
            i = i+1
 
    
    with open(save_BDV + BDV_file_name, 'w') as file:
        yaml.dump(BDV, file)

if __name__ == '__main__':
    main()
