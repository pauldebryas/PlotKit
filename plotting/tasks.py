#!/usr/bin/env python
import law
import luigi
import os
from plotting.Plotter import Plotter, load
import yaml
from law_customizations import Task, HTCondorWorkflow

class RunPlot(Task, HTCondorWorkflow, law.LocalWorkflow):
    '''
    Produce plots  
    '''
    period = luigi.Parameter()
    channel = luigi.Parameter()
    tag = luigi.Parameter()
    BCKestimation = luigi.Parameter(default='FakeRate')
    OutputTag = luigi.Parameter(default='')
    
    def load_args(self):
        if self.channel not in ['ttm', 'tmm', 'tte', 'tee', 'tem']:
            raise "channel not valid"

        if self.period not in ['2018','2017','2016','2016_HIPM']:
            raise "period not valid"

        if self.BCKestimation not in ['ABCD' ,'FakeRate' ,'GBReweighter', 'None']:
            raise "BCKestimation not valid"
        
        output_tag = self.tag + '_' + self.OutputTag

        self.input_folder = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', self.period, self.tag, self.channel, 'anatuple')
        #common config
        self.page_cfg = f'{os.getenv("RUN_PATH")}/config/shared/cms_stacked.yaml'
        self.page_cfg_custom = f'{os.getenv("RUN_PATH")}/config/shared/{self.period}.yaml'
        self.hist_maker_file = f'{os.getenv("RUN_PATH")}/hist_makers/hnl_maker_{self.BCKestimation}.py'
        #channel config
        self.hist_cfg = f'{os.getenv("RUN_PATH")}/config/hnl/hnl_{self.channel}/histograms.yaml'
        self.inputs_cfg = f'{os.getenv("RUN_PATH")}/config/hnl/hnl_{self.channel}/inputs_{self.BCKestimation}.yaml'
        self.output_dir = f'{os.getenv("RUN_PATH")}/output/{self.period}/{output_tag}/{self.channel}/'

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if self.channel == 'ttm':
            self.custom_title ="cat_text=#tau_{h}#tau_{h}#mu"
            self.var_to_plot = ['Tau1_pt', 'Tau2_pt', 'Muon_pt', 'pt_sum_ttm']

        if self.channel == 'tmm':
            self.custom_title ="cat_text=#tau_{h}#mu#mu"
            self.var_to_plot = ['Tau_pt', 'Muon1_pt', 'Muon2_pt', 'pt_sum_tmm', 'Tau_Jet_pt', 'Tau_decayMode']

        if self.channel == 'tte':
            self.custom_title ="cat_text=#tau_{h}#tau_{h} e"
            self.var_to_plot = ['Tau1_pt', 'Tau2_pt', 'Electron_pt', 'pt_sum_tte']

        if self.channel == 'tee':
            self.custom_title ="cat_text=#tau_{h} e e"
            self.var_to_plot = ['Tau_pt', 'Electron1_pt', 'Electron2_pt', 'pt_sum_tee', 'Tau_Jet_pt', 'Tau_decayMode']

        if self.channel == 'tem':
            self.custom_title ="cat_text=#tau_{h} e #mu"
            self.var_to_plot = ['Tau_pt', 'Electron_pt', 'Muon_pt', 'pt_sum_tem', 'Tau_Jet_pt']

        return
    
    def create_branch_map(self):
        self.load_args()
        branches = {}
        for i in range(len(self.var_to_plot)):
            branches[i] = os.path.join(self.output_dir, f'{self.var_to_plot[i]}.pdf')
        return branches
    
    def output(self):
        path = os.path.join(self.branch_data)
        return law.LocalFileTarget(path)
    
    def run(self):
        self.load_args()

        with open(self.inputs_cfg, 'r') as f:
            inputs = yaml.safe_load(f)
        for input in inputs:
            files = input.get('files')
            if files != None:
                for file in files:
                    if not os.path.isfile(os.path.join(self.input_folder, file)):
                        print(f'missing file: {file}')
                        raise(f'missing file: {file}')

        custom = None if self.custom_title is None else dict(item.split('=') for item in self.custom_title.split(','))
        plotter = Plotter(page_cfg=self.page_cfg, page_cfg_custom=self.page_cfg_custom, hist_cfg=self.hist_cfg, inputs_cfg=self.inputs_cfg)
        hist_maker = load(self.hist_maker_file)

        var = self.var_to_plot[self.branch]
        print("---------------------- Plotting "+var + "----------------------")
        hists = hist_maker.make_histograms(self.input_folder, hist_name=var, hist_cfg=plotter.hist_cfg, inputs_cfg=plotter.inputs_cfg, channel = self.channel)
        plotter.plot(var, hists, self.output_dir+f'{var}.pdf', custom=custom)

        return
