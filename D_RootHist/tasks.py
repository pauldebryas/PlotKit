#!/usr/bin/env python
import law
import luigi
import os
import yaml

from D_RootHist.plotting.Plotter import Plotter, load
from D_RootHist.hist_makers.helpers import make_root_file, load_root_file
from law_customizations import Task, HTCondorWorkflow
from common.helpers import get_hnl_masses 

def load_inputs(BCKestimationMethod, period, channel):
    
    inputs = []
    config_path_BCKG = os.path.join(os.getenv("RUN_PATH"),"common","config", "all", "inputs")

    #channel config files
    if BCKestimationMethod in ['FakeRate','ValidationMuFR', 'ValidationEleFR', 'ValidationTauFR']:
        #load TrueLepton background
        with open(os.path.join(config_path_BCKG, 'inputs_TrueLepton.yaml'), 'r') as f:
            MC_inputs = yaml.safe_load(f)
            if BCKestimationMethod in ['ValidationMuFR', 'ValidationEleFR']:
                MC_inputs[0]['files'].remove('ZHToTauTau_anatuple.root')
        inputs.append(MC_inputs[0])

        #load Fake Background
        with open(os.path.join(config_path_BCKG, 'inputs_FakeBackground.yaml'), 'r') as f:
            FB_inputs = yaml.safe_load(f)
        inputs.append(FB_inputs[0])

        #load signal
        if BCKestimationMethod == 'FakeRate':
            with open(os.path.join(config_path_BCKG, 'inputs_AllSignal.yaml'), 'r') as f:
                signal_inputs = yaml.safe_load(f)
            for i in range(len(signal_inputs)):
                inputs.append(signal_inputs[i])

        #load data
        with open(f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_{channel}/inputs_data.yaml', 'r') as f:
            data_inputs = yaml.safe_load(f)
        inputs.append(data_inputs[0])

    if BCKestimationMethod in ['MonteCarlo']:
        #load MC background
        with open(os.path.join(config_path_BCKG, f'inputs_MainMCbackground_{channel}.yaml'), 'r') as f:
            MC_inputs = yaml.safe_load(f)
        for i in range(len(MC_inputs)):
            inputs.append(MC_inputs[i])

        #load signal
        with open(os.path.join(config_path_BCKG, 'inputs_AllSignal.yaml'), 'r') as f:
            signal_inputs = yaml.safe_load(f)
        for i in range(len(signal_inputs)):
            inputs.append(signal_inputs[i])

        #load data
        with open(f'{os.getenv("RUN_PATH")}/common/config/{period}/hnl/hnl_{channel}/inputs_data.yaml', 'r') as f:
            data_inputs = yaml.safe_load(f)
        inputs.append(data_inputs[0])

    if len(inputs) == 0:
        raise (f'{BCKestimationMethod} not inplemented in load_inputs')
    
    return inputs

class MakeTH1Hist(Task, HTCondorWorkflow, law.LocalWorkflow):
    '''
    Produce TH1Hist for limit computation and plotting
    '''
    period = luigi.Parameter()
    channel = luigi.Parameter()
    tag = luigi.Parameter()
    BCKestimation = luigi.Parameter()
    PlotRegion = luigi.Parameter() # ['SignalRegion', 'InvertedBjetsVetoRegion', 'ttControlRegion', 'ttClosureRegion', 'DYClosureRegion', 'DYControlRegion']
    
    def load_args(self):

        if self.period not in ['2018','2017','2016','2016_HIPM']:
            raise "period parameter not valid"
        
        if self.channel not in ['ttm', 'tmm', 'tte', 'tee', 'tem', 'Zmu', 'Ze', 'tll', 'llmu', 'lle']:
            raise "channel parameter not valid"
        
        if self.channel in ['Zmu', 'Ze', 'llmu', 'lle']:
            anatuple_path = os.path.join('/eos/user/p/pdebryas/HNL_LLFF/anatuple', self.period)
        else:
            anatuple_path = os.path.join('/eos/user/p/pdebryas/HNL/anatuple', self.period)
        
        if self.tag not in os.listdir(anatuple_path):
            print(os.listdir(anatuple_path))
            raise "tag parameter not valid: must be within the list above"

        if f'hnl_maker_{self.BCKestimation}.py' not in os.listdir(f'{os.getenv("RUN_PATH")}/D_RootHist/hist_makers/'):
            print(os.listdir(f'{os.getenv("RUN_PATH")}/D_RootHist/hist_makers/'))
            raise "BCKestimation parameter not valid"
                                        
        #input folder where anatuple are stored
        self.input_folder = os.path.join(anatuple_path, self.tag, self.channel, 'anatuple')

        #BCKestimation method
        self.hist_maker_file = f'{os.getenv("RUN_PATH")}/D_RootHist/hist_makers/hnl_maker_{self.BCKestimation}.py'

        self.inputs = load_inputs(self.BCKestimation, self.period, self.channel)

        with open(f'{os.getenv("RUN_PATH")}/common/config/all/histograms/histograms_{self.channel}.yaml', 'r') as f:
            self.hist_cfg = yaml.safe_load(f)
        
        #output dir where to store TH1Hist
        self.output_dir = f'{os.getenv("RUN_PATH")}/D_RootHist/results/{self.period}/{self.tag}/{self.channel}/{self.BCKestimation}/{self.PlotRegion}/'
        os.makedirs(self.output_dir, exist_ok=True)

         # load var to plot
        if self.PlotRegion in ['SignalRegion', 'InvertedBjetsVetoRegion']:
            HNLMassRange = get_hnl_masses(self.period) #[20 ,30 ,40 ,50 ,60 ,70 ,75 ,85 ,100 ,125 ,150 ,200 ,250 ,300 ,350 ,400 ,450 ,500 ,600 ,700 ,800 ,900 ,1000]
            self.vars = []
            for var in list(self.hist_cfg.keys()):
                for HNLMass in HNLMassRange:
                    var_name = f'{var}_HNLMass{str(HNLMass)}'
                    self.vars.append(var_name)
        else:
            self.vars = list(self.hist_cfg.keys())

        return
    
    def create_branch_map(self):
        self.load_args()
        branches = {}
        for i in range(len(self.vars)):
            branches[i] = os.path.join(self.output_dir, f'TH1_{self.vars[i]}_HNL.root')
        #print(branches)
        return branches
    
    def output(self):
        path = os.path.join(self.branch_data)
        return law.LocalFileTarget(path)
    
    def run(self):
        self.load_args()
        #check no files is missing + remove unselected masses
        for input in self.inputs:
            files = input.get('files')
            if files != None:
                for file in files:
                    if not os.path.isfile(os.path.join(self.input_folder, file)):
                        print(f'missing file: {os.path.join(self.input_folder, file)}')
                        #raise(f'missing file')
        var = self.vars[self.branch]

        print("---------------------- Producing TH1Hists for "+ var + " in " + self.PlotRegion+ " ----------------------")
        hist_maker = load(self.hist_maker_file)
        hists = hist_maker.make_histograms(self.input_folder, hist_name=var, hist_cfg=self.hist_cfg, inputs_cfg=self.inputs, channel = self.channel, period=self.period, PlotRegion = self.PlotRegion, tag= self.tag)

        if self.PlotRegion == 'SignalRegion':
            nameDir = 'signal_region'
        else:
            nameDir = self.PlotRegion

        path_root_file = os.path.join(self.branch_data)

        make_root_file(path_root_file, nameDir, hists)

        return

class RunCMSPlot(Task, HTCondorWorkflow, law.LocalWorkflow):
    '''
    Produce CMS like plots  
    '''
    period = luigi.Parameter()
    channel = luigi.Parameter()
    tag = luigi.Parameter()
    BCKestimation = luigi.Parameter()
    PlotRegion = luigi.Parameter()
    UnblindData = luigi.BoolParameter(default=False)
    PlotSignalOff = luigi.BoolParameter(default=False)

    def workflow_requires(self):
        return { "TH1Hist": MakeTH1Hist.req(self, branch=self.branch) }

    def requires(self):
        return self.workflow_requires()
    
    def load_args(self):
        
        self.page_cfg = f'{os.getenv("RUN_PATH")}/common/config/all/cms_stacked.yaml'
        self.page_cfg_custom = f'{os.getenv("RUN_PATH")}/common/config/{self.period}/{self.period}_lumi.yaml'
        #channel config
        with open(f'{os.getenv("RUN_PATH")}/common/config/all/histograms/histograms_{self.channel}.yaml', 'r') as f:
            self.hist_cfg = yaml.safe_load(f)

        self.inputs = load_inputs(self.BCKestimation, self.period, self.channel)

        #output dir where to store CMShist
        self.output_dir = f'{os.getenv("RUN_PATH")}/D_RootHist/figures/{self.period}/{self.tag}/{self.channel}/{self.BCKestimation}/{self.PlotRegion}/'
        os.makedirs(self.output_dir, exist_ok=True)

         # load var to plot
        if self.PlotRegion in ['SignalRegion', 'InvertedBjetsVetoRegion']:
            HNLMassRange = HNLMassRange = get_hnl_masses(self.period)
            self.vars = []
            for var in list(self.hist_cfg.keys()):
                for HNLMass in HNLMassRange:
                    var_name = f'{var}_HNLMass{str(HNLMass)}'
                    self.vars.append(var_name)
        else:
            self.vars = list(self.hist_cfg.keys())
                
        if self.channel == 'ttm':
            self.custom_title ="cat_text=#tau_{h}#tau_{h}#mu"

        if self.channel == 'tmm':
            self.custom_title ="cat_text=#tau_{h}#mu#mu"

        if self.channel == 'tte':
            self.custom_title ="cat_text=#tau_{h}#tau_{h} e"

        if self.channel == 'tee':
            self.custom_title ="cat_text=#tau_{h} e e"

        if self.channel == 'tem':
            if self.PlotRegion == 'ttbarRegionFR':
                self.custom_title ="cat_text=CR_{#tau}^{t#bar{t}}"
            else:
                self.custom_title ="cat_text=#tau_{h} e #mu"

        if self.channel == 'Zmu':
            if self.PlotRegion == 'DYRegionFR':
                self.custom_title ="cat_text=CR_{#mu}^{DY}"
            elif self.PlotRegion == 'DYRegionValidation':
                self.custom_title ="cat_text=VR_{#mu}^{DY}"
            else:
                self.custom_title ="cat_text=Z#mu"

        if self.channel == 'Ze':
            if self.PlotRegion == 'DYRegionFR':
                self.custom_title ="cat_text=CR_{e}^{DY}"
            elif self.PlotRegion == 'DYRegionValidation':
                self.custom_title ="cat_text=VR_{e}^{DY}"
            else:
                self.custom_title ="cat_text=Ze"

        if self.channel == 'tll':
            if self.PlotRegion == 'ttbarRegionValidation':
                self.custom_title ="cat_text=VR_{#tau}^{t#bar{t}}"
            elif self.PlotRegion == 'DYRegionFR':
                self.custom_title ="cat_text=CR_{#tau}^{DY}"
            elif self.PlotRegion == 'DYRegionValidation':
                self.custom_title ="cat_text=VR_{#tau}^{DY}"
            else:
                self.custom_title ="cat_text=tll"

        if self.channel == 'llmu':
            if self.PlotRegion == 'ttbarRegionFR':
                self.custom_title ="cat_text=CR_{#mu}^{t#bar{t}}"
            elif self.PlotRegion == 'ttbarRegionValidation':
                self.custom_title ="cat_text=VR_{#mu}^{t#bar{t}}"
            else:
                self.custom_title ="cat_text=t#bar{t}"

        if self.channel == 'lle':
            if self.PlotRegion == 'ttbarRegionFR':
                self.custom_title ="cat_text=CR_{e}^{t#bar{t}}"
            elif self.PlotRegion == 'ttbarRegionValidation':
                self.custom_title ="cat_text=VR_{e}^{t#bar{t}}"
            else:
                self.custom_title ="cat_text=t#bar{t}"

        return
    
    def create_branch_map(self):
        self.load_args()
        branches = {}
        for i in range(len(self.vars)):
            branches[i] = os.path.join(self.output_dir, f'CMSPlot_{self.vars[i]}.pdf')
        #print(branches)
        return branches
    
    def output(self):
        path = os.path.join(self.branch_data)
        return law.LocalFileTarget(path)
    
    def run(self):
        self.load_args()
        var = self.vars[self.branch]
        if '_HNLMass' in var:
            MassHNL_Hyp = int(var.split('_HNLMass')[-1])
            self.hist_cfg[var] = self.hist_cfg[var.split('_HNLMass')[0]]
            if 'DNN' in var.split('_HNLMass')[0]:
                self.hist_cfg[var]['x_title'] = self.hist_cfg[var]['x_title'] + f'{MassHNL_Hyp}GeV'
        else:
            MassHNL_Hyp = '300'

        custom = None if self.custom_title is None else dict(item.split('=') for item in self.custom_title.split(','))

        plotter = Plotter(page_cfg=self.page_cfg, page_cfg_custom=self.page_cfg_custom, hist_cfg=self.hist_cfg, inputs_cfg=self.inputs)

        print("---------------------- Plotting " + var + " in " + self.PlotRegion+ " ----------------------")
        TH1Histfile = f'{os.getenv("RUN_PATH")}/D_RootHist/results/{self.period}/{self.tag}/{self.channel}/{self.BCKestimation}/{self.PlotRegion}/TH1_{var}_HNL.root'
        path_pdf_file = os.path.join(self.branch_data)

        if self.PlotRegion == 'SignalRegion':
            nameDir = 'signal_region'
        else:
            nameDir = self.PlotRegion
        
        hists = load_root_file(TH1Histfile, nameDir)
        plotter.plot(var, hists, path_pdf_file, custom=custom, HNL_mass=f'HNL{MassHNL_Hyp}', Unblind_data = self.UnblindData , PlotSignalOff = self.PlotSignalOff)
        return

