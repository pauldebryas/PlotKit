import os
import ROOT
import yaml
ROOT.gROOT.SetBatch(True)

python_dir = os.path.dirname(os.path.abspath(__file__))
header_dir = os.path.join(python_dir, "include")

ROOT.gInterpreter.Declare(f'#include "{header_dir}/PdfPrinter.h"')
ROOT.gInterpreter.Declare(f'#include "{header_dir}/StackedPlotDescriptor.h"')

def mk_smart_hist(hist):
  use_log_y = False
  max_y_sf = 1.5
  divide_by_bin_width = False
  return ROOT.root_ext.SmartHistogram('TH1D')(hist, use_log_y, max_y_sf, divide_by_bin_width)
inputFile = ROOT.TFile.Open('hh_ttbb_input.root', 'READ')
hist_signal = mk_smart_hist(inputFile.Get('hh_ttbb_2018_tauTau_0_13TeV/qqHH_CV_1_C2V_1_kl_1_ttbb'))
hist_DY = mk_smart_hist(inputFile.Get('hh_ttbb_2018_tauTau_0_13TeV/DY_0b_100JPt'))
hist_TT = mk_smart_hist(inputFile.Get('hh_ttbb_2018_tauTau_0_13TeV/TT'))
hist_data = mk_smart_hist(inputFile.Get('hh_ttbb_2018_tauTau_0_13TeV/data_obs'))

outputFile = 'output.pdf'



with open('config/cms_stacked.yaml', 'r') as f:
  config = yaml.safe_load(f)

def to_str(value):
  if isinstance(value, bool):
    return 'true' if value else 'false'
  if isinstance(value, list):
    return ','.join(map(to_str, value))
  return str(value)

def ToItem(entry):
  map = ROOT.std.map('std::string', 'std::string')()
  for key, value in entry.items():
    map[key] = to_str(value)
  item = ROOT.root_ext.draw_options.Item()
  item.properties = ROOT.analysis.PropertyList(map)
  return item

def LoadHistOptions(entry):
  item = ToItem(entry)
  return ROOT.root_ext.draw_options.Histogram(item)

def LoadPageOptions(entry):
  item = ToItem(entry)
  return ROOT.root_ext.draw_options.Page(item)

sgn_hist_opt = LoadHistOptions(config['sgn_hist'])
bkg_hist_opt = LoadHistOptions(config['bkg_hist'])
data_hist_opt = LoadHistOptions(config['data_hist'])
bkg_unc_hist_opt = LoadHistOptions(config['bkg_unc_hist'])

page = LoadPageOptions(config['page_setup'])

items = ROOT.root_ext.draw_options.ItemCollection()
for element in config['page_setup']['text_boxes'] + [ config['page_setup']['legend'] ]:
  items[element] = ToItem(config[element])

printer = ROOT.root_ext.PdfPrinter(outputFile, items, page)



desc = ROOT.root_ext.StackedPlotDescriptor(page, sgn_hist_opt, bkg_hist_opt, data_hist_opt, bkg_unc_hist_opt)
desc.AddSignalHistogram(hist_signal, 'HH', ROOT.root_ext.Color(ROOT.kRed), 100)
desc.AddBackgroundHistogram(hist_DY, 'DY', ROOT.root_ext.Color(ROOT.kBlue))
desc.AddBackgroundHistogram(hist_TT, 'TT', ROOT.root_ext.Color(ROOT.kGreen))
desc.AddDataHistogram(hist_data, 'Data')

printer.Print('MT2', desc, True)