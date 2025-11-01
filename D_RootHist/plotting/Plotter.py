import copy
import importlib
import os
import ROOT
import sys
import yaml

def mk_smart_hist(hist, hist_desc):
  use_log_y = hist_desc.get('use_log_y', False)
  max_y_sf = hist_desc.get('max_y_sf', 1.5)
  divide_by_bin_width = hist_desc.get('divide_by_bin_width', False)
  is_flow = hist_desc.get('flow', False)
  
  if is_flow:
    # Get number of bins
    nbins = hist.GetNbinsX()
    # Add overflow bin content to the last bin
    last_bin_content = hist.GetBinContent(nbins)
    overflow_content = hist.GetBinContent(nbins + 1)  # Overflow bin
    hist.SetBinContent(nbins, last_bin_content + overflow_content)
    hist.SetBinContent(nbins + 1, 0)  # Clear overflow bin'
    # Add underflow bin content to the 0 bin
    first_bin_content = hist.GetBinContent(1)
    underflow_content = hist.GetBinContent(0)  # underflow bin
    hist.SetBinContent(1, first_bin_content + underflow_content)
    hist.SetBinContent(0, 0)  # Clear underflow bin'

  smart_hist = ROOT.root_ext.SmartHistogram('TH1D')(hist, use_log_y, max_y_sf, divide_by_bin_width)
  if 'x_title' in hist_desc:
    smart_hist.GetXaxis().SetTitle(hist_desc['x_title'])
  if 'y_title' in hist_desc:
    smart_hist.GetYaxis().SetTitle(hist_desc['y_title'])

  return smart_hist

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

class Plotter(object):
  initialized = False

  def __init__(self, page_cfg, page_cfg_custom, hist_cfg, inputs_cfg):
    if not Plotter.initialized:
      ROOT.gROOT.SetBatch(True)
      header_dir = os.path.join(os.getenv("RUN_PATH"), "D_RootHist", "plotting", "include")
      ROOT.gInterpreter.Declare(f'#include "{header_dir}/PdfPrinter.h"')
      ROOT.gInterpreter.Declare(f'#include "{header_dir}/StackedPlotDescriptor.h"')
      Plotter.initialized = True

    with open(page_cfg, 'r') as f:
      self.page_cfg = yaml.safe_load(f)
    if page_cfg_custom:
      with open(page_cfg_custom, 'r') as f:
        self.page_cfg.update(yaml.safe_load(f))
    self.hist_cfg = hist_cfg
    #with open(hist_cfg, 'r') as f:
    #  self.hist_cfg = yaml.safe_load(f)
    #with open(inputs_cfg, 'r') as f:
    #  self.inputs_cfg = yaml.safe_load(f)
    self.inputs_cfg = inputs_cfg

    self.sgn_hist_opt = LoadHistOptions(self.page_cfg['sgn_hist'])
    self.bkg_hist_opt = LoadHistOptions(self.page_cfg['bkg_hist'])
    self.data_hist_opt = LoadHistOptions(self.page_cfg['data_hist'])
    self.bkg_unc_hist_opt = LoadHistOptions(self.page_cfg['bkg_unc_hist'])
    self.page = LoadPageOptions(self.page_cfg['page_setup'])


  def plot(self, hist_name, histograms, output_file, custom=None, HNL_mass='HNL300', Unblind_data=False, PlotSignalOff=False):

    page_cfg = copy.deepcopy(self.page_cfg)
    if custom:
      for key, value in custom.items():
        page_cfg[key]['text'] = value

    items = ROOT.root_ext.draw_options.ItemCollection()
    for element in page_cfg['page_setup']['text_boxes'] + [ page_cfg['page_setup']['legend'] ]:
      items[element] = ToItem(page_cfg[element])
    printer = ROOT.root_ext.PdfPrinter(output_file, items, self.page)

    desc = ROOT.root_ext.StackedPlotDescriptor(self.page, self.sgn_hist_opt, self.bkg_hist_opt, self.data_hist_opt,
                                               self.bkg_unc_hist_opt)
    smart_hists = {}
    for input in self.inputs_cfg:
      name = input['name']
      hist_type = input.get('type', 'background')
      if (hist_type == 'signal'):
        if name == HNL_mass:
          smart_hists[name] = mk_smart_hist(histograms[name], self.hist_cfg[hist_name])
          if not PlotSignalOff:
            desc.AddSignalHistogram(smart_hists[name], input['title'], ROOT.root_ext.Color.Parse(input['color']),
                                    input.get('scale', 1.))
      elif hist_type == 'background':
        smart_hists[name] = mk_smart_hist(histograms[name], self.hist_cfg[hist_name])
        desc.AddBackgroundHistogram(smart_hists[name], input['title'], ROOT.root_ext.Color.Parse(input['color']))
      elif hist_type == 'data':
        smart_hists[name] = mk_smart_hist(histograms[name], self.hist_cfg[hist_name])
        if Unblind_data:
          desc.AddDataHistogram(smart_hists[name], input['title'])
      else:
        raise RuntimeError(f'Unknown histogram type: {hist_type}')
    printer.Print(hist_name, desc, True)

def load(module_file):
  if not os.path.exists(module_file):
    raise RuntimeError(f"Cannot find path to {module_file}.")

  module_name, module_ext = os.path.splitext(module_file)
  spec = importlib.util.spec_from_file_location(module_name, module_file)
  module = importlib.util.module_from_spec(spec)
  sys.modules[module_name] = module
  spec.loader.exec_module(module)
  return module

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Plotting tool.')
  parser.add_argument('--page-cfg', required=True, type=str, help="Page style config file")
  parser.add_argument('--page-cfg-custom', required=False, default=None, type=str,
                      help="Additional page style customisations")
  parser.add_argument('--hist-cfg', required=True, type=str, help="Histogram style config file")
  parser.add_argument('--inputs-cfg', required=True, type=str, help="Inputs style config file")
  parser.add_argument('--hist-name', required=True, type=str, help="Histogram name")
  parser.add_argument('--custom', required=False, default=None, type=str, help="Customizations in format key1=value1,key2=value2,...")
  parser.add_argument('--output', required=True, type=str, help="Output pdf file.")
  parser.add_argument('--hist-maker', required=True, type=str, help="Path to histogram maker module.")
  parser.add_argument('--verbose', required=False, type=int, default=0, help="verbosity level")
  parser.add_argument('hist_maker_args', type=str, nargs='*', help="hist maker arguments")
  args = parser.parse_args()


  plotter = Plotter(page_cfg=args.page_cfg, page_cfg_custom=args.page_cfg_custom, hist_cfg=args.hist_cfg,
                    inputs_cfg=args.inputs_cfg)
  hist_maker = load(args.hist_maker)
  hists = hist_maker.make_histograms(*args.hist_maker_args, hist_name=args.hist_name, hist_cfg=plotter.hist_cfg)
  custom = None if args.custom is None else dict(item.split('=') for item in args.custom.split(','))
  plotter.plot(args.hist_name, hists, args.output, custom=custom)
