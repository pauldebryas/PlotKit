import hist
import os
import uproot
import ROOT

def ListToVector(list, type="string"):
  vec = ROOT.std.vector(f"{type}")()
  for item in list:
    vec.push_back(item)
  return vec

def ToRootHist(h):
  x_bins = ListToVector(h.axes[0].edges, type="double")
  root_hist = ROOT.TH1D('', '', len(x_bins) - 1, x_bins.data())
  for i in range(len(x_bins) - 1):
    root_hist.SetBinContent(i + 1, h.values()[i])
    root_hist.SetBinError(i + 1, h.variances()[i] ** 0.5)
  root_hist.SetDirectory(0)
  return root_hist

def make_histograms(input_dir, hist_name=None, hist_cfg=None):
  inputs = {
    'DY' : [
      'DYJetsToLL_M-50_anatuple.root',
      'DYJetsToLL_0J_anatuple.root',
      'DYJetsToLL_1J_anatuple.root',
      'DYJetsToLL_2J_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-0To50_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-50To100_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-100To250_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-250To400_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-400To650_anatuple.root',
      'DYJetsToLL_LHEFilterPtZ-650ToInf_anatuple.root',
    ],
    'TT' : [
      'TTTo2L2Nu_anatuple.root',
      'TTToSemiLeptonic_anatuple.root',
      'TTToHadronic_anatuple.root',
    ],
    'W': [
      'WJetsToLNu_anatuple.root',
      'W1JetsToLNu_anatuple.root',
      'W2JetsToLNu_anatuple.root',
      'W3JetsToLNu_anatuple.root',
      'W4JetsToLNu_anatuple.root',
      'WJetsToLNu_HT-70To100_anatuple.root',
      'WJetsToLNu_HT-100To200_anatuple.root',
      'WJetsToLNu_HT-200To400_anatuple.root',
      'WJetsToLNu_HT-400To600_anatuple.root',
      'WJetsToLNu_HT-600To800_anatuple.root',
      'WJetsToLNu_HT-800To1200_anatuple.root',
      'WJetsToLNu_HT-1200To2500_anatuple.root',
      'WJetsToLNu_HT-2500ToInf_anatuple.root',
    ],
    'HNL': [ 'HNL_tau_M-300_anatuple.root' ],
    'data': [ 'SingleMuon_2018_anatuple.root' ]
  }

  hists = {}
  hist_desc = hist_cfg[hist_name]
  branch = hist_name
  for input_name, input_files in inputs.items():
    hists[input_name] = hist.Hist.new.Variable(hist_desc['x_bins'], name='x').Weight()
    for file_name in input_files:
      file = uproot.open(os.path.join(input_dir, file_name))
      tree = file['Event']
      hists[input_name].fill(x=tree[branch].array(), weight=tree['genWeight'].array())
  root_hists = { hist_name : ToRootHist(h) for hist_name, h in hists.items() }
  return root_hists