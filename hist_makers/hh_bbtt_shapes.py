import ROOT

def Get(file, hist_name):
  hist = file.Get(hist_name)
  if not hist:
    raise RuntimeError(f'Cannot find histogram {hist_name}')
  hist = hist.Clone()
  hist.SetDirectory(0)
  return hist

def make_histograms(input_file, hist_name=None, hist_cfg=None):
  inputFile = ROOT.TFile.Open(input_file, 'READ')
  hists = {}
  hists['qqHH'] = Get(inputFile, 'hh_ttbb_2018_tauTau_0_13TeV/qqHH_CV_1_C2V_1_kl_1_ttbb')
  for n_b in range(3):
    for pt_bin in [ 0, 10, 30, 50, 100, 200 ]:
      hist = Get(inputFile, f'hh_ttbb_2018_tauTau_0_13TeV/DY_{n_b}b_{pt_bin}JPt')
      if 'DY' not in hists:
        hists['DY'] = hist
      else:
        hists['DY'].Add(hist)
  hists['W'] = Get(inputFile, 'hh_ttbb_2018_tauTau_0_13TeV/W')
  hists['TT'] = Get(inputFile, 'hh_ttbb_2018_tauTau_0_13TeV/TT')
  hists['QCD'] = Get(inputFile, 'hh_ttbb_2018_tauTau_0_13TeV/QCD')
  hists['data'] = Get(inputFile, 'hh_ttbb_2018_tauTau_0_13TeV/data_obs')
  return hists
