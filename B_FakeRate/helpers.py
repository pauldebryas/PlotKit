import hist
import correctionlib.convert
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt

from common.helpers import equalObs
from common.regions.regions import compute_region_mask

def AdaptBin(Lep_pt, n_bins):
  return np.array( equalObs( Lep_pt, n_bins) )

def print_nFake(inputs, channel, Lepton, Region):
    print('')
    print(f'for channel {channel}')
    Nevents_SR = 0
    N_fake = 0
    for MCsample in inputs[channel].keys():
        cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC', Region)
        Nevents_SR = Nevents_SR + np.sum(cut_region[Region])
        N_fake = N_fake + np.sum(cut_region[f'{Region}_Fake{Lepton}'])

    print(f'Nevents: {Nevents_SR}')
    print(f'N_fake: {N_fake}')
    print(f'(N_fake/Nevents_SR)*100: {(N_fake/Nevents_SR)*100}')
    return

def save_hist_schemav2(hist, desc, filename):
  cset = correctionlib.schemav2.CorrectionSet(
      schema_version=2,
      description=desc,
      corrections=[
          hist
      ],
  )
  with open(filename, "w") as fout:
      fout.write(cset.json(exclude_unset=True))
  return

def PrintContributionMC(name_lep, channels, inputs, ControlRegionName):

    #remove data from inputs
    for channel in channels:
        inputs[channel].pop('data')

    num_events = {}
    num_events_LepIsPromptLepton = {}
    # set num events per MC sample to 0
    for DataSample in inputs[channel].keys():
        num_events[DataSample] = 0
        num_events_LepIsPromptLepton[DataSample] = 0

    for channel in channels:
        for DataSample in inputs[channel].keys():
            cut_region = compute_region_mask(inputs[channel][DataSample], channel, 'MC', ControlRegionName)
            num_events[DataSample] += np.sum(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP']])
            num_events_LepIsPromptLepton[DataSample] += np.sum(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']])

    # Calculate total events
    total_events = sum(num_events.values())

    # Print the table header
    print(f"{'Category':<15} {'Total Events':<15} {'Fraction (%)':<15} {f'Events where {name_lep}IsPromptLepton':<30}")
    print("-" * 75)

    # Loop through the dictionary and print each category, event number, and fraction
    for category, events in num_events.items():
        fraction = (events / total_events) * 100
        print(f"{category:<15} {events:<15.2f} {fraction:<15.2f} {num_events_LepIsPromptLepton[category]:<30.2f}")

def create_histograms(dists, num, den, pt_bins, eta_bins, output_figures, pt_label=r'$p_t$ [GeV]'):
    err_num = np.where(num > 0, np.sqrt(num), 0)
    err_den = np.where(den > 0, np.sqrt(den), 0)

    # Set scale factor and its error
    sf = np.where((den > 0), num / den, 0)
    err_sf = np.where((num > 0) & (den > 0), sf * np.sqrt((err_num / num) ** 2 + (err_den / den) ** 2), 0)

    # Create histograms
    sfhist = hist.Hist(*dists.axes[1:], data=sf)
    err_sfhist = hist.Hist(*dists.axes[1:], data=err_sf)

    sfhist.name = "fake_rate"
    err_sfhist.name = "fake_rate_err"
    sfhist.label = "FR"
    err_sfhist.label = "FR_err"

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 10))
    sfhist.plot2d(ax=ax)

    pt_bins[0] = 10
    pt_bins[-1] = 100
    for i in range(len(eta_bins) - 1):
        for j in range(len(pt_bins) - 1):
            text_str = f"{round(sfhist.values().T[i, j], 2)} \n $\pm$ {round(err_sfhist.values().T[i, j], 2)}"
            x = pt_bins[j] + (pt_bins[j + 1] - pt_bins[j]) / 2
            y = eta_bins[i] + (eta_bins[i + 1] - eta_bins[i]) / 2
            ax.text(x, y, text_str, color="w", ha="center", va="center", fontweight="normal")

    ax.set(xlabel=pt_label, ylabel=r'$|\eta|$')
    ax.set_xscale('log')
    ax.set_xlim(left=10, right=100)
    plt.savefig(output_figures)

    return sfhist, err_sfhist

def make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures, ptlabel = r'$p_t$ [GeV]'):
    dists = (
        hist.Hist.new
        .StrCat(["num", "den"], name="selection", growth=True)
        .Variable(pt_bins, name= "pt")
        .Variable(eta_bins, name="abs_eta")
        .Weight()
        .fill(
            selection="num",
            pt=num_pt,
            abs_eta=np.abs(num_eta),
            weight=num_weight
        )
        .fill(
            selection="den",
            pt=den_pt,
            abs_eta=np.abs(den_eta),
            weight=den_weights
        )
    )
    num = np.array(dists["num", :, :].values())
    den = np.array(dists["den", :, :].values())

    return create_histograms(dists, num, den, pt_bins, eta_bins, output_figures, pt_label = ptlabel)

def compute_FakeRate(name_lep, inputs , bins, output_results, output_figures, channels, region):

  den_pt = []
  den_eta = []
  den_weights = []
  num_pt = []
  num_eta = []
  num_weight = []
  if len(channels) == 1 :
    channel = channels[0]
    name_fig = f'P_fake_TL_{name_lep}_{channel}'
    if f'{name_lep}_pt' in ak.fields(inputs[channel][list(inputs[channel].keys())[0]]):
      for MCsample in inputs[channel].keys():
        cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC', region)
        den_pt.append(inputs[channel][MCsample][f'{name_lep}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']]) 
        den_eta.append(inputs[channel][MCsample][f'{name_lep}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']]) 
        den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']])
        num_pt.append(inputs[channel][MCsample][f'{name_lep}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
        num_eta.append(inputs[channel][MCsample][f'{name_lep}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
        num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
    else:
      for i in ['1','2']:      
        for MCsample in inputs[channel].keys():
            cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC', region)
            den_pt.append(inputs[channel][MCsample][f'{name_lep}{i}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']]) 
            den_eta.append(inputs[channel][MCsample][f'{name_lep}{i}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']]) 
            den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']])
            num_pt.append(inputs[channel][MCsample][f'{name_lep}{i}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
            num_eta.append(inputs[channel][MCsample][f'{name_lep}{i}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
            num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
  else:
    name_fig = f'P_fake_TL_{name_lep}_MC'
    for channel in channels:
      if f'{name_lep}_pt' in ak.fields(inputs[channel][list(inputs[channel].keys())[0]]):
        for MCsample in inputs[channel].keys():
          cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC', region)
          den_pt.append(inputs[channel][MCsample][f'{name_lep}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']]) 
          den_eta.append(inputs[channel][MCsample][f'{name_lep}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']]) 
          den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}']])
          num_pt.append(inputs[channel][MCsample][f'{name_lep}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
          num_eta.append(inputs[channel][MCsample][f'{name_lep}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
          num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}_passTight']])
      else:
        for i in ['1','2']:      
          for MCsample in inputs[channel].keys():
              cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC', region)
              den_pt.append(inputs[channel][MCsample][f'{name_lep}{i}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']]) 
              den_eta.append(inputs[channel][MCsample][f'{name_lep}{i}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']]) 
              den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}']])
              num_pt.append(inputs[channel][MCsample][f'{name_lep}{i}_pt'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
              num_eta.append(inputs[channel][MCsample][f'{name_lep}{i}_eta'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
              num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region[f'{region}_PassLooseWP_Fake{name_lep}{i}_passTight']])
      
  den_pt = np.concatenate(den_pt)
  den_eta = np.concatenate(den_eta)
  den_weights = np.concatenate(den_weights)
  num_pt = np.concatenate(num_pt)
  num_eta = np.concatenate(num_eta)
  num_weight = np.concatenate(num_weight)

  den_pt = np.array(den_pt).flatten()
  den_eta = np.array(den_eta).flatten()
  den_weights = np.array(den_weights).flatten()
  num_pt = np.array(num_pt).flatten()
  num_eta = np.array(num_eta).flatten()
  num_weight = np.array(num_weight).flatten()

  output_figures = output_figures + name_fig + '.pdf'
  output_json = output_results + name_fig + '.json'
  output_json_err = output_results + name_fig + '_err.json'

  if isinstance(bins, int):
    pt_bins = np.array(equalObs(den_pt, bins))
    eta_bins = np.array(equalObs(np.abs(den_eta), bins))
  elif len(bins) == 2:
    pt_bins = np.array(bins[0])
    eta_bins = np.array(bins[1])
  else:
     print('pb with nbins')

  sfhist, err_sfhist = make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures)

  fake_rate = correctionlib.convert.from_histogram(sfhist)
  fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

  # set overflow bins behavior (default is to raise an error when out of bounds)
  fake_rate.data.flow = "clamp"
  fake_rate_err.data.flow = "clamp"

  save_hist_schemav2(fake_rate, f"fake rate {name_lep}", output_json)
  save_hist_schemav2(fake_rate_err, f"fake rate {name_lep} err", output_json_err)
  return

def compute_FakeRateData(name_lep, inputs ,bins, output_results, output_figures, channels, ControlRegionName):
    # computed in DY or ttbar region so in single Tau channel
    den_pt = []
    den_eta = []
    den_weights = []
    num_pt = []
    num_eta = []
    num_weight = []
    name_fig = f'P_fake_TL_{name_lep}'
    for channel in channels:
      if f'{name_lep}_pt' in ak.fields(inputs[channel][list(inputs[channel].keys())[0]]):
        for DataSample in inputs[channel].keys():
          if DataSample == 'data':
            cut_region = compute_region_mask(inputs[channel][DataSample], channel, 'data', ControlRegionName)
            den_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
            den_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
            den_weights.append(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP']])
            num_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP']])
            num_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP']])
            num_weight.append(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP']])
          else:
            cut_region = compute_region_mask(inputs[channel][DataSample], channel, 'MC', ControlRegionName)
            den_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']]) 
            den_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']]) 
            den_weights.append(-1*inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']])
            num_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])
            num_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])
            num_weight.append(-1*inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])            
      else:
        print(f'missing {name_lep}_pt in fields')
      
    den_pt = np.concatenate(den_pt)
    den_eta = np.concatenate(den_eta)
    den_weights = np.concatenate(den_weights)
    num_pt = np.concatenate(num_pt)
    num_eta = np.concatenate(num_eta)
    num_weight = np.concatenate(num_weight)

    den_pt = np.array(den_pt).flatten()
    den_eta = np.array(den_eta).flatten()
    den_weights = np.array(den_weights).flatten()
    num_pt = np.array(num_pt).flatten()
    num_eta = np.array(num_eta).flatten()
    num_weight = np.array(num_weight).flatten()

    output_figures = output_figures + name_fig + '.pdf'
    output_json = output_results + name_fig + '.json'
    output_json_err = output_results + name_fig + '_err.json'

    if isinstance(bins, int):
      pt_bins = np.array(equalObs(den_pt, bins))
      eta_bins = np.array(equalObs(np.abs(den_eta), bins))
    elif len(bins) == 2:
      pt_bins = np.array(bins[0])
      eta_bins = np.array(bins[1])
    else:
      print('pb with nbins')
      
    sfhist, err_sfhist = make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures)

    fake_rate = correctionlib.convert.from_histogram(sfhist)
    fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

    # set overflow bins behavior (default is to raise an error when out of bounds)
    fake_rate.data.flow = "clamp"
    fake_rate_err.data.flow = "clamp"

    save_hist_schemav2(fake_rate, f"fake rate {name_lep}", output_json)
    save_hist_schemav2(fake_rate_err, f"fake rate {name_lep} err", output_json_err)
    return


def compute_FakeRateMC(name_lep, inputs, pt_bins, eta_bins, output_results, output_figures, channels, ControlRegionName):
    # computed in DY or ttbar region so in single Tau channel

    den_pt = []
    den_eta = []
    den_weights = []
    num_pt = []
    num_eta = []
    num_weight = []
    name_fig = f'P_fake_TL_{name_lep}_MC'
    for channel in channels:
      if f'{name_lep}_pt' in ak.fields(inputs[channel][list(inputs[channel].keys())[0]]):
        for DataSample in inputs[channel].keys():
          if DataSample != 'data':
            cut_region = compute_region_mask(inputs[channel][DataSample], channel, 'MC', ControlRegionName)

            den_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
            den_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
            den_weights.append(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP']])
            num_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP']])
            num_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP']])
            num_weight.append(inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP']])

            den_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']]) 
            den_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']]) 
            den_weights.append(-1*inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP_TauIsPromptLepton']])
            num_pt.append(inputs[channel][DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])
            num_eta.append(inputs[channel][DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])
            num_weight.append(-1*inputs[channel][DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP_TauIsPromptLepton']])            
      else:
        print(f'missing {name_lep}_pt in fields')
      
    den_pt = np.concatenate(den_pt)
    den_eta = np.concatenate(den_eta)
    den_weights = np.concatenate(den_weights)
    num_pt = np.concatenate(num_pt)
    num_eta = np.concatenate(num_eta)
    num_weight = np.concatenate(num_weight)

    den_pt = np.array(den_pt).flatten()
    den_eta = np.array(den_eta).flatten()
    den_weights = np.array(den_weights).flatten()
    num_pt = np.array(num_pt).flatten()
    num_eta = np.array(num_eta).flatten()
    num_weight = np.array(num_weight).flatten()

    output_figures = output_figures + name_fig + '.pdf'
    output_json = output_results + name_fig + '.json'
    output_json_err = output_results + name_fig + '_err.json'
    sfhist, err_sfhist = make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures)

    fake_rate = correctionlib.convert.from_histogram(sfhist)
    fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

    # set overflow bins behavior (default is to raise an error when out of bounds)
    fake_rate.data.flow = "clamp"
    fake_rate_err.data.flow = "clamp"

    save_hist_schemav2(fake_rate, f"fake rate {name_lep}", output_json)
    save_hist_schemav2(fake_rate_err, f"fake rate {name_lep} err", output_json_err)
    return

def compute_FakeRateDataLL(name_lep, inputs ,bins, output_results, output_figures, channel, ControlRegionName, CorrFactor):
    #here channel is Zmu or Ze

    den_pt = []
    den_pt_data = []
    den_eta = []
    den_eta_data = []
    den_weights = []
    num_pt = []
    num_eta = []
    num_weight = []
    name_fig = f'P_fake_TL_{name_lep}'
    
    if f'{name_lep}_pt' in ak.fields(inputs[list(inputs.keys())[0]]):
        for DataSample in inputs.keys():
            if DataSample == 'data':
              cut_region = compute_region_mask(inputs[DataSample], channel, 'data', ControlRegionName)
              #pt_parton = inputs[DataSample][f'{name_lep}_ConeCorrectedPt']-0.15*inputs[DataSample][f'{name_lep}_pt']
              pt_parton = np.concatenate(inputs[DataSample][f'{name_lep}_ConeCorrectedPt'])*CorrFactor[name_lep]
              #p_t_corr = np.concatenate(inputs[DataSample][f'{name_lep}_ConeCorrectedPt']*CorrFactor[name_lep])
              #mask_iso = np.array(inputs[DataSample][f'{name_lep}_pfRelIso03_all'] <= 0.15)
              #p_t_corr[mask_iso] = np.array(inputs[DataSample][f'{name_lep}_pt'])[mask_iso]
              p_t_corr = np.where(cut_region[f'{ControlRegionName}_PassLooseNotTightWP'], pt_parton, inputs[DataSample][f'{name_lep}_pt'])

              den_pt.append(p_t_corr[cut_region[f'{ControlRegionName}_PassLooseWP']])
              den_pt_data.append(p_t_corr[cut_region[f'{ControlRegionName}_PassLooseWP']])
              den_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
              den_eta_data.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP']])
              den_weights.append(inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP']])
              num_pt.append(inputs[DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP']])
              num_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP']])
              num_weight.append(inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP']])  
            else:
              cut_region = compute_region_mask(inputs[DataSample], channel, 'MC', ControlRegionName)
              #pt_parton = inputs[DataSample][f'{name_lep}_ConeCorrectedPt']-0.15*inputs[DataSample][f'{name_lep}_pt']
              pt_parton = np.concatenate(inputs[DataSample][f'{name_lep}_ConeCorrectedPt'])*CorrFactor[name_lep]
              #mask_iso = np.array(inputs[DataSample][f'{name_lep}_pfRelIso03_all'] <= 0.15)
              #p_t_corr[mask_iso] = np.array(inputs[DataSample][f'{name_lep}_pt'])[mask_iso]
              p_t_corr = np.where(cut_region[f'{ControlRegionName}_PassLooseNotTightWP'], pt_parton, inputs[DataSample][f'{name_lep}_pt'])

              den_pt.append(p_t_corr[cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']])
              den_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']]) 
              den_weights.append(-1*inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']])
              num_pt.append(inputs[DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])
              num_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])
              num_weight.append(-1*inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])          
    else:
        print(f'missing {name_lep}_pt in fields')
      
    den_pt = np.concatenate(den_pt)
    den_pt_data = np.concatenate(den_pt_data)
    den_eta = np.concatenate(den_eta)
    den_eta_data = np.concatenate(den_eta_data)
    den_weights = np.concatenate(den_weights)
    num_pt = np.concatenate(num_pt)
    num_eta = np.concatenate(num_eta)
    num_weight = np.concatenate(num_weight)

    den_pt = np.array(den_pt).flatten()
    den_pt_data = np.array(den_pt_data).flatten()
    den_eta = np.array(den_eta).flatten()
    den_eta_data = np.array(den_eta_data).flatten()
    den_weights = np.array(den_weights).flatten()
    num_pt = np.array(num_pt).flatten()
    num_eta = np.array(num_eta).flatten()
    num_weight = np.array(num_weight).flatten()

    output_figures = output_figures + name_fig + '.pdf'
    output_json = output_results + name_fig + '.json'
    output_json_err = output_results + name_fig + '_err.json'

    if isinstance(bins, int):
      pt_bins = np.array(equalObs(den_pt, bins))
      eta_bins = np.array(equalObs(np.abs(den_eta), bins))
    elif len(bins) == 2:
      pt_bins = np.array(bins[0])
      eta_bins = np.array(bins[1])
    else:
      print('pb with nbins')

    sfhist, err_sfhist = make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures, ptlabel = r'$p_T^{corr}$ [GeV]')
    fake_rate = correctionlib.convert.from_histogram(sfhist)
    fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

    # set overflow bins behavior (default is to raise an error when out of bounds)
    fake_rate.data.flow = "clamp"
    fake_rate_err.data.flow = "clamp"

    save_hist_schemav2(fake_rate, f"fake rate {name_lep}", output_json)
    save_hist_schemav2(fake_rate_err, f"fake rate {name_lep} err", output_json_err)
    return

def compute_FakeRateMCLL(name_lep, inputs, pt_bins, eta_bins, output_results, output_figures, channel, ControlRegionName, CorrFactor):
    # here channel is Zmu or Ze

    den_pt = []
    den_eta = []
    den_weights = []
    num_pt = []
    num_eta = []
    num_weight = []
    name_fig = f'P_fake_TL_{name_lep}_MC'

    if f'{name_lep}_pt' in ak.fields(inputs[list(inputs.keys())[0]]):
      for DataSample in inputs.keys():
        if DataSample != 'data':
          cut_region = compute_region_mask(inputs[DataSample], channel, 'MC', ControlRegionName)
          #pt_parton = inputs[DataSample][f'{name_lep}_ConeCorrectedPt']-0.15*inputs[DataSample][f'{name_lep}_pt']
          pt_parton = np.concatenate(inputs[DataSample][f'{name_lep}_ConeCorrectedPt']*CorrFactor[name_lep])
          p_t_corr = np.where(cut_region[f'{ControlRegionName}_PassLooseNotTightWP'], pt_parton, inputs[DataSample][f'{name_lep}_pt'])

          den_pt.append(p_t_corr[cut_region[f'{ControlRegionName}_PassLooseWP']]) 
          den_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP']]) 
          den_weights.append(inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP']])
          num_pt.append(inputs[DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP']])
          num_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP']])
          num_weight.append(inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP']])

          den_pt.append(p_t_corr[cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']]) 
          den_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']]) 
          den_weights.append(-1*inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassLooseWP_{name_lep}IsPromptLepton']])
          num_pt.append(inputs[DataSample][f'{name_lep}_pt'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])
          num_eta.append(inputs[DataSample][f'{name_lep}_eta'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])
          num_weight.append(-1*inputs[DataSample]['genWeight'][cut_region[f'{ControlRegionName}_PassTightWP_{name_lep}IsPromptLepton']])            
    else:
      print(f'missing {name_lep}_pt in fields')
      
    den_pt = np.concatenate(den_pt)
    den_eta = np.concatenate(den_eta)
    den_weights = np.concatenate(den_weights)
    num_pt = np.concatenate(num_pt)
    num_eta = np.concatenate(num_eta)
    num_weight = np.concatenate(num_weight)

    den_pt = np.array(den_pt).flatten()
    den_eta = np.array(den_eta).flatten()
    den_weights = np.array(den_weights).flatten()
    num_pt = np.array(num_pt).flatten()
    num_eta = np.array(num_eta).flatten()
    num_weight = np.array(num_weight).flatten()

    output_figures = output_figures + name_fig + '.pdf'
    output_json = output_results + name_fig + '.json'
    output_json_err = output_results + name_fig + '_err.json'
    sfhist, err_sfhist = make_hist_sf(den_pt, den_eta, den_weights, num_pt, num_eta, num_weight, pt_bins, eta_bins, output_figures)

    fake_rate = correctionlib.convert.from_histogram(sfhist)
    fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

    # set overflow bins behavior (default is to raise an error when out of bounds)
    fake_rate.data.flow = "clamp"
    fake_rate_err.data.flow = "clamp"

    save_hist_schemav2(fake_rate, f"fake rate {name_lep}", output_json)
    save_hist_schemav2(fake_rate_err, f"fake rate {name_lep} err", output_json_err)
    return

lepton = {
  "tem": {
    "lepton1" : {"name": 'Tau', "type":'Tau'},
    "lepton2" : {"name": 'Electron', "type":'Electron'},
    "lepton3" : {"name": 'Muon', "type":'Muon'}
  },
  "ttm": {
    "lepton1" : {"name": 'Tau1', "type":'Tau'},
    "lepton2" : {"name": 'Tau2', "type":'Tau'},
    "lepton3" : {"name": 'Muon', "type":'Muon'}
  },
  "tte": {
    "lepton1" : {"name": 'Tau1', "type":'Tau'},
    "lepton2" : {"name": 'Tau2', "type":'Tau'},
    "lepton3" : {"name": 'Electron', "type":'Electron'}
  },
  "tmm": {
    "lepton1" : {"name": 'Tau', "type": 'Tau'},
    "lepton2" : {"name": 'Muon1', "type":'Muon'},
    "lepton3" : {"name": 'Muon2', "type":'Muon'}
  },
  "tee": {
    "lepton1" : {"name": 'Tau', "type": 'Tau'},
    "lepton2" : {"name": 'Electron1', "type":'Electron'},
    "lepton3" : {"name": 'Electron2', "type":'Electron'}
  }
}

# Define the function to calculate the ratio
def ratio_func(numerator, denominator, prompt, data):
    mask_valid =np.array(prompt.values() / data.values()) < 0.5
    ratio = (numerator.values() - denominator.values()) / denominator.values()
    return np.where(mask_valid, np.array(ratio), np.array(ratio/ratio)*0.3)

'''
def compute_FRweight(FR, W_DY_ttbar, branches, Lepton, cut):
  Lepton_pt = np.array(branches['data'][f'{Lepton["name"]}_pt'][cut]).flatten()
  Lepton_abseta = np.array(np.abs(branches['data'][f'{Lepton["name"]}_eta'][cut])).flatten()
  weights_DY = FR[f'{Lepton["type"]}_DY']['fake_rate'].evaluate(Lepton_pt, Lepton_abseta)
  weights_ttbar = FR[f'{Lepton["type"]}_ttbar']['fake_rate'].evaluate(Lepton_pt, Lepton_abseta)
  weights = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar
  return weights
'''

def compute_FRweight(FR, W_DY_ttbar, branches, Lepton, cut, CorrFactor):
  if Lepton["type"] == 'Tau':
    Lepton_pT = np.array(branches['data'][f'{Lepton["name"]}_pt'][cut]).flatten()
    Lepton_eta = np.array(np.abs(branches['data'][f'{Lepton["name"]}_eta'][cut])).flatten()

    weights_DY = FR['Tau_DY']['fake_rate'].evaluate(Lepton_pT, Lepton_eta)
    weights_ttbar = FR['Tau_ttbar']['fake_rate'].evaluate(Lepton_pT, Lepton_eta)
    weights = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar
  else:
    Lepton_pTParton = np.array(branches['data'][f'{Lepton["name"]}_ConeCorrectedPt'][cut]).flatten()*CorrFactor[Lepton["type"]]
    Lepton_pT = np.array(branches['data'][f'{Lepton["name"]}_pt'][cut]).flatten()
    Lepton_iso = np.array(branches['data'][f'{Lepton["name"]}_pfRelIso03_all'][cut]).flatten()
    Lepton_pTCorr = np.where(Lepton_iso> 0.15, Lepton_pTParton, Lepton_pT)
    Lepton_eta = np.array(np.abs(branches['data'][f'{Lepton["name"]}_eta'][cut])).flatten()

    weights_DY = FR[f'{Lepton["type"]}_DY']['fake_rate'].evaluate(Lepton_pTCorr, Lepton_eta)
    weights_ttbar = FR[f'{Lepton["type"]}_ttbar']['fake_rate'].evaluate(Lepton_pTCorr, Lepton_eta)
    weights = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar
  return weights

def compute_FR(weight):
  return weight/(1-weight)

def apply_FR_method(FR, W_DY_ttbar, branches, channel, Lepton, RegionName, CorrFactor):
  application_region = ['FFF','FFP','FPF','FPP','PFF','PFP','PPF']
  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  #Set values to 0 ...
  sumw_MC = {}
  N_MC = {}
  for XXX in application_region:
    sumw_MC[f'AppRegion{XXX}_TrueLepton'] = 0
    N_MC[f'AppRegion{XXX}_TrueLepton'] = 0
  # ... and fill values
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', RegionName)
  for XXX in application_region:
    sumw_MC[f'AppRegion{XXX}_TrueLepton'] += np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{RegionName}_AppRegion{XXX}_TrueLepton']])
    N_MC[f'AppRegion{XXX}_TrueLepton'] += np.sum(cut_region[f'{RegionName}_AppRegion{XXX}_TrueLepton'])

  Fakes = []
  Weights = []

  cut_region = compute_region_mask(branches['data'], channel, 'data', RegionName)
  for XXX in application_region:
    cut = cut_region[f'{RegionName}_AppRegion{XXX}']
    N_events_XXX = np.sum(cut)

    if XXX == 'FFF':
      weights1 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton1"], cut, CorrFactor)
      weights2 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton2"], cut, CorrFactor)
      weights3 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton3"], cut, CorrFactor)
      weights = compute_FR(weights1)*compute_FR(weights2)*compute_FR(weights3)

    if XXX == 'FFP':
      weights1 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton1"], cut, CorrFactor)
      weights2 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton2"], cut, CorrFactor)
      weights = (-1)*compute_FR(weights1)*compute_FR(weights2)

    if XXX == 'FPF':
      weights1 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton1"], cut, CorrFactor)
      weights3 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton3"], cut, CorrFactor)
      weights = (-1)*compute_FR(weights1)*compute_FR(weights3)

    if XXX == 'PFF':
      weights2 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton2"], cut, CorrFactor)
      weights3 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton3"], cut, CorrFactor)
      weights = (-1)*compute_FR(weights2)*compute_FR(weights3)

    if XXX == 'FPP':
      weights1 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton1"], cut, CorrFactor)
      weights = compute_FR(weights1)

    if XXX == 'PFP':
      weights2 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton2"], cut, CorrFactor)
      weights = compute_FR(weights2)

    if XXX == 'PPF':
      weights3 = compute_FRweight(FR, W_DY_ttbar, branches, lepton[channel]["lepton3"], cut, CorrFactor)
      weights = compute_FR(weights3)

    FakesProp = (N_events_XXX - sumw_MC[f'AppRegion{XXX}_TrueLepton'])/N_events_XXX

    Fakes = np.concatenate( [Fakes,  np.array(branches['data'][f'{Lepton}_pt']).flatten()[cut]] )
    Weights = np.concatenate( [Weights, weights*FakesProp])

  return Fakes, Weights

def apply_FR_methodLL(FR, W_DY_ttbar, branches, channel, Lepton_name, PlotRegion, CorrFactor):

  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  sumw_MC = np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassLooseNotTightWP_{Lepton_name}IsPromptLepton']])
  #print(f' - Prompt {Lepton_name} events in AppRegion = {sumw_MC} \n')

  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  cut = cut_region[f'{PlotRegion}_PassLooseNotTightWP']
  N_events = np.sum(cut)
  #print(f' - N_events in AppRegion= {N_events} \n')

  pt_parton = np.concatenate(branches['data'][f'{Lepton_name}_ConeCorrectedPt'])*CorrFactor[Lepton_name]
  ConeCorrectedPt = np.where(cut, pt_parton, branches['data'][f'{Lepton_name}_pt'])

  weights_DY = FR[f'{Lepton_name}_DY']['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{Lepton_name}_eta'][cut])).flatten())
  weights_ttbar = FR[f'{Lepton_name}_ttbar']['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{Lepton_name}_eta'][cut])).flatten())
  weight = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar

  weights = compute_FR(weight)

  FakesProp = (N_events - sumw_MC)/N_events
  if FakesProp < 0:
      FakesProp = 0

  #print(f' - FakesProp = {round(FakesProp,1)} \n')
  Fakes = compute_ptcorr(branches['data'], Lepton_name, CorrFactor)[cut]
  Weights = weights*FakesProp
  #print(f' - FR weights = {round(np.sum(Weights),1)} \n')

  return Fakes, Weights

def apply_FR_methodTau(FR, W_DY_ttbar, branches, channel, Lepton_name, PlotRegion):

  #Compute sumw/N of TrueLepton in MC for FakesProp computation
  cut_region = compute_region_mask(branches['TrueLepton'], channel, 'MC', PlotRegion)
  sumw_MC = np.sum(branches['TrueLepton']['genWeight'][cut_region[f'{PlotRegion}_PassLooseNotTightWP_{Lepton_name}IsPromptLepton']])
  #print(f' - Prompt {Lepton_name} events in AppRegion = {sumw_MC} \n')

  cut_region = compute_region_mask(branches['data'], channel, 'data', PlotRegion)
  cut = cut_region[f'{PlotRegion}_PassLooseNotTightWP']
  N_events = np.sum(cut)
  #print(f' - N_events in AppRegion= {N_events} \n')

  ConeCorrectedPt = branches['data'][f'{Lepton_name}_pt']

  weights_DY = FR[f'{Lepton_name}_DY']['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{Lepton_name}_eta'][cut])).flatten())
  weights_ttbar = FR[f'{Lepton_name}_ttbar']['fake_rate'].evaluate(np.array(ConeCorrectedPt[cut]).flatten(), np.array(np.abs(branches['data'][f'{Lepton_name}_eta'][cut])).flatten())
  weight = W_DY_ttbar['DY'][0]*weights_DY + W_DY_ttbar['ttbar'][0]*weights_ttbar

  weights = compute_FR(weight)

  FakesProp = (N_events - sumw_MC)/N_events
  if FakesProp < 0:
      FakesProp = 0

  #print(f' - FakesProp = {round(FakesProp,1)} \n')
  Fakes = branches['data'][f'{Lepton_name}_pt'][cut]
  Weights = weights*FakesProp
  #print(f' - FR weights = {round(np.sum(Weights),1)} \n')

  return Fakes, Weights

def compute_corrFactor(branches, channels, leptonname, RegionName, selection, nbins, output_figures):
    font_size = 20
    plt.rcParams.update({
        'axes.labelsize': font_size, 
        'axes.titlesize': font_size, 
        'xtick.labelsize': font_size, 
        'ytick.labelsize': font_size,
        'legend.fontsize': font_size
    })

    # --- Setup ---
    selection = f'{RegionName}_{selection}'

    # --- Binning ---
    bins = np.linspace(0, 0.5, nbins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    pfRelIso03_Data = []
    ConeCorrectedPt_Data = []
    LeptonPt_Data = []

    for channel in channels:
        cut_region_data = compute_region_mask(branches[channel], channel, 'data', RegionName)
        mask_data = cut_region_data[f'{RegionName}_{selection.split("_")[1]}']

        pfRelIso03_Data.append(np.array(branches[channel][mask_data][f'{leptonname}_pfRelIso03_all']))
        ConeCorrectedPt_Data.append(np.array(branches[channel][mask_data][f'{leptonname}_ConeCorrectedPt']))
        LeptonPt_Data.append(np.array(branches[channel][mask_data][f'{leptonname}_pt']))

    pfRelIso03_Data = np.concatenate(pfRelIso03_Data)
    ConeCorrectedPt_Data = np.concatenate(ConeCorrectedPt_Data)
    ConeCorrectedPt_Data = np.concatenate(ConeCorrectedPt_Data)
    LeptonPt_Data = np.concatenate(LeptonPt_Data)
    weights_Data = np.ones_like(ConeCorrectedPt_Data)

    ConeCorrectedPt_Data_original = ConeCorrectedPt_Data.copy()

    # --- Compute Lepton factor ---
    mask_low = (pfRelIso03_Data >= 0.10) & (pfRelIso03_Data < 0.15)
    mask_high = (pfRelIso03_Data >= 0.15) & (pfRelIso03_Data < 0.20)

    avg_low = np.average(LeptonPt_Data[mask_low], weights=weights_Data[mask_low])
    avg_high = np.average(ConeCorrectedPt_Data[mask_high], weights=weights_Data[mask_high])
    Lepton_factor_Data = avg_low / avg_high if avg_high != 0 else 1.0
    print(f"Computed corrected factor for {leptonname} in Data = {Lepton_factor_Data:.3f}")

    ConeCorrectedPt_Data = ConeCorrectedPt_Data * Lepton_factor_Data
    mask_iso_data = pfRelIso03_Data <= 0.15
    ConeCorrectedPt_Data[mask_iso_data] = LeptonPt_Data[mask_iso_data]
    ConeCorrectedPt_Data_original[mask_iso_data] = LeptonPt_Data[mask_iso_data]

    indices_Data = np.digitize(pfRelIso03_Data, bins) - 1
    avg_ConeCorrectedPt_Data = np.zeros(nbins)
    avg_ConeCorrectedPt_Data_orig = np.zeros(nbins)
    for i in range(nbins):
        mask = indices_Data == i
        if np.any(mask):
            avg_ConeCorrectedPt_Data[i] = np.average(ConeCorrectedPt_Data[mask], weights=weights_Data[mask])
            avg_ConeCorrectedPt_Data_orig[i] = np.average(ConeCorrectedPt_Data_original[mask], weights=weights_Data[mask])
        else:
            avg_ConeCorrectedPt_Data[i] = np.nan
            avg_ConeCorrectedPt_Data_orig[i] = np.nan

    # --- Plotting ---
    plt.figure(figsize=(10, 6))
    plt.errorbar(bin_centers, avg_ConeCorrectedPt_Data_orig, fmt='o', capsize=5, label='Data without correction factor', color='grey')
    plt.errorbar(bin_centers, avg_ConeCorrectedPt_Data, fmt='o', capsize=5, label=f'Data with correction factor (f = {Lepton_factor_Data:.3f})', color='black')
    plt.axvline(0.15, linestyle='--', color='gray', alpha=0.6)
    plt.text(0.152, plt.ylim()[1]*0.95, 'cut = 0.15', fontsize=12, color='gray')
    plt.xlabel(f'{leptonname}_pfRelIso03_all', fontsize=18)
    plt.ylabel(r'Average $p_T^{corr}$ [GeV]', fontsize=18)
    plt.title('Avg ConeCorrectedPt vs pfRelIso03')
    plt.grid(True)
    plt.legend()
    plt.legend(fontsize="10", loc ="upper left")
    plt.tight_layout()
    plt.savefig(output_figures + f"correctionFactorsLL_{leptonname}.pdf")
    #plt.show()
    
    return float(Lepton_factor_Data)

def compute_ptcorr(branch, Lepton_name, CorrFactor):
    lepton_pt = branch[f'{Lepton_name}_pt']
    parton_pt = np.concatenate(branch[f'{Lepton_name}_ConeCorrectedPt'])*CorrFactor[Lepton_name]
    pt_corr = np.where(branch[f'{Lepton_name}_pfRelIso03_all'] < 0.15, lepton_pt, parton_pt)
    return np.array(pt_corr).flatten()