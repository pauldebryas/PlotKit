import hist
import os
import correctionlib.convert
import numpy as np
import awkward as ak
from hist_makers.helpers import equalObs
from hist_makers.regions.FakeRate_regions import compute_region_mask
import matplotlib.pyplot as plt

def make_hist_sf(den_tau_pt, den_tau_eta, den_weights, num_tau_pt, num_tau_eta, num_weight, n_bins, name_fig = None):
  pt_bins = np.array( equalObs( den_tau_pt, n_bins) )
  eta_bins = np.array( equalObs( np.abs( den_tau_eta) , n_bins) )
  dists = (
      hist.Hist.new
      .StrCat(["num", "den"], name="selection", growth=True)
      .Variable(pt_bins, name="pt")
      .Variable(eta_bins, name="abs_eta")
      .Weight()
      .fill(
          selection="num",
          pt=num_tau_pt,
          abs_eta=np.abs(num_tau_eta),
          weight=num_weight
      )
      .fill(
          selection="den",
          pt=den_tau_pt,
          abs_eta=np.abs(den_tau_eta),
          weight=den_weights
      )
  )
  num = np.array(dists["num", :, :].values())
  den = np.array(dists["den", :, :].values())
  err_num = np.where( num > 0, np.sqrt(num), 0)
  err_den =  np.where( den > 0, np.sqrt(den), 0)

  #We will set it to 0 anywhere we run out of statistics for the correction, to avoid divide by zero issues.
  sf = np.where((den > 0), num / den, 0)
  err_sf = np.where((num > 0) & (den > 0), sf*np.sqrt( (err_num/num)**2 + (err_den/den)**2 ), 0) 

  sfhist = hist.Hist(*dists.axes[1:], data=sf)
  err_sfhist = hist.Hist(*dists.axes[1:], data=err_sf)

  # without a name, the resulting object will fail validation
  sfhist.name = "fake_rate"
  err_sfhist.name = 'fake_rate_err'
  sfhist.label = "FR"
  err_sfhist.label = 'FR_err'

  #plotting fake rate
  if name_fig != None:
    fig, ax = plt.subplots(figsize=(10, 10))
    sfhist.plot2d(ax=ax)

    for i in range(len(eta_bins)-1):
        for j in range(len(pt_bins)-1):
            str = f"{round(sfhist.values().T[i,j],2)} \n $\pm$ {round(err_sfhist.values().T[i,j],2)}" 
            ax.text(pt_bins[j]+((pt_bins[j+1]-pt_bins[j])/2),eta_bins[i]+((eta_bins[i+1]-eta_bins[i])/2), str, color="w", ha="center", va="center", fontweight="normal")

    ax.set(xlabel=r'$p_t$ [GeV]', ylabel = r'$|\eta|$')
    ax.set_xscale('log')
    plt.savefig(f'{os.getenv("RUN_PATH")}/FakeRate/figures/{name_fig}.pdf')

  return sfhist, err_sfhist

def save_hist_schemav2(hist, desc, filename):
  cset = correctionlib.schemav2.CorrectionSet(
      schema_version=2,
      description=desc,
      corrections=[
          hist
      ],
  )
  with open(f'{os.getenv("RUN_PATH")}/FakeRate/results/{filename}.json', "w") as fout:
      fout.write(cset.json(exclude_unset=True))
  return

def compute_FakeRate(inputs, MCsample, channel ,n_bins):
  den_tau_pt = []
  den_tau_eta = []
  den_weights = []
  num_tau_pt = []
  num_tau_eta = []
  num_weight = []

  if channel != 'all':
    name_fig = f'P_fake_TL_{MCsample}_{channel}'
    cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC')
    den_tau_pt.append(inputs[channel][MCsample]['Tau_Jet_pt'][cut_region['ControlRegion_fake']]) 
    den_tau_eta.append(inputs[channel][MCsample]['Tau_Jet_eta'][cut_region['ControlRegion_fake']]) 
    den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region['ControlRegion_fake']])
    num_tau_pt.append(inputs[channel][MCsample]['Tau_Jet_pt'][cut_region['ControlRegion_pass_fake']])
    num_tau_eta.append(inputs[channel][MCsample]['Tau_Jet_eta'][cut_region['ControlRegion_pass_fake']])
    num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region['ControlRegion_pass_fake']])

  else:
    name_fig = f'P_fake_TL_{MCsample}'
    for channel in ['tem','tmm','tee']:
      cut_region = compute_region_mask(inputs[channel][MCsample], channel, 'MC')
      den_tau_pt.append(inputs[channel][MCsample]['Tau_Jet_pt'][cut_region['ControlRegion_fake']])
      den_tau_eta.append(inputs[channel][MCsample]['Tau_Jet_eta'][cut_region['ControlRegion_fake']])
      den_weights.append(inputs[channel][MCsample]['genWeight'][cut_region['ControlRegion_fake']])
      num_tau_pt.append(inputs[channel][MCsample]['Tau_Jet_pt'][cut_region['ControlRegion_pass_fake']])
      num_tau_eta.append(inputs[channel][MCsample]['Tau_Jet_eta'][cut_region['ControlRegion_pass_fake']])
      num_weight.append(inputs[channel][MCsample]['genWeight'][cut_region['ControlRegion_pass_fake']])

  den_tau_pt = np.concatenate(den_tau_pt)
  den_tau_eta = np.concatenate(den_tau_eta)
  den_weights = np.concatenate(den_weights)
  num_tau_pt = np.concatenate(num_tau_pt)
  num_tau_eta = np.concatenate(num_tau_eta)
  num_weight = np.concatenate(num_weight)

  den_tau_pt = np.array(den_tau_pt).flatten()
  den_tau_eta = np.array(den_tau_eta).flatten()
  num_tau_pt = np.array(num_tau_pt).flatten()
  num_tau_eta = np.array(num_tau_eta).flatten()

  sfhist, err_sfhist = make_hist_sf(den_tau_pt, den_tau_eta, den_weights, num_tau_pt, num_tau_eta, num_weight, n_bins, name_fig)

  fake_rate = correctionlib.convert.from_histogram(sfhist)
  fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

  # set overflow bins behavior (default is to raise an error when out of bounds)
  fake_rate.data.flow = "clamp"
  #fake_rate_err.flow = "clamp"

  save_hist_schemav2(fake_rate, f"fake rate", name_fig)
  save_hist_schemav2(fake_rate_err, f"fake rate err", 'err_' + name_fig)
  return

def compute_FakePropInSideband(inputs, channel, n_bins):
  cut_region_data = compute_region_mask(inputs[channel]['data'], channel, 'data')

  # Combined all MC backgrounds
  list = []
  for key in inputs[channel].keys():
     if key not in ['data','HNL']:
        list.append(inputs[channel][key])
  allMCbackground = ak.concatenate(list)

  cut_region_allMCbackground = compute_region_mask(allMCbackground, channel, 'MC')

  # Fake in sideband: data - True (from MC)
  num_tau_pt = np.concatenate([inputs[channel]["data"]['Tau_pt'][cut_region_data['ControlRegion_fail']], allMCbackground['Tau_pt'][cut_region_allMCbackground['ControlRegion_fail_true']]])
  num_tau_eta = np.abs(np.concatenate([inputs[channel]["data"]['Tau_eta'][cut_region_data['ControlRegion_fail']], allMCbackground['Tau_eta'][cut_region_allMCbackground['ControlRegion_fail_true']]]))
  num_weight = np.concatenate([1.*inputs[channel]["data"]['genWeight'][cut_region_data['ControlRegion_fail']], -1.*allMCbackground['genWeight'][cut_region_allMCbackground['ControlRegion_fail_true']]])

  den_tau_pt = inputs[channel]["data"]['Tau_pt'][cut_region_data['ControlRegion_fail']]
  den_tau_eta = inputs[channel]["data"]['Tau_eta'][cut_region_data['ControlRegion_fail']]
  den_weights = 1.*inputs[channel]["data"]['genWeight'][cut_region_data['ControlRegion_fail']]

  sfhist, err_sfhist = make_hist_sf(den_tau_pt, den_tau_eta, den_weights, num_tau_pt, num_tau_eta, num_weight, n_bins, f'P_FakeInSideband_{channel}')

  fake_rate = correctionlib.convert.from_histogram(sfhist)
  fake_rate_err = correctionlib.convert.from_histogram(err_sfhist)

  # set overflow bins behavior (default is to raise an error when out of bounds)
  fake_rate.data.flow = "clamp"
  #fake_rate_err.flow = "clamp"

  save_hist_schemav2(fake_rate, f"fake in sideband {channel}", f'P_FakeInSideband_{channel}')
  save_hist_schemav2(fake_rate_err, f"fake in sideband err {channel}", f'P_FakeInSideband_err_{channel}')
  return 
