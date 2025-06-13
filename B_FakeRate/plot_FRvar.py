import json
import matplotlib.pyplot as plt
import numpy as np

'''
Plot the variation of FR as function of pt and eta + projection on pt and eta axis
'''
# Parameters -----------------------------------------------------------------------------------------------------
# period = period of data taking 
period = '2017'
# tag = tag used for FR production
tag = 'LightLeptFFV2'
#Lepton name (Tau, Electron or Muon)
LeptonName = 'Muon'
# Region = region where we have computed the FR
Region = 'ttbarRegion'
# Name of the json file where FR are (without .json)
json_file = 'P_fake_TL_Muon'
# Label for figure
x_label = r'$p_T^{corr}$ [GeV]'
# ----------------------------------------------------------------------------------------------------------------

path_to_json =   f'B_FakeRate/results/{tag}/{period}/FakeFactors{LeptonName}/{Region}/'
output_figures = f'B_FakeRate/figures/{tag}/{period}/FakeFactors{LeptonName}/{Region}/'

# Load the JSON file directly
with open(path_to_json + f'{json_file}.json') as f:
    data_tl = json.load(f)

with open(path_to_json + f'{json_file}_err.json') as f:
    data_tl_err = json.load(f)

# Extract the bin edges and contents from the JSON files
pt_bins = np.array(data_tl['corrections'][0]['data']['edges'][0])  # X-axis bins (pt bins)
eta_bins = np.array(data_tl['corrections'][0]['data']['edges'][1])  # Y-axis bins (eta bins)

if pt_bins[0] < 10:
    pt_bins[0] = 10.

# Reshape the content to 2D arrays for both data and MC corrections
sfhist_tl = np.array(data_tl['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))

# Also reshape the errors for both data and MC corrections
err_sfhist_tl = np.array(data_tl_err['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))

# Start plotting
fig = plt.figure(figsize=(8, 8))

# Create the main 2D plot
ax_main = plt.subplot2grid((4, 4), (1, 0), colspan=3, rowspan=3)
ax_main.set_xscale('log')
ax_x = plt.subplot2grid((4, 4), (0, 0), colspan=3, sharex=ax_main)
ax_y = plt.subplot2grid((4, 4), (1, 3), rowspan=3, sharey=ax_main)

# Plot the 2D difference
#c = ax_main.imshow(diff_sfhist.T, origin='lower', cmap='viridis', extent=[pt_bins[0], pt_bins[-1], eta_bins[0], eta_bins[-1]])
#c = ax_main.imshow(diff_sfhist.T, origin='lower', cmap='viridis', extent=[np.log10(pt_bins[0]), pt_bins[-1], eta_bins[0], eta_bins[-1]])
#c = ax_main.imshow(diff_sfhist.T, origin='lower', cmap='viridis', extent=[pt_bins[0], pt_bins[-1], eta_bins[0], eta_bins[-1]], aspect='auto')
c = ax_main.pcolormesh(pt_bins, eta_bins, sfhist_tl.T, cmap='viridis', shading='auto')

# Loop through the bins and place the text in the center of each bin
for i in range(len(eta_bins) - 1):
    for j in range(len(pt_bins) - 1):
        # Calculate the center of the bin
        pt_center = (pt_bins[j] + pt_bins[j+1]) / 2
        eta_center = (eta_bins[i] + eta_bins[i+1]) / 2
        
        # Get the bin value (difference) to display
        bin_value = sfhist_tl.T[i, j]
        
        # Display the value inside the bin
        ax_main.text(pt_center, eta_center, f'{bin_value:.1e}', 
                     ha='center', va='center', color='white', fontsize=8, rotation=90)


# Set labels for the main plot
ax_main.set_xlabel(x_label)
ax_main.set_ylabel(r'$|\eta|$')

# Plot the projection on the X axis (pt)
pt_projection = np.mean(sfhist_tl, axis=1)
pt_projection_err = np.sqrt(np.sum(err_sfhist_tl**2, axis=1))  # Uncertainty for pt projection
ax_x.plot(pt_bins[:-1] + np.diff(pt_bins)/2, pt_projection, color='blue')
ax_x.fill_between(pt_bins[:-1] + np.diff(pt_bins)/2, pt_projection - pt_projection_err, pt_projection + pt_projection_err, color='blue', alpha=0.3)  # Error band
ax_x.axhline(np.mean(pt_projection), color='black', linestyle='--')  # Line at y = mean(sf)
ax_x.set_ylabel(r'Proj. FR (mean)')
ax_x.set_xscale('log')

# Plot the projection on the Y axis (eta)
eta_projection = np.mean(sfhist_tl, axis=0)
eta_projection_err = np.sqrt(np.sum(err_sfhist_tl**2, axis=0))  # Uncertainty for eta projection
ax_y.plot(eta_projection, eta_bins[:-1] + np.diff(eta_bins)/2, color='blue')
ax_y.fill_betweenx(eta_bins[:-1] + np.diff(eta_bins)/2, eta_projection - eta_projection_err, eta_projection + eta_projection_err, color='blue', alpha=0.3)  # Error band
ax_y.axvline(np.mean(eta_projection), color='black', linestyle='--')  # Line at y = mean(sf)
ax_y.set_xlabel(r'Proj. FR (mean)')

# Tight layout for better spacing
plt.tight_layout()

# Save the figure
plt.savefig(output_figures + f'FR_variation.pdf')

