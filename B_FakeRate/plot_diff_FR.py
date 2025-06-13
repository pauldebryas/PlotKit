import json
import matplotlib.pyplot as plt
import numpy as np

'''
Plot diff between 2 sets of fake rates
'''

# Parameters -----------------------------------------------------------------------------------------------------
# First FR to compare
path_to_json_FR1 = 'B_FakeRate/results/LightLeptFFV2/2017/FakeFactorsMuon/DYRegion/' 
json_file_FR1 =    'P_fake_TL_Muon'
# Second FR to compare
path_to_json_FR2 = 'B_FakeRate/results/LightLeptFFV2/2017/FakeFactorsMuon/ttbarRegion/'
json_file_FR2 = 'P_fake_TL_Muon'
# Name of the out figure
output_figures =   'B_FakeRate/figures/LightLeptFFV2/2017/FakeFactorsMuon/diff_FR_DY_ttbar.pdf'
# Label for figure
x_label = r'$p_T^{corr}$ [GeV]'
y_label = r'$\Delta FF (DY-t\bar{t})$'
# ----------------------------------------------------------------------------------------------------------------

# Load the JSON file 
with open(path_to_json_FR1 + f'{json_file_FR1}.json') as f:
    FR1 = json.load(f)

with open(path_to_json_FR2 + f'{json_file_FR2}.json') as f:
    FR2 = json.load(f)

with open(path_to_json_FR1 + f'{json_file_FR1}_err.json') as f:
    FR1_err = json.load(f)

with open(path_to_json_FR2 + f'{json_file_FR2}_err.json') as f:
    FR2_err = json.load(f)

# Extract the bin edges and contents from the JSON files
pt_bins = np.array(FR1['corrections'][0]['data']['edges'][0])  # X-axis bins (pt bins)
eta_bins = np.array(FR1['corrections'][0]['data']['edges'][1])  # Y-axis bins (eta bins)

if pt_bins[0] < 10:
    pt_bins[0] = 10.

# Reshape the content to 2D arrays for both data and MC corrections
sfhist_tl_FR1 = np.array(FR1['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))
sfhist_tl_FR2 = np.array(FR2['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))

# Also reshape the errors for both data and MC corrections
err_sfhist_tl_FR1 = np.array(FR1_err['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))
err_sfhist_tl_FR2 = np.array(FR2_err['corrections'][0]['data']['content']).reshape((len(eta_bins)-1, len(pt_bins)-1))

# Compute the difference in the values between data and MC
diff_sfhist = sfhist_tl_FR1 - sfhist_tl_FR2

# Compute the combined errors
err_diff_sfhist = np.sqrt(err_sfhist_tl_FR1**2 + err_sfhist_tl_FR2**2)

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
c = ax_main.pcolormesh(pt_bins, eta_bins, diff_sfhist.T, cmap='viridis', shading='auto')

# Loop through the bins and place the text in the center of each bin
for i in range(len(eta_bins) - 1):
    for j in range(len(pt_bins) - 1):
        # Calculate the center of the bin
        pt_center = (pt_bins[j] + pt_bins[j+1]) / 2
        eta_center = (eta_bins[i] + eta_bins[i+1]) / 2
        
        # Get the bin value (difference) to display
        bin_value = diff_sfhist.T[i, j]
        
        # Display the value inside the bin
        ax_main.text(pt_center, eta_center, f'{bin_value:.1e}', 
                     ha='center', va='center', color='white', fontsize=8, rotation=90)


# Move the colorbar to the left side
#cbar_ax = fig.add_axes([0., 0.15, 0.0, 0.7])  # Adjust this to position and size the colorbar on the left
#fig.colorbar(c, cax=cbar_ax)
        
# Set labels for the main plot
ax_main.set_xlabel(x_label)
ax_main.set_ylabel(r'$|\eta|$')

# Plot the projection on the X axis (pt)
pt_projection = np.mean(diff_sfhist, axis=1)
pt_projection_err = np.sqrt(np.mean(err_diff_sfhist**2, axis=1))
ax_x.plot(pt_bins[:-1] + np.diff(pt_bins)/2, pt_projection, color='blue')
ax_x.fill_between(pt_bins[:-1] + np.diff(pt_bins)/2, pt_projection - pt_projection_err, pt_projection + pt_projection_err, color='blue', alpha=0.3)  # Error band
ax_x.axhline(0, color='black', linestyle='--')  # Line at y = 0
ax_x.set_ylabel(y_label)
ax_x.set_xscale('log')

# Plot the projection on the Y axis (eta)
eta_projection = np.mean(diff_sfhist, axis=0)
eta_projection_err = np.sqrt(np.mean(err_diff_sfhist**2, axis=0))  # Uncertainty for eta projection
ax_y.plot(eta_projection, eta_bins[:-1] + np.diff(eta_bins)/2, color='blue')
ax_y.fill_betweenx(eta_bins[:-1] + np.diff(eta_bins)/2, eta_projection - eta_projection_err, eta_projection + eta_projection_err, color='blue', alpha=0.3)  # Error band
ax_y.axvline(0, color='black', linestyle='--')  # Line at x = 0
ax_y.set_xlabel(y_label)

# Tight layout for better spacing
plt.tight_layout()

# Save the figure
plt.savefig(output_figures)

