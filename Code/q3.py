from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy import stats
import matplotlib.colors as mcolors
from scipy import optimize
from sklearn.metrics import mean_squared_error, r2_score

# Open the FITS file
fits_file = "nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits"
hdul = fits.open(fits_file)
data = hdul[1].data
x = data['x']
y = data['y']
z = data['z']
RGal = np.sqrt(x**2 + y**2 + z**2)  # Galactocentric radius
AO = data['A_O'] 

# Set up the figure and axes
fig, axs = plt.subplots(1, 2, figsize=(18, 6))  # Two panels side by side

# Panel 1: Logarithmic density plot
ax = axs[0]
h = ax.hist2d(RGal, AO, bins=100, cmap='plasma', norm=LogNorm())
fig.colorbar(h[3], ax=ax, label='Density')

# Linear fit
slope, intercept, r_value, p_value, std_err = stats.linregress(RGal, AO)
ax.plot(RGal, intercept + slope * RGal, color='white', label=f'Fit: y={slope:.2f}x + {intercept:.2f}')
ax.set_xlabel('RGal (kpc)')
ax.set_ylabel('A(O)')
ax.set_title('Logarithmic Density Plot')
ax.legend()

# Panel 2: Residuals plot
ax = axs[1]
residuals = AO - (intercept + slope * RGal)
h_res = ax.hist2d(RGal, residuals, bins=100, cmap='plasma', norm=LogNorm())
fig.colorbar(h_res[3], ax=ax, label='Density')
ax.axhline(0, color='r', linestyle='--')
ax.set_xlabel('RGal (kpc)')
ax.set_ylabel('Residuals ΔA(O)')
ax.set_title('Residuals of the Fit')

# Adjust layout
plt.tight_layout()

# Save the figure to a file before displaying it
plt.savefig('Figures/log_density_and_residuals_fit.png', dpi=200)

# Display the plot
plt.show()

# Define a linear function
def linear_func(x, m, b):
    return m * x + b

# Use curve_fit from scipy.optimize to fit a linear model
params, covariance = optimize.curve_fit(linear_func, RGal, AO)

# Extracting slope and intercept and calculating uncertainties from covariance matrix
slope, intercept = params
slope_err, intercept_err = np.sqrt(np.diag(covariance))

# Print the results
print(f"Slope: {slope} ± {slope_err}")
print(f"Intercept: {intercept} ± {intercept_err}")

# Predicted values
AO_pred = linear_func(RGal, slope, intercept)

# Calculate residuals
residuals = AO - AO_pred

# Compute overall RMSE and R²
rmse_total = np.sqrt(mean_squared_error(AO, AO_pred))
r2_total = r2_score(AO, AO_pred)

# Print the overall RMSE and R²
print(f"Overall RMSE: {rmse_total}")
print(f"Overall R²: {r2_total}")

# Split the data into regions for further analysis (e.g., inner and outer regions)
# Example: Inner region RGal < 10 kpc, Outer region RGal > 10 kpc
inner_region = RGal < 10
outer_region = RGal >= 10

# Compute RMSE for inner and outer regions
rmse_inner = np.sqrt(mean_squared_error(AO[inner_region], AO_pred[inner_region]))
rmse_outer = np.sqrt(mean_squared_error(AO[outer_region], AO_pred[outer_region]))

# Print RMSE for different regions
print(f"RMSE (RGal < 10 kpc): {rmse_inner}")
print(f"RMSE (RGal >= 10 kpc): {rmse_outer}")

# Fit the data using curve_fit
params, _ = optimize.curve_fit(linear_func, RGal, AO)
slope, intercept = params

# Fitted A(O)
AO_fitted = linear_func(RGal, slope, intercept)

# Residuals
AO_residuals = AO - AO_fitted

# Define the 2D bins for x and y
x_bins = np.linspace(min(x), max(x), 80)
y_bins = np.linspace(min(y), max(y), 80)

# (a) 2D-histogram of the median simulated A(O)
simulated_AO_hist, xedges, yedges, binnumber = stats.binned_statistic_2d(x, y, AO, statistic='median', bins=[x_bins, y_bins])

# (b) 2D-histogram of the median fitted A(O)
fitted_AO_hist, _, _, _ = stats.binned_statistic_2d(x, y, AO_fitted, statistic='median', bins=[x_bins, y_bins])

# (c) 2D-histogram of the median residuals ∆A(O)
residuals_AO_hist, _, _, _ = stats.binned_statistic_2d(x, y, AO_residuals, statistic='median', bins=[x_bins, y_bins])

# Plot the 3-panel figure
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# Plot (a) 2D-histogram of the median simulated A(O)
im1 = axs[0].imshow(simulated_AO_hist.T, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', cmap='plasma')
axs[0].set_title('Simulated A(O)')
axs[0].set_xlabel('x (kpc)')
axs[0].set_ylabel('y (kpc)')
fig.colorbar(im1, ax=axs[0])

# Plot (b) 2D-histogram of the median fitted A(O)
im2 = axs[1].imshow(fitted_AO_hist.T, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', cmap='plasma')
axs[1].set_title('Fitted A(O)')
axs[1].set_xlabel('x (kpc)')
axs[1].set_ylabel('y (kpc)')
fig.colorbar(im2, ax=axs[1])

# Plot (c) 2D-histogram of the median residuals ∆A(O)
im3 = axs[2].imshow(residuals_AO_hist.T, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', cmap='coolwarm')
axs[2].set_title('Residuals ∆A(O)')
axs[2].set_xlabel('x (kpc)')
axs[2].set_ylabel('y (kpc)')
fig.colorbar(im3, ax=axs[2])

plt.savefig('Figures/3 panel figure for the x vs y plane of 2Dhistogram of the median A(O)', dpi=200)

plt.show()