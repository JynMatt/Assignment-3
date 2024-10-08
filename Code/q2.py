from astropy.io.votable import parse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Parse the VOTable file
votable = parse("A3-result-2.vot")
table = votable.get_first_table()

# Convert the data into a Pandas DataFrame 
df = table.to_table().to_pandas()

# Quality cut 1: Filter out stars with bad 2MASS photometry (ph_qual not 'AAA')
df_good_photometry = df[df['ph_qual'] == 'AAA']

# Quality cut 2: Filter out stars with non-positive parallaxes (parallax <= 0)
df_final = df_good_photometry[df_good_photometry['parallax'] > 0]
df_final_count = df[df['parallax'] < 0]
# Output how many stars remain after applying both quality cuts
remaining_stars_count = df_final.shape[0]

# Print the result
print(f"Number of stars on intial query: {df.shape[0]}")
print(f"Number of stars with bad 2MASS photometry: {df.shape[0]-df_good_photometry.shape[0]}")
print(f"Number of stars with with negative (or non-positive) parallaxes in the Gaia data: {df_final_count.shape[0]}")
print(f"Number of stars remaining after cuts: {remaining_stars_count}")

# Plotting 
df_final = df_final.copy()

# Calculate distance in parsecs
df_final.loc[:, 'distance_pc'] = 1000.0 / df_final['parallax']

# Calculate absolute G magnitude
df_final.loc[:, 'G_abs'] = df_final['phot_g_mean_mag'] - 5 * np.log10(df_final['distance_pc'] / 10)

plt.figure(figsize=(12, 6))
# Panel (a) - Color-Magnitude Diagram (CMD)
# Subplot 1: Gaia BP-RP vs. Absolute G magnitude
plt.subplot(1, 2, 1)
plt.scatter(df_final['bp_rp'], df_final['G_abs'], color='blue', s=10)
plt.gca().invert_yaxis()  # Invert y-axis since magnitudes are lower for brighter stars
plt.xlabel("Gaia BP-RP Color")
plt.ylabel("Absolute G Magnitude")
plt.title("Gaia BP-RP vs Absolute G Magnitude (CMD)")

# Panel (b) - 2MASS J-Ks vs. apparent K magnitude diagram
# Subplot 2: 2MASS J-Ks vs Apparent Ks magnitude
df_final.loc[:, 'J_minus_Ks'] = df_final['j_m'] - df_final['ks_m']  # J-Ks color
plt.subplot(1, 2, 2)
plt.scatter(df_final['J_minus_Ks'], df_final['ks_m'], color='green', s=10)
plt.gca().invert_yaxis() 
plt.xlabel("2MASS J-Ks Color")
plt.ylabel("Apparent Ks Magnitude")
plt.title("2MASS J-Ks vs Apparent Ks Magnitude")

# Show the figure
plt.tight_layout()
plt.savefig('Figures/cmds_M67.png', dpi=200)
print("cmds M67.png has been save to your local machine")
plt.show()
