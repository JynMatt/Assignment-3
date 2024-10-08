# Assignment-3
This repo contains the data, python code and figures for the question.

I first create this repo at the github. Next, I do this at my local terminal.

echo "# Assignment-3" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/JynMatt/Assignment-3.git
git push -u origin main

This repo contain Code which has py files for the corresponding question, Figures where the expected output from the py file and laslty the data is the files that extracted to use for the py files.

gaia_Messier_67-result.vot is the file for q2 which is downloaded from gaia with this respective query. 
SELECT gaia.source_id, 
       gaia.ra, 
       gaia.dec, 
       gaia.phot_g_mean_mag, 
       gaia.bp_rp, 
       tmass.designation, 
       tmass.j_m, 
       tmass.h_m, 
       tmass.ks_m
FROM gaiadr3.gaia_source AS gaia
JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xmatch 
  ON gaia.source_id = xmatch.source_id
JOIN gaiadr1.tmass_original_valid AS tmass 
  ON xmatch.original_ext_source_id = tmass.designation
WHERE 1 = CONTAINS(
    POINT('ICRS', gaia.ra, gaia.dec),
    CIRCLE('ICRS', 132.825, 11.8, 1.0)
)
AND gaia.phot_g_mean_mag < 14

nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits is file downloaded from: https://github.com/svenbuder/astr4004_2024_week7/blob/main/ data/nihao_uhd_simulation_g8.26e11_xyz_positions_and_oxygen_ao.fits
