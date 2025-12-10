# Download sample RMM & RMM-C46 files for tt, X->HH, X->SH samples 
Download from https://cernbox.cern.ch/s/XRoG379lvSCfFVE 

# Compare performance with supervised MLP training:
python mlp_compare_rmm_supervised_all.py

# Compare performance with unsupervised AE training:
python ae_compare_rmm_all.py

# Compare performance with unsupervised VAE training:
python vae_compare_rmm_all.py

# Display RMM-C46 plots (heatmap):
python ../visualize/display_RMM_C46.py



