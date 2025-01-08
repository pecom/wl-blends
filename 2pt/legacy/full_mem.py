import pandas as pd
import numpy as np
import treecorr

import matplotlib.pyplot as plt

verbose = True

# Configuration for shear-shear correlation
nbins = 20
min_sep = 6 / 60
max_sep = 150 / 60
bin_slop = 0.1
ss_config = {
    "nbins": nbins,
    "min_sep": min_sep,
    "max_sep": max_sep,
    "bin_type": "Log",
    "bin_size": np.log(max_sep / min_sep) / nbins,
    "bin_slop": bin_slop,
    "flip_g1": False,
    "flip_g2": True,
}

cat_config = {
	'ra_col': 1,
	'dec_col': 2,
	'g1_col': 9,
	'g2_col': 10,
	'k_col': 11,
	'ra_units': 'degrees',
	'dec_units': 'degrees'
}

zbin = 4

# Template for patch file path
patch_file = (
    "/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst/patch_centers_n{npatch}.fits"
)

patch_file = "/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst/patch_centers_n200.fits"

# Configuration for covariance calculation
cov_config = {
    "var_method": "jackknife",
    "save_patch_lens_dir": None,
    "save_patch_source_dir": "/scratch/gpfs/pa2937/patches/",
    "save_patch_rand_dir": None,
    "patch_file": patch_file,
    "verbose": 2,
}
if verbose: print("Loading file...")

# bz = pd.read_csv('/home/pa2937/rf_data/data/buzzard/subset1.csv')
# bz = pd.read_csv('/home/pa2937/rf_data/data/buzzard/subset2.csv')
# bz = pd.read_csv(f'/scratch/gpfs/pa2937/pz/pz_bin{zbin}.csv')
if verbose: print("File loaded...")
# if verbose: print(f"We have {len(bz)} rows")

# Flag to control whether to remake patches
remake_patches = False
ctype = "ss"  # correlation type, used in `save_patch_source_dir`

# ra = bz['ra']
# dec = bz['dec']
# g1 = bz['gamma1']
# g2 = bz['gamma2']
# k = bz['kappa']
# w = np.ones_like(k)

cat = treecorr.Catalog(
	bz_file,
	cat_config,
	flip_g1=ss_config["flip_g1"],
	flip_g2=ss_config["flip_g2"],
	npatch = 200,
	save_patch_dir=cov_config["save_patch_source_dir"],
	patch_centers = patch_file,
	file_type='ASCII',
	delimiter=',',
)

if verbose: print("Catalog constructed ..")
# del ra, dec, g1, g2, w

# Set up and calculate the correlation
gg = treecorr.GGCorrelation(
	min_sep=ss_config["min_sep"],
	max_sep=ss_config["max_sep"],
	nbins=ss_config["nbins"],
	bin_type=ss_config["bin_type"],
	bin_slop=ss_config["bin_slop"],
	var_method='jackknife',
	sep_units="degrees",
)

if verbose: print("Calculating shear-shear autocorrelation...")
gg.process(cat, low_mem=True)

if verbose: print("Writing output")
gg.write(f'/home/pa2937/rf_data/output/buzzard/gg_{zbin}_test.out')

if verbose: print("Done! Creating plots")

#  r = np.exp(gg.meanlogr)
#  xip = gg.xip
#  xim = gg.xim
#  sig = np.sqrt(gg.varxip)
#  
#  plt.plot(r, xip, color='blue')
#  plt.plot(r, -xip, color='blue', ls=':')
#  plt.errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='blue', lw=0.1, ls='')
#  plt.errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='blue', lw=0.1, ls='')
#  lp = plt.errorbar(-r, xip, yerr=sig, color='blue')
#  
#  plt.plot(r, xim, color='green')
#  plt.plot(r, -xim, color='green', ls=':')
#  plt.errorbar(r[xim>0], xim[xim>0], yerr=sig[xim>0], color='green', lw=0.1, ls='')
#  plt.errorbar(r[xim<0], -xim[xim<0], yerr=sig[xim<0], color='green', lw=0.1, ls='')
#  lm = plt.errorbar(-r, xim, yerr=sig, color='green')
#  
#  plt.xscale('log')
#  plt.yscale('log', nonpositive='clip')
#  plt.xlabel(r'$\theta$ (arcmin)')
#  
#  plt.legend([lp, lm], [r'$\xi_+(\theta)$', r'$\xi_-(\theta)$'])
#  # plt.xlim( [1,200] )
#  plt.ylabel(r'$\xi_{+,-}$')
#  plt.savefig('/home/pa2937/rf_data/output/graphs/xi_test.png')
#  #plt.show()

