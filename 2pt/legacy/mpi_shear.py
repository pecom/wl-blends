import pandas as pd
import numpy as np
import treecorr
from mpi4py import MPI
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()

verbose = True

# Configuration for shear-shear correlation
nbins = 20
min_sep = 6 / 60
max_sep = 150 / 60
bin_slop = 0.0
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
	'dec_units': 'degrees',
	'first_row':2
}

zbin = 4

# Template for patch file path
# patch_file = (
#     "/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst/patch_centers_n{npatch}.fits"
# )

# patch_file = "/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst/patch_centers_n200.fits"
patch_file = "/scratch/gpfs/pa2937/patches/patch_centers_pz4.fits"
# Pre-computed 200 patches

# Configuration for covariance calculation
cov_config = {
    "var_method": "jackknife",
    "save_patch_lens_dir": None,
    "save_patch_source_dir": "/scratch/gpfs/pa2937/tmp/",
    "save_patch_rand_dir": None,
    "patch_file": patch_file,
    "verbose": 2,
}
#    "save_patch_source_dir": "/scratch/gpfs/pa2937/patches/",
if verbose: print("Loading file...")

bz_file = f'/scratch/gpfs/pa2937/pz/pz_bin{zbin}.csv'

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
	save_patch_dir=cov_config["save_patch_source_dir"],
	patch_centers = patch_file,
	file_type='ASCII',
	delimiter=',',
)

if verbose: print("Catalog constructed ..")

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
gg.process(cat, comm=comm, low_mem=True)

comm.Barrier()
if verbose: print("Writing output")
if rank==0:
	gg.write(f'/home/pa2937/rf_data/output/buzzard/gg_{zbin}.out')
