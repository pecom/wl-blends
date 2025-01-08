import numpy as np
import healpy as hp
import pandas as pd
import os
import treecorr
from mpi4py import MPI
import argparse


pdir = os.getenv("PSCRATCH")
hdir = os.getenv("HOME")
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

comm.Barrier()
print(f"Saying hi from: {comm.Get_rank()}/{comm.Get_size()}", flush=True)

# Configuration for shear-shear correlation
nbins = 20
min_sep = 6 / 60
max_sep = 150 / 60
ss_config = {
    "nbins": nbins,
    "min_sep": min_sep,
    "max_sep": max_sep,
    "bin_type": "Log",
    "bin_size": np.log(max_sep / min_sep) / nbins,
    "bin_slop": 0.1,
    "flip_g1": False,
    "flip_g2": True,
}

zsb_dir = '/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst_old/zsb/r3/cats'

patch_files=pdir + '/patches/test_patches.fits'
zsb_file = pdir + '/buzzard/zsb_spec3.csv'
# zsb_file = pdir + '/buzzard/zsb_12.csv'


cat = treecorr.Catalog(zsb_file, ra_col=1, dec_col=2, g1_col=3, g2_col=4, k_col=5, w_col=8,
                       ra_units='deg', dec_units='deg', delimiter=',', patch_centers=patch_files)

# cat = treecorr.Catalog(zsb_file, 
#                        ra_col='ra', dec_col='dec',
#                        g1_col='gamma1', g2_col='gamma2',
#                        ra_units='deg', dec_units='deg',
#                        k_col=5, w_col=8,
#                        delimiter=',', flip_g1=False,
#                        flip_g2=True, patch_centers=patch_files
#                        )


if rank==0:
    print("Created catalog...", flush=True)

gg = treecorr.GGCorrelation(
    min_sep=ss_config["min_sep"],
    max_sep=ss_config["max_sep"],
    nbins=ss_config["nbins"],
    bin_type=ss_config["bin_type"],
    bin_slop=ss_config["bin_slop"],
    var_method="jackknife",
    sep_units="degrees",
    verbose=2,
    num_threads=40
)

if rank==0:
    print("Processing!", flush=True)

gg.process(cat, comm=comm, low_mem=False)

comm.Barrier()
if rank==0:
    print("All done! Writing results!", flush=True)

if rank==0:
    gg.write(f'{pdir}/output/zsb_gg4_spec.out')
