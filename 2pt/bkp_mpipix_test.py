import numpy as np
import healpy as hp
import pandas as pd
import os
import treecorr
from mpi4py import MPI

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
# zsbs = np.load(hdir + '/rf_data/data/buzzard/zsb_nums.npy')

# buzzard = []
# 
# for i, zsb in enumerate(zsbs):
#     if i%2==0:
#         print(f"Done with index {i}/{len(zsbs)}")
#     zsb_pd = pd.read_pickle(zsb_dir + f'/zsb.{zsb}_r3.pickle')
#     zsb_pd = zsb_pd.filter(['ra', 'dec', 'gamma1', 'gamma2', 'kappa', 'zs', 'ei'])
#     buzzard.append(zsb_pd)
# 
# kale = pd.concat(buzzard)

# if rank==0:
#     print("Concated pixels...")
# 
# pz_bin4 = (kale['zs'] > .9) * (kale['zs'] < 1.2)
# i_snr = 1/(10**(kale['ei']/2.5) - 1)
# snr_filt = i_snr > 10
# kale_pz = kale[pz_bin4 * snr_filt]
# 
# kale_weights = 1/.21**2 * np.ones(len(kale_pz))

patch_files=pdir + '/patches/test_patches.fits'
# zsb_file = pdir + '/buzzard/zsb_12.csv'
zsb_file = pdir + '/buzzard/zsb_zs4.csv'


cat = treecorr.Catalog(zsb_file,
                       ra_col='ra', dec_col='dec',
                       g1_col='gamma1', g2_col='gamma2',
                       ra_units='deg', dec_units='deg',
                       k_col=5, w_col=8,
                       delimiter=',', flip_g1=False,
                       flip_g2=True, patch_centers=patch_files
                       )

print(cat.g2)

# cat = treecorr.Catalog(ra=kale_pz['ra'],
#                        dec=kale_pz['dec'],
#                        g1=kale_pz['gamma1'],
#                        g2=kale_pz['gamma2'],
#                        k=kale_pz['kappa'],
#                        w=kale_weights,
#                        ra_units='degrees',
#                        dec_units='degrees',
#                        flip_g1=False,
#                        flip_g2=True,
#                        patch_centers=patch_files,
#                        save_patch_dir=pdir+'/patches/')

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
)

if rank==0:
    print("Processing!", flush=True)

gg.process(cat, comm=comm, low_mem=False)

comm.Barrier()
if rank==0:
    print("All done! Writing results!", flush=True)

if rank==0:
    gg.write(f'{pdir}/output/zsb_gg4.out')
