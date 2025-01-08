import numpy as np
import healpy as hp
import pandas as pd
import os
import treecorr

pdir = os.getenv("PSCRATCH")
hdir = os.getenv("HOME")
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

# zsbs = [95, 96, 97]
zsb_dir = '/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst_old/zsb/r3/cats'
zsbs = np.load(hdir + '/rf_data/data/buzzard/zsb_nums.npy')

buzzard = []

for i, zsb in enumerate(zsbs):
    print(f"Done with index {i}/{len(zsbs)}")
    zsb_pd = pd.read_pickle(zsb_dir + f'/zsb.{zsb}_r3.pickle')
    zsb_pd = zsb_pd.filter(['ra', 'dec', 'gamma1', 'gamma2', 'kappa', 'zs', 'ei'])
    zsb_sub = zsb_pd.iloc[::100, :]
    del zsb_pd
    buzzard.append(zsb_sub)

kale = pd.concat(buzzard)

pz_bin4 = (kale['zs'] > .9) * (kale['zs'] < 1.2)
i_snr = 1/(10**(kale['ei']/2.5) - 1)
snr_filt = i_snr > 10
kale_pz = kale[pz_bin4 * snr_filt]

kale_weights = 1/.21**2 * np.ones(len(kale_pz))
patch_files=pdir + '/patches/test_patches.fits'

print("Creating catalog...")
cat = treecorr.Catalog(ra=kale_pz['ra'],
                       dec=kale_pz['dec'],
                       g1=kale_pz['gamma1'],
                       g2=kale_pz['gamma2'],
                       k=kale_pz['kappa'],
                       w=kale_weights,
                       ra_units='degrees',
                       dec_units='degrees',
                       flip_g1=False,
                       flip_g2=True,
                       npatch=200
                       )

print("Writing centers...")
cat.write_patch_centers(patch_files)
