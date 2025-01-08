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

def zsb_filt(pix):
    plen = len(pix)

    i_snr = 1/(10**(pix['ei']/2.5) - 1)
    snr_filt = i_snr > 10
    print(f"SNR filter: {np.sum(snr_filt)/plen}")

    e1_noisy = pix['e1_lensed_convolved_noisy']
    e2_noisy = pix['e2_lensed_convolved_noisy']
    enoisy = np.sqrt(e1_noisy**2 + e2_noisy**2) 
    e_filt = enoisy < .8
    print(f"Ellips. filter: {np.sum(e_filt)/plen}")
    
    e1 = pix['e1_lensed_only']
    e2 = pix['e2_lensed_only']
    etrue = np.sqrt(e1**2 + e2**2)

    deltae = np.abs(enoisy - etrue)
    pix_weight = 1/(.21**2 + deltae**2)
    pix['weight'] = pix_weight
    deltae_filt = deltae < .25
    print(f"Shape noise filter: {np.sum(deltae_filt)/plen}")
    
    total_filt = snr_filt * e_filt * deltae_filt 
    print(f"Total filter: {np.sum(total_filt)/plen}")
    pix_clean = pix[total_filt]
    
    return pix_clean

def ez_filt(pix, blend=False):
    plen = len(pix)
    
    deltae = pix['delta_e']
    pix_weight = 1/(.21**2 + deltae**2)
    pix['weight'] = pix_weight

    weaklen_filt = (pix['bad_shape'] == 0)
    if blend:
        weaklen_filt = (pix['bad_shape'] == 0) * (pix['bflag'] < 2)
    pix_clean = pix[weaklen_filt]
    return pix_clean
    

# zsbs = [95, 96, 97]
zsb_dir = '/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst_old/zsb/r3/cats'
# zsbs = np.load(hdir + '/rf_data/data/buzzard/zsb_nums.npy')
zsbs = np.load(hdir + '/rf_data/data/buzzard/zsb_deconvolved.npy')


buzzard = []

for i, zsb in enumerate(zsbs):
    print(f"Done with index {i}/{len(zsbs)}", flush=True)
    zsb_pd = pd.read_pickle(zsb_dir + f'/zsb.{zsb}_r3_deconvolved.pickle')
    clean_zsb = ez_filt(zsb_pd, True)
    zsb_good = clean_zsb.filter(['ra', 'dec', 'gamma1', 'gamma2', 'kappa', 'zs', 'ei', 'weight'])
    buzzard.append(zsb_good)

kale = pd.concat(buzzard)
print("Made kale...")
kale.to_csv(f'{pdir}/buzzard/zsb_full.csv', index=False)

zbins = [.3,.5,.7,.9,1.2]

for kz_ndx in range(4):
    print(f"Creating kz {kz_ndx}")
    spec_filt = (kale['zs'] > zbins[kz_ndx]) * (kale['zs'] < zbins[kz_ndx+1])
    kale_pz = kale[spec_filt]
    kale_pz.to_csv(f'{pdir}/buzzard/zsb_spec{kz_ndx}_noblends.csv', index=False)

