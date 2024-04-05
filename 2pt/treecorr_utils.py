import numpy as np
import treecorr
import os

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
    "bin_slop": 0.0,
    "flip_g1": False,
    "flip_g2": True,
}

# Template for patch file path
patch_file = (
    "/scratch/gpfs/en7908/blending/buzzard_v2.0.0_lsst/patch_centers_n{npatch}.fits"
)

# Configuration for covariance calculation
cov_config = {
    "var_method": "jackknife",
    "save_patch_lens_dir": None,
    "save_patch_source_dir": None,
    "save_patch_rand_dir": None,
    "patch_file": patch_file,
    "verbose": 2,
}

# Flag to control whether to remake patches
remake_patches = True


def shear_shear_corr(
    pos1=None,
    pos2=None,
    shear1=None,
    shear2=None,
    ij=None,
    k1=None,
    k2=None,
    w1=None,
    w2=None,
    patch_centers=None,
    num_threads=0,
):
    """
    Calculate shear-shear correlation.

    Parameters
    ----------
    pos1 : array_like or None
        Positions (RA, Dec) of the first shear catalog.
    pos2 : array_like or None
        Positions (RA, Dec) of the second shear catalog. Only needed for cross-correlation.
    shear1 : tuple of array_like
        Shear components (g1, g2) of the first shear catalog.
    shear2 : tuple of array_like
        Shear components (g1, g2) of the second shear catalog. Only needed for cross-correlation.
    ij : tuple of int
        Indices identifying the tomographic bin pair.
    k1 : array_like
        Convergence values for the first shear catalog.
    k2 : array_like
        Convergence values for the second shear catalog. Only needed for cross-correlation.
    w1 : array_like
        Weights for the first shear catalog.
    w2 : array_like
        Weights for the second shear catalog. Only needed for cross-correlation.
    patch_centers : array_like
        Centers of patches for spatial jackknife error estimation.
    num_threads : int
        Number of threads for parallel processing. 0 means using all available threads.

    Returns
    -------
    treecorr.GGCorrelation
        An object containing the shear-shear correlation data.

    Notes
    -----
    This function calculates the auto-correlation (if ij[0] == ij[1]) or
    cross-correlation (if ij[0] != ij[1]) of shear data. The function uses
    TreeCorr for efficient calculation.
    """
    ctype = "ss"  # correlation type, used in `save_patch_source_dir`
    auto = ij[0] == ij[1]  # check if it is an auto-correlation

    # Handle auto-correlation
    if auto:
        ra, dec = pos1  # either pos1 or pos2 works
        g1, g2 = shear1
        k = k1
        w = w1
        # Manage saving and using source patches
        tomoidx = ij[0]
        save_patch_source_dir_compiled = (
            cov_config["save_patch_source_dir"].format(**locals())
            if cov_config["save_patch_source_dir"] != "None"
            else None
        )
        # Check and manage existing patch directories
        if cov_config["save_patch_source_dir"] != "None":
            if remake_patches and os.path.exists(save_patch_source_dir_compiled):
                os.system(f"rm -r {save_patch_source_dir_compiled}")
                print(
                    f"Removed old source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} to make new ones."
                )
            else:
                if os.path.exists(save_patch_source_dir_compiled):
                    print(
                        f"Using old source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx}."
                    )
                else:
                    print(
                        f"Creating new source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} since they do not exist."
                    )
        else:
            print(
                f"{ctype}_{ij[0]}_{ij[1]}  ->  We are not saving source patches and we are not using saved source patches."
            )

        # Prepare catalog for TreeCorr
        cat = treecorr.Catalog(
            g1=g1,
            g2=g2,
            k=k,
            ra=ra,
            dec=dec,
            w=w,
            ra_units="degrees",
            dec_units="degrees",
            flip_g1=ss_config["flip_g1"],
            flip_g2=ss_config["flip_g2"],
            patch_centers=patch_centers,
            save_patch_dir=save_patch_source_dir_compiled,
        )
        del ra, dec, g1, g2, w
    # Handle cross-correlation
    else:
        # Prepare first catalog
        ra1, dec1 = pos1
        ra2, dec2 = pos2
        g1_1st, g2_1st = shear1
        g1_2nd, g2_2nd = shear2

        # Manage saving and using source patches for first catalog
        tomoidx = ij[0]
        save_patch_source_dir_compiled = (
            cov_config["save_patch_source_dir"].format(**locals())
            if cov_config["save_patch_source_dir"] != "None"
            else None
        )
        # Check and manage existing patch directories for first catalog
        if cov_config["save_patch_source_dir"] != "None":
            if remake_patches and os.path.exists(save_patch_source_dir_compiled):
                os.system(f"rm -r {save_patch_source_dir_compiled}")
                print(
                    f"Removed old source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} to make new ones."
                )
            else:
                if os.path.exists(save_patch_source_dir_compiled):
                    print(
                        f"Using old source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx}."
                    )
                else:
                    print(
                        f"Creating new source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} since they do not exist."
                    )
        else:
            print(
                f"{ctype}_{ij[0]}_{ij[1]}  ->  We are not saving source patches and we are not using saved source patches."
            )

        # Prepare second catalog
        tomoidx = ij[1]
        save_patch_source_dir_compiled = (
            cov_config["save_patch_source_dir"].format(**locals())
            if cov_config["save_patch_source_dir"] != "None"
            else None
        )
        # Check and manage existing patch directories for second catalog
        if cov_config["save_patch_source_dir"] != "None":
            if remake_patches and os.path.exists(save_patch_source_dir_compiled):
                os.system(f"rm -r {save_patch_source_dir_compiled}")
                print(
                    f"Removed old source data patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} to make new ones."
                )
            else:
                if os.path.exists(save_patch_source_dir_compiled):
                    print(
                        f"Using old source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx}."
                    )
                else:
                    print(
                        f"Creating new source patches for {ctype}_{ij[0]}_{ij[1]} in tomographic bin {tomoidx} since they do not exist."
                    )
        else:
            print(
                f"{ctype}_{ij[0]}_{ij[1]}  ->  We are not saving source patches and we are not using saved source patches."
            )

        cat1 = treecorr.Catalog(
            g1=g1_1st,
            g2=g2_1st,
            k=k1,
            ra=ra1,
            dec=dec1,
            w=w1,
            ra_units="degrees",
            dec_units="degrees",
            flip_g1=ss_config["flip_g1"],
            flip_g2=ss_config["flip_g2"],
            patch_centers=patch_centers,
            save_patch_dir=save_patch_source_dir_compiled,
        )

        cat2 = treecorr.Catalog(
            g1=g1_2nd,
            g2=g2_2nd,
            k=k2,
            ra=ra2,
            dec=dec2,
            w=w2,
            ra_units="degrees",
            dec_units="degrees",
            flip_g1=ss_config["flip_g1"],
            flip_g2=ss_config["flip_g2"],
            patch_centers=patch_centers,
            save_patch_dir=save_patch_source_dir_compiled,
        )

        del ra1, ra2, dec1, dec2, g1_1st, g2_1st, g1_2nd, g2_2nd

    del pos1, pos2, shear1, shear2, k1, k2, w1, w2

    # Set up and calculate the correlation
    gg = treecorr.GGCorrelation(
        min_sep=ss_config["min_sep"],
        max_sep=ss_config["max_sep"],
        nbins=ss_config["nbins"],
        bin_type=ss_config["bin_type"],
        bin_slop=ss_config["bin_slop"],
        var_method=cov_config["var_method"],
        sep_units="degrees",
    )

    if auto:
        gg.process(cat, num_threads=num_threads)
    else:
        gg.process(cat1, cat2, num_threads=num_threads)

    return gg
