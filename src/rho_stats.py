"""
Functions to calculate rho statistics
PSFs and stars should be Maker-type objects
"""

# Global imports
import numpy as np
import treecorr
import os
import pdb

# Local imports
from src.rho_stats_plotter import plot_rho_stats

def run_rho_stats(
    psfs, stars, psf_type, rho_params=None, vb=False, outdir='./'):
    """
    Method to obtain rho-statistics for current PSF fit & save plots
    Requires StarMaker to be provided through 'stars' parameter
    """

    # If rho params not set, use these defaults:
    if not rho_params:
        rho_params={'min_sep':200,'max_sep':12000,'nbins':12}

    # First do calculations
    rho1, rho2, rho3, rho4, rho5 = _run_rho_stats(
        psfs, stars, rho_params, psf_type, outdir, vb
    )

    # Set output name
    outname = os.path.join(
        outdir,'_'.join([psf_type,'rho_stats.pdf'])
    )

    # Make plots
    plot_rho_stats(
        rho1, rho2, rho3, rho4, rho5,
        pixel_scale=stars.pixel_scale, outname=outname
    )

    print("Finished rho stat computation & plotting")

def _run_rho_stats(psfs, stars, rho_params, psf_type, outdir, vb):

    min_sep = rho_params['min_sep']
    max_sep = rho_params['max_sep']
    nbins = rho_params['nbins']

    # Define quantities to be used
    wg = (psfs.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)
    if len(wg.nonzero()[0])<1:
        print("too many stars failed, exiting")

    x = stars.x; y = stars.y
    star_g1 = stars.hsm_g1[wg]
    star_g2 = stars.hsm_g2[wg]
    psf_g1 = psfs.hsm_g1[wg]
    psf_g2 = psfs.hsm_g2[wg]

    dg1 = star_g1 - psf_g1
    dg2 = star_g2 - psf_g2

    T  = 2.0*(psfs.hsm_sig[wg]**2)
    Tpsf = 2.0*(stars.hsm_sig[wg]**2)
    dT = T-Tpsf; dTT = dT/T

    # Stars & size-residual-scaled stars
    starcat = treecorr.Catalog(x=x[wg], y=y[wg], g1=star_g1, g2=star_g2)
    rs_starcat = treecorr.Catalog(x=x[wg], y=y[wg], g1=star_g1*dTT, g2=star_g2*dTT)

    # PSFs & size-residual-scaled PSFs
    psfcat = treecorr.Catalog(x=x[wg], y=y[wg], g1=psf_g1, g2=psf_g2)
    rs_psfcat = treecorr.Catalog(x=x[wg], y=y[wg], g1=psf_g1*dTT, g2=psf_g2*dTT)

    # PSF Resids
    psf_resid_cat = treecorr.Catalog(x=x[wg], y=y[wg], g1=dg1, g2=dg2)

    # rho-1: psf_ellip residual autocorrelation
    rho1 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    rho1.process(psf_resid_cat)
    rho1.write(os.path.join(outdir,'_'.join([str(psf_type),'rho_1.txt'])))
    if vb:
        print('bin_size = %.6f' % rho1.bin_size)
        print('mean rho1 = %.4e median = %.4e std = %.4e' %
            (np.mean(rho1.xip),np.median(rho1.xip),
            np.std(rho1.xip)))

    # rho-2: psf_ellip x psf_ellip residual correlation
    rho2 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    rho2.process(starcat,psf_resid_cat)
    rho2.write(os.path.join(outdir,'_'.join([str(psf_type),'rho_2.txt'])))
    if vb:
        print('mean rho2 = %.4e median = %.4e std = %.4e' %
            (np.mean(rho2.xip),np.median(rho2.xip),
            np.std(rho2.xip)))

    # My *guess* at rho-3: psf_ellip x vignet_size residual correlation
    rho3 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    rho3.process(rs_starcat)
    rho3.write(os.path.join(outdir,'_'.join([str(psf_type),'rho_3.txt'])))
    if vb:
        print('mean rho3 = %.4e median = %.4e std = %.4e' %
            (np.mean(rho3.xip),np.median(rho3.xip),
            np.std(rho3.xip)))

    # My *guess* at rho-4: psf ellip resid x (psf ellip *size resid)
    rho4 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    rho4.process(psf_resid_cat,rs_starcat)
    rho4.write(os.path.join(outdir,'_'.join([str(psf_type),'rho_4.txt'])))
    if vb:
        print('mean rho4 = %.4e median = %.4e std = %.4e' %
        (np.mean(rho4.xip),np.median(rho4.xip), np.std(rho4.xip)))

    # My *guess* at rho-4: psf ellip  x (psf ellip *size resid)
    rho5 = treecorr.GGCorrelation(min_sep=min_sep, max_sep=max_sep, nbins=nbins)
    rho5.process(starcat,rs_starcat)
    rho5.write(os.path.join(outdir,'_'.join([str(psf_type),'rho_5.txt'])))
    if vb:
        print('mean rho5 = %.4e median = %.4e std = %.4e' %
        (np.mean(rho5.xip),np.median(rho5.xip), np.std(rho5.xip)))

    return rho1, rho2, rho3, rho4, rho5
