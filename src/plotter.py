import psfex
import galsim,galsim.des
import treecorr
import piff
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import os, re
import pdb, ipdb
from astropy.table import Table
import fitsio

from src.utils import set_rc_params


def size_mag_plots(im_cat, star_cat, plot_name, filter_name):
    '''
    Save to file a size-magnitude plot with the stellar locus highlighted
    Inputs:
        im_cat: exposure catalog, either file name or Table() object
        star_cat : star catalog passed to PIFF, either file name or Table() object
        plot_name : plot file name
        filter_name : filter name (no support for 30mas/60mas rn)
    '''

    set_rc_params(fontsize=14)

    if type(im_cat) == str:
        image_catalog = Table.read(im_cat)
    else:
        image_catalog = im_cat

    if type(im_cat) == str:
        star_catalog = Table.read(star_cat)
    else:
        star_catalog = star_cat

    fig, axs = plt.subplots(1,2, tight_layout=True, figsize=(12, 6.5))

    # First, do FWHM
    axs[0].plot(
        image_catalog['MAG_AUTO'], image_catalog['FWHM_WORLD']*3600, '.',
        label='all objects', markersize=3
    )
    axs[0].plot(
        star_catalog['MAG_AUTO'], star_catalog['FWHM_WORLD']*3600, '.',
        label='selected stars', markersize=3
    )

    #axs[0].set_xlabel(r'\texttt{MAG_AUTO}', fontsize=16)
    axs[0].set_ylabel(r'\texttt{FWHM_WORLD} (arcsec)')
    axs[0].set_xlabel(r'\texttt{MAG_AUTO}')
    axs[0].set_ylim(-0.05, 1.1)
    axs[0].set_xlim(16, 29)
    axs[0].grid(True)
    axs[0].legend(markerscale=3, fontsize=14, loc='upper left')

    # Then, flux_radius
    axs[1].plot(
        image_catalog['MAG_AUTO'], image_catalog['FLUX_RADIUS'], '.',
        label='all objects', markersize=3
    )
    axs[1].plot(
        star_catalog['MAG_AUTO'], star_catalog['FLUX_RADIUS'], '.',
        label='selected stars', markersize=3
    )

    axs[1].set_xlabel(r'\texttt{MAG_AUTO}')
    axs[1].set_ylabel(r'\texttt{FLUX_RADIUS} (pix)')
    axs[1].set_ylim(0.8, 13)
    axs[1].set_xlim(16, 29)
    axs[1].grid(True)
    axs[1].legend(markerscale=3, fontsize=14, loc='upper left')

    fig.savefig(plot_name)

    return


def make_resid_plot(psf, stars, outname='star_psf_resid.png', vb=False):
    '''
    make figures of average stars, psf renderings,
    and residuals between the two

    :avg_psf : should be an instance of PSFMaker()
    :avg_star: should be an instance of StarMaker()
    '''

    set_rc_params()
    fontsize = 16
    vmin = 0.0001
    vmax = 75

    avg_stars = np.nanmean(stars.stamps,axis=0)
    avg_psfim = np.nanmean(psf.stamps,axis=0)
    avg_resid = np.nanmean(psf.resids,axis=0)

    if vb==True:
        print("avg_star total flux = %.3f" % np.sum(avg_stars))
        print("avg_psf total flux = %.3f" % np.sum(avg_psfim))

    # Calculate average sizes to display in image
    wg = (psf.hsm_sig > -9999.) & (stars.hsm_sig > -9999.)
    psf_fwhm  = np.nanmean(psf.fwhm[wg])
    psf_sigma = np.nanmean(psf.hsm_sig[wg])
    star_fwhm = np.nanmean(stars.fwhm[wg])
    star_sigma = np.nanmean(stars.hsm_sig[wg])

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15,7), tight_layout=True)

    f1 = axs[0].imshow(avg_stars, cmap=plt.cm.bwr_r,
                        norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    axs[0].set_title('avg star HSM sigma = %.4f\ngs.calculateFWHM() = %.4f'
                        % (star_sigma,star_fwhm), fontsize=16)
    axs[0].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[0].axhline((avg_stars.shape[1]-1)*0.5,color='black')
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(f1, cax=cax)

    f2 = axs[1].imshow(avg_psfim, cmap=plt.cm.bwr_r,
                        norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    axs[1].set_title('avg PSF HSM sigma = %.4f\ngs.calculateFWHM() = %.4f'
                        % (psf_sigma,psf_fwhm), fontsize=16)
    axs[1].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[1].axhline((avg_stars.shape[1]-1)*0.5,color='black')

    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right",size="5%", pad=0.05)
    plt.colorbar(f2, cax=cax)

    resid_norm = colors.TwoSlopeNorm(0, vmin=0.9*np.min(avg_resid),
                                        vmax=0.9*np.max(avg_resid)
                                        )
    f3 = axs[2].imshow(avg_resid,norm=resid_norm, cmap=plt.cm.seismic_r)
    axs[2].set_title('sum(mean resid)= %.3f\nmean=%.2e std=%.2e' %
                        (np.nansum(avg_resid), np.nanmean(avg_resid),
                            np.nanstd(avg_resid)), fontsize=16)
    axs[2].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[2].axhline((avg_stars.shape[1]-1)*0.5,color='black')
    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(f3, cax=cax)

    plt.savefig(outname)

    # Write fits files out to file too
    outdir = os.path.dirname(outname)
    fout = fitsio.FITS(
            os.path.join(outdir, 'mean_star_model_resid_ims.fits'), 'rw')
    fout.write([avg_stars, avg_psfim, avg_resid],
                names=['STARS', 'MODELS', 'RESIDS']
                )
