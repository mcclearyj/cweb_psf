import psfex
import galsim,galsim.des
import treecorr
import piff
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
rc('font',**{'family':'serif'})
rc('text', usetex=True)
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import os, re


def set_rc_params():
    '''
    Set figure parameters
    This should be a config one day
    '''
    plt.rcParams.update({'figure.facecolor':'w'})
    plt.rcParams.update({'axes.linewidth': 1.3})
    plt.rcParams.update({'xtick.labelsize': 16})
    plt.rcParams.update({'ytick.labelsize': 16})
    plt.rcParams.update({'xtick.major.size': 8})
    plt.rcParams.update({'xtick.major.width': 1.3})
    plt.rcParams.update({'xtick.minor.visible': True})
    plt.rcParams.update({'xtick.minor.width': 1.})
    plt.rcParams.update({'xtick.minor.size': 6})
    plt.rcParams.update({'xtick.direction': 'out'})
    plt.rcParams.update({'ytick.major.width': 1.3})
    plt.rcParams.update({'ytick.major.size': 8})
    plt.rcParams.update({'ytick.minor.visible': True})
    plt.rcParams.update({'ytick.minor.width': 1.})
    plt.rcParams.update({'ytick.minor.size':6})
    plt.rcParams.update({'ytick.direction':'out'})

    return

def size_mag_plot(im_cat, star_cat, plot_name, filter_name):
    '''
    Save to file a size-magnitude plot with the stellar locus highlighted
    Inputs:
        im_cat: exposure catalog, either file name or Table() object
        star_cat : star catalog passed to PIFF, either file name or Table() object
        plot_name : plot file name
        filter_name : filter name (no support for 30mas/60mas rn)
    '''

    set_rc_params()

    if type(im_cat) == str:
        image_catalog = Table.read(im_cat)
    else:
        image_catalog = im_cat

    if type(im_cat) == str:
        star_catalog = Table.read(star_cat)
    else:
        star_catalog = star_cat

    fig, ax = plt.subplots(1,1, tight_layout=True, figsize=(10,8))

    ax.plot(image_catalog['MAG_AUTO'], image_catalog['FWHM_WORLD']*3600, '.', \
            label='all objects', markersize=2)
    ax.plot(star_catalog['MAG_AUTO'], star_catalog['FWHM_WORLD']*3600, '*', \
            label='selected stars', markersize=2)

    ax.set_xlabel('MAG_AUTO (zpt=30.0)', fontsize=14)
    ax.set_ylabel('FWHM_WORLD (arcsec)', fontsize=14)
    ax.set_ylim(-0.05, 0.5)
    ax.set_xlim(16, 35)

    ax.set_title(f'{str(filter_name)} SExtractor catalog', fontsize=14)
    ax.legend(markerscale=5, fontsize=14)
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
    fontsize=16

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

    fig,axs = plt.subplots(nrows=1, ncols=3, figsize=(21,7))

    norm = colors.TwoSlopeNorm(0)
    f1 = axs[0].imshow(avg_stars, norm=colors.TwoSlopeNorm(0),
                        cmap=plt.cm.seismic_r)
    axs[0].set_title('avg star HSM sigma = %.5f\ngs.calculateFWHM() = %.5f'
                        % (star_sigma,star_fwhm), fontsize=14)
    axs[0].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[0].axhline((avg_stars.shape[1]-1)*0.5,color='black')
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(f1, cax=cax)

    f2 = axs[1].imshow(avg_psfim,norm=colors.TwoSlopeNorm(0),
                        cmap=plt.cm.seismic_r)
    axs[1].set_title('avg PSF HSM sigma = %.5f\ngs.calculateFWHM() = %.5f'
                        % (psf_sigma,psf_fwhm), fontsize=14)
    axs[1].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[1].axhline((avg_stars.shape[1]-1)*0.5,color='black')

    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right",size="5%", pad=0.05)
    plt.colorbar(f2, cax=cax)

    f3 = axs[2].imshow(avg_resid,norm=colors.TwoSlopeNorm(0), cmap=plt.cm.seismic_r)
    axs[2].set_title('sum(mean resid)= %.3f\nmean=%.3e std=%.3e' %
                        (np.nansum(avg_resid), np.nanmean(avg_resid),
                            np.nanstd(avg_resid)), fontsize=14)
    axs[2].axvline((avg_stars.shape[0]-1)*0.5,color='black')
    axs[2].axhline((avg_stars.shape[1]-1)*0.5,color='black')
    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(f3, cax=cax)

    fig.suptitle(outname.split('.')[0])
    plt.tight_layout()
    plt.savefig(outname)

    return 0


def make_quiverplot(psf=None, stars=None, outname='quiverplot.png'):
    '''
    Take a table or whatever, make quiverplot with it
    Filter out failed fits first!
    '''

    set_rc_params()
    fontsize=16

    fscl = 2.355*psf.pixel_scale # convert sigma in pixels --> FWHM in arcsec

    if (psf==None) or (stars==None):
        print("Can't make quiverplot without both star and PSF HSM fits")
        sys.exit()

    wg = (psf.hsm_g1 > -9999) & (stars.hsm_g1 > -9999)

    x = stars.x[wg]; y = stars.y[wg]
    star_g1 = stars.hsm_g1[wg]
    star_g2 = stars.hsm_g2[wg]
    psf_g1  = psf.hsm_g1[wg]
    psf_g2  = psf.hsm_g2[wg]
    psf_sig = psf.hsm_sig[wg] * fscl
    star_sig = stars.hsm_sig[wg] * fscl

    mean_star_g = np.median(np.sqrt(star_g1**2+star_g2**2))
    mean_psf_g  = np.median(np.sqrt(psf_g1**2+psf_g2**2))

    mean_star_sig = np.median(star_sig)
    mean_psf_sig  = np.median(psf_sig)

    gdiff1 = star_g1 - psf_g1
    gdiff2 = star_g2 - psf_g2
    gdiff = np.sqrt(gdiff1**2+gdiff2**2)
    sig_diff = star_sig - psf_sig

    norm = colors.TwoSlopeNorm(np.median(star_sig))
    div_norm = colors.TwoSlopeNorm(np.median(sig_diff))

    fig,[ax1,ax2,ax3] = plt.subplots(
                                nrows=1, ncols=3, sharey=True,
                                figsize=[21,6], tight_layout=True)

    q_star = ax1.quiver(y, x, star_g1,
                    star_g2, star_sig, cmap='cool',
                    units='xy', angles='uv', pivot='mid',
                    headaxislength=0, headwidth=0, headlength=0,
                    norm=norm)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(q_star, cax=cax)
    ax1.set_title('avg star HSM ellip = %.5f fwhm = %.3f' %
                    (mean_star_g, mean_star_sig), fontsize=14)

    q_res = ax2.quiver(y, x, psf_g1,
                    psf_g2, psf_sig, cmap='cool',
                    units='xy', angles='uv', pivot='mid',
                    headaxislength=0, headwidth=0, headlength=0,
                    scale=q_star.scale, norm=norm)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(q_res, cax=cax)
    ax2.set_title('avg psf ellip = %.5f fwhm = %.3f' %
                    (mean_psf_g,mean_psf_sig), fontsize=14)

    q_diff = ax3.quiver(y, x, gdiff1, gdiff2,
                    sig_diff, units='xy', angles='uv',
                    pivot='mid', headaxislength=0, headwidth=0,
                    headlength=0, cmap='cool', scale=q_star.scale,
                    norm=div_norm)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(q_diff, cax=cax)
    ax3.set_title('avg psf g_hsm resid = %.5f avg psf sig_hsm diff = %.3f' %
                    (np.median(gdiff), np.median(sig_diff)), fontsize=14)

    plt.savefig(outname)

    return

def plot_rho_stats(rho1, rho2, rho3, rho4, rho5, pixel_scale, outname=None):
    ##
    ## rho1 correlation: dg x dg
    ##

    set_rc_params()
    fontsize = 16

    fig,axes=plt.subplots(nrows=2,ncols=1,figsize=[12,8], sharex=True, tight_layout=True)

    r = np.exp(rho1.meanlogr) * pixel_scale / 60
    xip = np.abs(rho1.xip)
    sig = np.sqrt(rho1.varxip)

    lab1 = r'$\rho_1(\theta)$'
    lp1 = axes[0].plot(r, xip, color='tab:blue',marker='o',ls='-',label=lab1)
    axes[0].plot(r, -xip, color='tab:blue', marker='o',ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:blue', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:blue', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:blue', capsize=5)

    #axes[0].set_xlabel(r'$\theta$ (arcmin)', fontsize=fontsize)
    axes[0].set_ylabel(r'$\xi_+(\theta)$', fontsize=fontsize)
    #axes[0].set_xscale('log')
    axes[0].set_yscale('log', nonpositive='clip')

    ##
    ## rho3 correlation: dg x dg
    ##
    r = np.exp(rho3.meanlogr) * pixel_scale / 60
    xip = np.abs(rho3.xip)
    sig = np.sqrt(rho3.varxip)

    lab3 = r'$\rho_3(\theta)$'
    lp3 = axes[0].plot(r, xip, color='tab:orange',marker='o',ls='-',label=lab3)
    axes[0].plot(r, -xip, color='tab:orange', marker='o',ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:orange', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:orange', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:orange', capsize=5)

    ##
    ## rho4 correlation: dg x dg
    ##
    r = np.exp(rho4.meanlogr) * pixel_scale / 60
    xip = np.abs(rho4.xip)
    sig = np.sqrt(rho4.varxip)

    lab4 = r'$\rho_4(\theta)$'
    lp4 = axes[0].plot(r, xip, color='tab:green',marker='o',ls='-',label=lab4)
    axes[0].plot(r, -xip, color='tab:green', marker='o',ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:green', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:green', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:green', capsize=5)

    axes[0].legend([lp1, lp3, lp4], fontsize=14)
    axes[0].legend(fontsize=14)

    ##
    ## rho 2 correlation: g x dg
    ##
    r = np.exp(rho2.meanlogr) * pixel_scale / 60
    xip = np.abs(rho2.xip)
    sig = np.sqrt(rho2.varxip)

    lp2 = axes[1].plot(r, xip, color='magenta',marker='o', ls='-', label=r'$\rho_2(\theta)$')
    axes[1].plot(r, -xip, color='magenta', marker='o', ls=':')
    axes[1].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='magenta', ls='', capsize=5)
    axes[1].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='magenta', ls='', capsize=5)
    axes[1].errorbar(-r, xip, yerr=sig, color='magenta', capsize=5)

    axes[1].set_xlabel(r'$\theta$ (arcmin)', fontsize=fontsize)
    axes[1].set_ylabel(r'$\xi_+(\theta)$', fontsize=fontsize)
    axes[1].set_xscale('log')
    axes[1].set_yscale('log', nonpositive='clip')

    ##
    ## rho5 correlation
    ##
    r = np.exp(rho5.meanlogr) * pixel_scale / 60.
    xip = rho5.xip
    sig = np.sqrt(rho5.varxip)

    lp5 = axes[1].plot(r, xip, color='darkblue',marker='o', ls='-', label=r'$\rho_5(\theta)$')
    axes[1].plot(r, -xip, color='darkblue', marker='o', ls=':',)
    axes[1].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='darkblue', ls='', capsize=5)
    axes[1].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='darkblue', ls='', capsize=5)
    axes[1].errorbar(-r, xip, yerr=sig, color='darkblue', capsize=5)

    axes[1].legend([lp2,lp5], fontsize=14)
    #axes[1].set_ylim([2e-7,5e-3])
    plt.legend(fontsize=14)

    fig.savefig(outname)


    return
