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

# Local imports
from src.utils import set_rc_params

def plot_rho_stats(rho1, rho2, rho3, rho4, rho5, pixel_scale, outname=None):

    fontsize = 16
    set_rc_params(fontsize)

    fig, axes = plt.subplots(
    nrows=2, ncols=1, figsize=[9,7], sharex=True, tight_layout=True
    )

    ##
    ## rho1 correlation: dg x dg
    ##
    r = np.exp(rho1.meanlogr) * pixel_scale / 60
    xip = rho1.xip
    sig = np.sqrt(rho1.varxip)

    lab1 = r'$\rho_1(\theta)$'
    lp1 = axes[0].plot(r, xip, color='tab:blue',marker='o', ls='-', lw=1.5, label=lab1)
    axes[0].plot(r, -xip, color='tab:blue', marker='o',ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:blue', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:blue', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:blue', lw=1.5, capsize=5)
    axes[0].set_ylabel(r'$\rho(\theta)$', fontsize=fontsize)
    axes[0].set_yscale('log', nonpositive='clip')

    ##
    ## rho3 correlation: dg x dg
    ##
    r = np.exp(rho3.meanlogr) * pixel_scale / 60
    xip = rho3.xip
    sig = np.sqrt(rho3.varxip)

    lab3 = r'$\rho_3(\theta)$'
    lp3 = axes[0].plot(r, xip, color='tab:orange',marker='o', ls='-', lw=1.5, label=lab3)
    axes[0].plot(r, -xip, color='tab:orange', marker='o', ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:orange', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:orange', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:orange', lw=1.5, capsize=5)

    ##
    ## rho4 correlation: dg x dg
    ##
    r = np.exp(rho4.meanlogr) * pixel_scale / 60
    xip = rho4.xip
    sig = np.sqrt(rho4.varxip)

    lab4 = r'$\rho_4(\theta)$'
    lp4 = axes[0].plot(r, xip, color='tab:green', marker='o', ls='-', lw=1.5, label=lab4)
    #axes[0].plot(r, -xip, color='tab:green', marker='o', markerfacecolor='white', ls=':', markersize=8)
    axes[0].plot(r, -xip, color='tab:green', marker='o', ls=':')
    axes[0].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:green', ls='', capsize=5)
    #axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], marker='o', \markerfacecolor='white', markersize=8, color='tab:green', ls='', capsize=5)
    axes[0].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], marker='o', color='tab:green', ls='', capsize=5)
    axes[0].errorbar(-r, xip, yerr=sig, color='tab:green', lw=1.5, capsize=5)

    axes[0].legend([lp1, lp3, lp4], fontsize=fontsize, loc='upper right')
    #axes[0].legend(fontsize=fontsize, loc='upper right')
    axes[0].minorticks_on()

    ##
    ## rho 2 correlation: g x dg
    ##
    r = np.exp(rho2.meanlogr) * pixel_scale / 60
    xip = rho2.xip
    sig = np.sqrt(rho2.varxip)

    lab2=r'$\rho_2(\theta)$'
    #lp2 = axes[1].plot(r, abs(xip), color='tab:cyan',marker='o', ls=':')
    lp2 = axes[1].plot(r, xip, color='tab:cyan',marker='o', ls='-', lw=1.5, label=lab2)
    axes[1].plot(r, -xip, color='tab:cyan', marker='o', ls=':')
    axes[1].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:cyan', ls='', capsize=5)
    axes[1].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:cyan', ls='', capsize=5)
    axes[1].errorbar(-r, xip, yerr=sig, color='tab:cyan', lw=1.5, capsize=5)

    axes[1].set_xlabel(r'$\theta$ (arcmin)', fontsize=fontsize)
    axes[1].set_ylabel(r'$\rho(\theta)$', fontsize=fontsize)
    axes[1].set_xscale('log')
    axes[1].set_yscale('log', nonpositive='clip')

    ##
    ## rho5 correlation
    ##
    r = np.exp(rho5.meanlogr) * pixel_scale / 60.
    xip = rho5.xip
    sig = np.sqrt(rho5.varxip)

    lp5 = axes[1].plot(r, xip, color='tab:purple',marker='o', ls='-', lw=1.5, label=r'$\rho_5(\theta)$')
    axes[1].plot(r, -xip, color='tab:purple', marker='o', ls=':',)
    axes[1].errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='tab:purple', ls='', capsize=5)
    axes[1].errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='tab:purple', ls='', capsize=5)
    axes[1].errorbar(-r, xip, yerr=sig, color='tab:purple', lw=1.5, capsize=5)

    axes[1].legend([lp2, lp5], fontsize=fontsize)
    plt.legend(loc='upper right')
    axes[1].minorticks_on()

    fig.savefig(outname)
    fig.savefig(outname.replace('png', 'pdf'))

    return

## This function probably doesn't work anymore...
def compare_rho_stats(prefix, pixel_scale, file_path='./', rho_files=None):
    '''
    Make rho ratio plots for different PSF types.
    Note that the master_psf_diagnostics.py file nomenclature is assumed:
    [psf_type]_rho_[1-5].txt with psf_type={'epsfex', 'gpsfex', 'piff'}
    '''

    set_rc_params(fontsize=15)
    plt.style.use('dark_background')
    plt.rcParams.update({'figure.facecolor':'w'})

    print(f'Looking for rho-stat files in {file_path}')

    for i in range(1,6):

        try:
            pexn=os.path.join(file_path,''.join(['epsfex_rho_',str(i),'.txt']))
            pex=Table.read(pexn, format='ascii', header_start=1)
        except:
            print('no pex found')
            pex=None
        try:
            gpsfn = os.path.join(file_path,''.join(['gpsfex_rho_',str(i),'.txt']))
            gpsf=Table.read(gpsfn, format='ascii', header_start=1)
        except:
            print('no gpsf found')
            gpsf=None
        try:
            piffn = os.path.join(file_path,''.join(['piff_rho_',str(i),'.txt']))
            piff=Table.read(piffn, format='ascii', header_start=1)
        except:
            print('no piff found')
            piff=None

        try:
            singlen = os.path.join(file_path,''.join(['single_rho_',str(i),'.txt']))
            single=Table.read(singlen, format='ascii', header_start=1)
        except:
            print('no single_model found')
            single=None

        try:
            webbpsfn = os.path.join(file_path,''.join(['webbpsf_rho_',str(i),'.txt']))
            webbpsf=Table.read(webbpsfn, format='ascii', header_start=1)
        except:
            print('no single_model found')
            webbpsf=None

        savename = os.path.join(file_path,'rho_{}_comparisons.png'.format(i))

        fig,ax=plt.subplots(nrows=1,ncols=1,figsize=[10,6])
        ax.set_xscale('log')
        ax.set_yscale('log', nonpositive='clip')
        ax.set_ylabel(r'$\rho_{}(\theta)$'.format(i))
        ax.set_xlabel(r'$\theta$ (arcsec)')
        legends = []

        if pex is not None:

            r = pex['meanr'] * pixel_scale / 60 # from pixels --> arcminutes
            pex_xip = pex['xip']
            pex_sig = pex['sigma_xip']

            lp = ax.plot(r, pex_xip, color='C3', marker='o', ls='-', lw=2,
                            label='pex')
            ax.plot(r, -pex_xip, color='C3',  marker='o',ls=':')
            ax.errorbar(r[pex_xip>0], pex_xip[pex_xip>0], color='C3',
                            yerr=pex_sig[pex_xip>0], capsize=5, ls='')
            ax.errorbar(r[pex_xip<0], -pex_xip[pex_xip<0], color='C3',
                            yerr=pex_sig[pex_xip<0], capsize=5, ls='')
            ax.errorbar(-r, pex_xip, yerr=pex_sig, capsize=5, color='C3')

            legends.append(lp)

        if gpsf is not None:

            r = gpsf['meanr'] * pixel_scale / 60 # from pixels --> arcminutes
            gpsf_xip = gpsf['xip']
            gpsf_sig = gpsf['sigma_xip']

            lp2 = ax.plot(r, gpsf_xip, color='C4', marker='o', lw=2,
                            ls='-', label='gpsf')
            ax.plot(r, -gpsf_xip, color='C4', lw=2, marker='o',ls=':')
            ax.errorbar(r[gpsf_xip>0], gpsf_xip[gpsf_xip>0], yerr=gpsf_sig[gpsf_xip>0],
                            capsize=5, color='C4', ls='')
            ax.errorbar(r[gpsf_xip<0], -gpsf_xip[gpsf_xip<0], yerr=gpsf_sig[gpsf_xip<0],
                            capsize=5, color='C4', ls='')
            ax.errorbar(-r, gpsf_xip, yerr=gpsf_sig, capsize=5, color='C4')

            legends.append(lp2)

        if piff is not None:

            r = piff['meanr'] * pixel_scale / 60 # from pixels --> arcminutes
            piff_xip = piff['xip']
            piff_sig = piff['sigma_xip']

            lp3 = ax.plot(r, piff_xip, color='C5', marker='o',lw=2,
                            ls='-', label='piff')
            ax.plot(r, -piff_xip, color='C5',  marker='o', lw=2, ls=':')
            ax.errorbar(r[piff_xip>0], piff_xip[piff_xip>0], color='C5',
                            yerr=piff_sig[piff_xip>0], capsize=5, ls='')
            ax.errorbar(r[piff_xip<0], -piff_xip[piff_xip<0], color='C5',
                            yerr=piff_sig[piff_xip<0], capsize=5,  ls='')
            ax.errorbar(-r, piff_xip, yerr=piff_sig, capsize=5, color='C5')

            legends.append(lp3)

        if single is not None:

            r = single['meanr'] * pixel_scale / 60 # from pixels --> arcminutes
            single_xip = single['xip']
            single_sig = single['sigma_xip']

            lp4 = ax.plot(r, single_xip, color='C6', marker='o',lw=2,
                            ls='-', label='single')
            ax.plot(r, -single_xip, color='C6',  marker='o', lw=2, ls=':')
            ax.errorbar(r[single_xip>0], single_xip[single_xip>0], color='C6',
                            yerr=single_sig[single_xip>0], capsize=5, ls='')
            ax.errorbar(r[single_xip<0], -single_xip[single_xip<0], color='C6',
                            yerr=single_sig[single_xip<0], capsize=5,  ls='')
            ax.errorbar(-r, single_xip, yerr=single_sig, capsize=5, color='C6')

            legends.append(lp4)

        if webbpsf is not None:

            r = webbpsf['meanr'] * pixel_scale / 60 # from pixels --> arcminutes
            webbpsf_xip = np.abs(webbpsf['xip'])
            webbpsf_sig = webbpsf['sigma_xip']

            lp5 = ax.plot(r, webbpsf_xip, color='C7', marker='o',lw=2,
                            ls='-', label='webbpsf')
            ax.plot(r, -webbpsf_xip, color='C7',  marker='o', lw=2, ls=':')
            ax.errorbar(r[webbpsf_xip>0], webbpsf_xip[webbpsf_xip>0], color='C7',
                            yerr=webbpsf_sig[webbpsf_xip>0], capsize=5, ls='')
            ax.errorbar(r[webbpsf_xip<0], -webbpsf_xip[webbpsf_xip<0], color='C7',
                            yerr=webbpsf_sig[webbpsf_xip<0], capsize=5,  ls='')
            ax.errorbar(-r, webbpsf_xip, yerr=webbpsf_sig, capsize=5, color='C7')

            legends.append(lp5)


        plt.legend()
        fig.savefig(savename)

    print("rho comparison plots done")
