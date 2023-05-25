import psfex
import galsim,galsim.des
import treecorr
import piff
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.table import Table
import pdb,ipdb
from argparse import ArgumentParser
import yaml

from src.starmaker import StarMaker, StampBackground
from src.psfmaker import PSFMaker
from src.plotter import compare_rho_stats

def parse_args():

    parser = ArgumentParser()

    # I/O file names
    #parser.add_argument('config',type=str, default='psf_diagnostic_config.yaml'
    #                    help='Configuration file')
    parser.add_argument('basedir',type=str,
                        help='Directory containing star catalogs & images')
    parser.add_argument('star_cat',type=str,
                        help='Star catalog to use for PSF diagnostic')
    parser.add_argument('-min_snr',type=float, default=None,
                        help='Optional S/N cut for star catalog [default=None]')
    parser.add_argument('-pix_scale',type=float, default=None,
                        help='Image/PSF pixel scale [default=0.03]')
    parser.add_argument('-vignet_size',type=float, default=None,
                        help='Image/PSF pixel scale [default=26]')
    parser.add_argument('-outdir',type=str, default=None,
                        help='Output directory for diagnostics [default=./psf_diagnostics]')
    parser.add_argument('-psfex_name',type=str, default=None,
                        help='PSFEx model filename')
    parser.add_argument('-im_name',type=str, default=None,
                        help='FITS image filename for GalSim PSFEx diagnostic')
    parser.add_argument('-piff_name',type=str, default=None,
                        help='PIFF psf model filename')
    # Select which diagnostics to run
    parser.add_argument('--epsfex',action='store_true', default=False,
                        help='Run esheldon psfex diagnostic')
    parser.add_argument('--gpsfex',action='store_true', default=False,
                        help='Run galsim.des_psfex diagnostic')
    parser.add_argument('--piff',action='store_true', default=False,
                        help='Run PIFF diagnostic')
    parser.add_argument('--noisefree',action='store_true',default=False,
                        help='Disable adding noise to PSF stamps')
    parser.add_argument('--verbose','-v',action='store_true', default=False,
                        help='Verbosity')

    return parser.parse_args()

def get_config_file(config_file):
    "Read config yaml file"
    with open(config_file, 'r') as stream:
        try:
            config_data = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config_data


def make_output_table(makers, prefix,
                        outfile='hsm_fit_result.fits'):
    '''
    Concatenate arbitrary number of Maker() objects with HSM fits
    into an output FITS table & save to file

    : data :   list of Maker() instances
    : prefix : list of prefixes for column names
    '''

    # Bit of sanity checking
    assert type(makers) == list
    assert type(prefix) == list
    assert type(makers[0]) in [PSFMaker,StarMaker]

    mtab = {}
    mtab['x'] = makers[0].x
    mtab['y'] = makers[1].y

    # First, go through and make sub_tables:
    for i,maker in enumerate(makers):
        mtab['_'.join([prefix[i],'hsm_sig'])] = maker.hsm_sig
        mtab['_'.join([prefix[i],'hsm_g1'])] = maker.hsm_g1
        mtab['_'.join([prefix[i],'hsm_g2'])] = maker.hsm_g2
        mtab['_'.join([prefix[i],'fwhm'])] = maker.fwhm

    t = Table(mtab)
    t.write(outfile,format='fits',overwrite=True)

    return t


def main(args):
    '''
    config_file = args.config
    config = get_config_file("config.yaml")
    '''

    basedir = args.basedir
    outdir = args.outdir
    star_cat = args.star_cat
    min_snr = args.min_snr
    pix_scale = args.pix_scale
    vignet_size = args.vignet_size
    im_name = args.im_name
    psf_name = args.psfex_name
    piff_name = args.piff_name
    run_piff = args.piff
    run_gpsf = args.gpsfex
    run_pex  = args.epsfex
    noisefree = args.noisefree
    vb = args.verbose

    rho_params={'min_sep':600,'max_sep':20000,'nbins':12}

    if outdir is None:
        outdir = './psf_diagnostics_plots'
    if not os.path.isdir(outdir):
        cmd = 'mkdir -p %s' % outdir
        os.system(cmd)
    if pix_scale is None:
        pix_scale = 0.03
        print(f'Using default image/PSF pixel scale {pix_scale}')

    star_cat = Table.read(os.path.join(basedir, star_cat), hdu='LDAC_OBJECTS')
    #star_cat = Table.read(os.path.join(basedir, star_cat))

    if min_snr is not None:
        print(f"selecting S/N > {min_snr:.1f} stars")
        wg = star_cat['SNR_WIN'] > min_snr
        star_cat = star_cat[wg]

    # Calculate star stamp background
    cs = StampBackground(star_cat)
    sky_bg, sky_std = cs.calc_star_bkg(vb=vb)

    # Do star HSM fits
    sm = StarMaker(star_cat=star_cat,
                    pix_scale=pix_scale,
                    vignet_size=vignet_size,
                    vb=vb
                    )
    sm.run(bg_obj=cs, vb=vb)

    # makers stores the StarMaker and PSFMaker objects
    # prefix stores the prefix to use for saving output files
    makers = []; makers.append(sm)
    prefix = []; prefix.append('star')

    # Render PSFs, do HSM fits, save diagnostics to file
    if run_pex==True:
        pex = psfex.PSFEx(
                    os.path.join(basedir, psf_name)
                    )
        psf_pex = PSFMaker(psf_file=pex,
                            psf_type='epsfex',
                            pix_scale=pix_scale,
                            vignet_size=vignet_size,
                            noisefree=noisefree,
                            rho_params=rho_params
                            )
        psf_pex.run_all(stars=sm, vb=vb, outdir=outdir)
        makers.append(psf_pex); prefix.append('pex')

    if run_gpsf==True:
        psfex_des = galsim.des.DES_PSFEx(
                    os.path.join(basedir, psf_name),
                    os.path.join(basedir, im_name)
                    )
        psf_des = PSFMaker(psf_file=psfex_des,
                            psf_type='gpsfex',
                            pix_scale=pix_scale,
                            vignet_size=vignet_size,
                            noisefree=noisefree,
                            rho_params=rho_params
                            )
        psf_des.run_all(stars=sm, vb=vb, outdir=outdir)
        makers.append(psf_des); prefix.append('gpsf')

    if run_piff==True:
        piff_psf = piff.read(os.path.join(basedir, piff_name))
        psf_piff = PSFMaker(psf_file=piff_psf,
                                psf_type='piff',
                                pix_scale=pix_scale,
                                noisefree=noisefree,
                                rho_params=rho_params
                                )
        psf_piff.run_all(stars=sm, vb=vb, outdir=outdir)
        makers.append(psf_piff); prefix.append('piff')

    # Write star & psf HSM fits to file
    outfile = os.path.join(outdir, 'star+psf_HSMfit.fits')
    make_output_table(makers, prefix, outfile=outfile)

    # Compare rho stats for each type of PSF model fit
    compare_rho_stats(prefix, pixel_scale=pix_scale, file_path=outdir)

    return 0

if __name__ == "__main__":

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('master_psf_diagnostic.py has completed succesfully')
    else:
        print(f'master_psf_diagnostic.py has failed w/ rc={rc}')
