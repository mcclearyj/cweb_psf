###
### This is a quick and dirty script to get PSF models
### for JWST data in "i2d" format and output an augmented catalog with PSF
### renderings and ERR cutouts for each galaxy
###

import numpy as np
import os
from astropy.io import fits
import pdb
from astropy.table import Table, Column
import glob
from esutil import htm
from argparse import ArgumentParser
import ipdb, pdb
import psfex
from src.plotter import size_mag_plots
from src.utils import read_yaml, make_outdir
from src.box_cutter import BoxCutter

import fitsio


def parse_args():

    parser = ArgumentParser()

    parser.add_argument('images', nargs='+',
                        help = 'Images to process (wildcards OK)')
    parser.add_argument('-config', default=None,
                        help = 'Configuration file for this script')
    parser.add_argument('--vb', action='store_true', default=True,
                        help = 'Print detailed messages [does nothing for now]')

    return parser.parse_args()


def run_sextractor(image_file, weight_file, star_params, run_config):
    '''
    Run Source Extractor, including the appropriate star FWHM for
    preliminary star identification.

    Inputs
        image_file: image for SExtractor
        weight_file: weight file for SExtractor
        star_params: table with appx. size of stars in image
        run_config: configuration file with config parameters for this run
    '''

    sci = run_config['sci_image']['hdu']
    wht = run_config['weight_image']['hdu']
    outdir = run_config['outdir']
    configdir = run_config['astro_config_dir']

    filter_name = fits.getval(image_file, 'FILTER', ext=0)
    filter_params = star_params['filter_params']
    star_fwhm = filter_params[filter_name]['fwhm']

    img_basename = os.path.basename(image_file)

    cat_name = os.path.join(
                outdir, img_basename.replace('.fits','.cat.fits')
                )
    bkg_sub  = os.path.join(
                    outdir, img_basename.replace('.fits','.sub.fits')
                    )
    sgm_name = os.path.join(
                outdir, img_basename.replace('.fits','.sgm.fits')
                )

    image_arg  = f'"{image_file}[{sci}]"'
    seeing_arg = f'-SEEING_FWHM {star_fwhm}'
    weight_arg = f'-WEIGHT_IMAGE "{weight_file}[{wht}]"'
    name_arg   = f'-CATALOG_NAME {cat_name}'
    check_arg  = f'-CHECKIMAGE_NAME  {bkg_sub},{sgm_name}'
    param_arg  = '-PARAMETERS_NAME ' + os.path.join(configdir, 'sextractor.param')
    nnw_arg    = '-STARNNW_NAME ' + os.path.join(configdir,'default.nnw')
    filter_arg = '-FILTER_NAME ' +  os.path.join(configdir,'gauss_3.0_5x5.conv')
    config_arg = '-c ' + os.path.join(configdir, 'sextractor_hiSN.config')

    cmd = ' '.join([
            'sex', image_arg, weight_arg, name_arg, check_arg, param_arg,
                nnw_arg, filter_arg, seeing_arg, config_arg
                ])

    print("sex cmd is " + cmd)
    os.system(cmd)

    return cat_name


def make_starcat(image_file, cat_file, star_params, run_config):
    '''
    Create a star catalog from the SExtractor image catalog using cuts on
    the SExtractor CLASS_STAR parameter and the supplied table of star properties.
    Also plot size-magnitude diagram for image and star catalogs.

    Inputs:
        image_file : image file
        star_params : table with stellar locus parameters
        thresh : CLASS_STAR threshold for accepting a star candidate
        cat_file : input sextractor catalog

    Outputs:
        star_cat_file: name of the star catalog (saved to file)

    '''

    img_basename = os.path.basename(image_file)
    outdir = run_config['outdir']
    imcat_name = cat_file
    filter_name = fits.getval(image_file, 'FILTER', ext=0)

    starcat_name = os.path.join(
                outdir, img_basename.replace('.fits','_starcat.fits')
                )
    plot_name = os.path.join(
                outdir, img_basename.replace('.fits','_sizemag.pdf')
                )

    if not os.path.exists(imcat_name):
        raise Exception(f'\n\ncould not find im_cat file {imcat_name}\n\n')
    else:
        im_cat_fits = fits.open(imcat_name)
        im_cat = im_cat_fits['LDAC_OBJECTS'].data

    if filter_name not in star_params['filter_names']:
        raise Exception(f'No stellar locus parameters defined for filter {filter_name}')

    filter_params = star_params['filter_params']

    min_size = filter_params[filter_name]['min_size']
    max_size = filter_params[filter_name]['max_size']
    max_mag = filter_params[filter_name]['max_mag']
    min_mag = filter_params[filter_name]['min_mag']
    thresh = star_params['class_star_thresh']

    star_selec = (im_cat['CLASS_STAR'] > thresh) \
                    & (im_cat['MAG_AUTO'] < max_mag) \
                    & (im_cat['MAG_AUTO'] > min_mag) \
                    & (im_cat['FLUX_RADIUS'] > min_size) \
                    & (im_cat['FLUX_RADIUS'] < max_size)

    cat_ind = np.arange(len(star_selec))[star_selec]

    star_cat = im_cat[cat_ind]

    im_cat_fits['LDAC_OBJECTS'].data = star_cat
    im_cat_fits.writeto(starcat_name, overwrite=True)

    # Make size-mag plot
    size_mag_plots(im_cat, star_cat, plot_name, filter_name)

    return starcat_name


def add_err_cutout(image_file, cat_file, boxcut, run_config, ext='ERR'):
    '''
    Wrapper to call BoxCutter.grab_boxes and add an extra ERR stamp (or other!)
    for chi2 calculations. Adding the extra column to the FITS HDU List is
    roundabout, but I couldn't figure out a better way to do it.

    Inputs
        grabber: should be a BoxCutter instance
        image_file: the error file to add to catalog
        cat_file: catalog that's going to get an extra vignet/stamp
        ext: what extension are we reading in?
    '''

    cat_hdu = run_config['input_catalog']['hdu']
    box_size = run_config['box_size']

    # Read in the fits file so that we can add the column ourselves
    sc_fits = fitsio.FITS(cat_file, 'rw')
    imcat = sc_fits[cat_hdu].read()

    # We need box_size to be same as star size for chi2 calc
    if (box_size != imcat['VIGNET'][0].shape[0]):
        print(f'supplied box_size={box_size} and vignet size={imcat["VIGNET"][0].shape[0]} differ!!')
        print(f'overriding supplied box_size to {imcat["VIGNET"][0].shape[0]}')
        box_size = imcat['VIGNET'][0].shape[0]
        boxcut.box_size = box_size

    # Call to grab_boxes method
    boxes = boxcut.grab_boxes(image_file=image_file, cat_file=imcat)

    data = np.zeros(len(boxes), dtype=[(f'{ext}_VIGNET', \
                        'f4', (box_size, box_size))])
    for i, box in enumerate(boxes):
        data[f'{ext}_VIGNET'][i] = box
    sc_fits[cat_hdu].insert_column(f'{ext}_VIGNET', data[f'{ext}_VIGNET'])
    sc_fits.close()

    return


def run_psfex(image_file, starcat_file, run_config):
    '''
    Run PSFEx, creating an output directory if one doesn't already exist.
    Default is to create one directory per exposure.
    '''

    astro_config_dir = run_config['astro_config_dir']
    base_name = os.path.basename(image_file)
    outcat_name  = starcat_file.replace('cat.fits', 'psfex_cat.fits')
    psfex_outdir = os.path.join(run_config['outdir'], \
        'psfex-output', base_name.split('.')[0])
    if not os.path.isdir(psfex_outdir):
        os.system(f'mkdir -p {psfex_outdir}')

    psfex_config_arg = '-c ' + os.path.join(astro_config_dir,'psfex.config')
    outcat_arg = f'-OUTCAT_NAME {outcat_name}'
    outdir_arg = f'-PSF_DIR {psfex_outdir}'

    cmd = ' '.join(['psfex',
                starcat_file, psfex_config_arg, outcat_arg, outdir_arg]
            )

    print(f'psfex cmd is {cmd}')
    os.system(cmd)

    # Again, if there's a better way to do this, I don't know it.
    cleanup_cmd = ' '.join([f'mv *_{os.path.basename(starcat_file)}',
        f'*_{os.path.basename(starcat_file).replace(".fits", ".pdf")} *.xml',
        psfex_outdir])

    os.system(cleanup_cmd)

    # And for convenience... save a PSFEx stars-only file
    pexcat = Table.read(outcat_name, hdu=2)
    starcat = Table.read(starcat_file, hdu=2)

    pexstar_name = os.path.join(run_config['outdir'],
                    base_name.replace('.fits', '_pex_stars.fits'))

    pex_stars = pexcat['FLAGS_PSF']==0
    starcat[pex_stars].write(pexstar_name, format='fits', overwrite=True)

    return


def render_psf(image_file, cat_file, run_config):
    '''
    Place a PSFEx rendering into the galaxy catalog
    '''

    # Preliminaries
    basename = os.path.basename(image_file)
    psf_name = basename.replace('.fits', '_starcat.psf')
    psfex_outdir = os.path.join(run_config['outdir'], \
        'psfex-output', basename.split('.')[0])
    psf_path = os.path.join(psfex_outdir, psf_name)

    # Load in the PSF
    psf = psfex.PSFEx(psf_path)

    # Read in the fits file so that we can add the column ourselves
    gc_fits = fitsio.FITS(cat_file, 'rw')
    hdu = run_config['input_catalog']['hdu']
    gal_catalog = gc_fits[hdu].read()

    # Grab X and Y arrays
    x = gal_catalog[run_config['input_catalog']['x_tag']]
    y = gal_catalog[run_config['input_catalog']['y_tag']]

    # Create empty list that will hold the PSF renderings
    psf_images = []

    # Loop prevents weird memory overflow errors for very large catalogs
    for xi,yi in zip(x, y):
        psf_images.append(psf.get_rec(yi, xi))

    # Create PSF image array and add to catalog
    psf_im_arr = np.zeros(len(psf_images),
                    dtype=[('PSF_VIGNET', 'f4', (psf_images[0].shape))])

    for i, psf_im in enumerate(psf_images): psf_im_arr['PSF_VIGNET'][i] = psf_im

    # Add new column to galaxy catalog, save to file
    gc_fits[hdu].insert_column('PSF_VIGNET', psf_im_arr['PSF_VIGNET'])
    gc_fits.close()

    return


def main(args):

    i2d_images = args.images

    run_config = read_yaml(args.config)

    # Adds an outdir parameter to config if it was missing
    make_outdir(run_config)

    # Get astromatic configs
    astro_config_dir = run_config['astro_config_dir']

    # Stellar locus parameters
    star_params = read_yaml(run_config['star_param_file'])

    # Create a BoxCutter instance
    boxcut = BoxCutter(config_file=args.config)

    # Process exposures
    for j, i2d in enumerate(i2d_images):

        print(f'Working on file {i2d}...\n\n')

        image_file = i2d
        weight_file = i2d


        cat_file = run_sextractor(image_file=image_file,
                        weight_file=weight_file, star_params=star_params,
                        run_config=run_config)

        add_err_cutout(image_file=image_file, cat_file=cat_file,
                        boxcut=boxcut, run_config=run_config)


        cat_file = os.path.join(run_config['outdir'],
                    os.path.basename(image_file).replace('.fits', '.cat.fits'))

        starcat_file = make_starcat(image_file=image_file, cat_file=cat_file,
                        star_params=star_params, run_config=run_config)

        run_psfex(image_file=image_file, starcat_file=starcat_file,
                        run_config=run_config)

        render_psf(image_file=image_file, cat_file=cat_file,
                        run_config=run_config)

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
