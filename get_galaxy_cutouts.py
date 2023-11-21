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
import psfex
import fitsio
from src.plotter import size_mag_plots
from src.utils import read_yaml, make_outdir
from src.box_cutter import BoxCutter
from src.psf_renderer import PSFRenderer
import copy

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('images', nargs='+',
                        help = 'Images to process (wildcards OK)')
    parser.add_argument('-config', default=None,
                        help = 'Configuration file for this script')
    parser.add_argument('--vb', action='store_true', default=True,
                        help = 'Print detailed messages [does nothing for now]')

    return parser.parse_args()


def run_sextractor(image_file, weight_file, star_config, run_config):
    """
    Run Source Extractor, including the appropriate star FWHM for
    preliminary star identification.

    Inputs
        image_file: image for SExtractor
        weight_file: weight file for SExtractor
        star_params: table with appx. size of stars in image
        run_config: configuration file with config parameters for this run
    """

    sci = run_config['sci_image']['hdu']
    wht = run_config['weight_image']['hdu']
    wht_type = run_config['weight_image']['type']
    outdir = run_config['outdir']
    configdir = run_config['astro_config_dir']

    filter_name = fits.getval(image_file, 'FILTER', ext=0)
    filter_params = star_config['filter_params']
    star_fwhm = filter_params[filter_name]['fwhm']

    img_basename = os.path.basename(image_file)

    cat_name = os.path.join(
                outdir, img_basename.replace('.fits','.cat.fits')
                )
    bkg_sub  = os.path.join(
                    outdir, img_basename.replace('.fits','.sub.fits')
                    )
    aper_name = os.path.join(
                outdir, img_basename.replace('.fits','.aper.fits')
                )

    image_arg  = f'"{image_file}[{sci}]"'
    seeing_arg = f'-SEEING_FWHM {star_fwhm}'
    weight_arg = f'-WEIGHT_IMAGE "{weight_file}[{wht}]" -WEIGHT_TYPE {wht_type}'
    name_arg   = f'-CATALOG_NAME {cat_name}'
    check_arg  = f'-CHECKIMAGE_NAME  {bkg_sub},{aper_name}'
    param_arg  = '-PARAMETERS_NAME ' + os.path.join(configdir, 'sextractor.param')
    nnw_arg    = '-STARNNW_NAME ' + os.path.join(configdir,'default.nnw')
    filter_arg = '-FILTER_NAME ' +  os.path.join(configdir,'gauss_2.5_5x5.conv')
    config_arg = '-c ' + os.path.join(configdir, 'sextractor.config')

    cmd = ' '.join([
            'sex', image_arg, weight_arg, name_arg, check_arg, param_arg,
                nnw_arg, filter_arg, seeing_arg, config_arg
                ])

    print("sex cmd is " + cmd)
    os.system(cmd)

    return cat_name

def _select_stars_for_psf(imcat, star_config, run_config, filter_name):
    """
    Method to obtain stars from SExtractor catalog using the truth catalog
    Inputs
        imcat: input catalog from which to select stars. Should be a file
        star_config: star configuration file
    """

    if filter_name not in star_config['filter_names']:
        raise Exception('No stellar locus parameters defined ',
                        f'for filter {filter_name}')
    else:
        filt_params = star_config['filter_params'][filter_name]

    # Are we using a reference star catalog?
    if star_config['truthcat']['use'] == True:
        truth_config = star_config['truthcat']
        print('using truth catalog %s' % truth_config['filename'])

        try:
            ref_stars = Table.read(truth_config['filename'],
                        format='fits', hdu=truth_config['hdu'])
        except:
            ref_stars = Table.read(truth_config['filename'],
                        format=truth_config['format'])

        # match imcat against truth star catalog
        print(f'Selecting stars from ref star cat {truth_config["filename"]} ',
              f'containing {len(ref_stars)} stars')

        star_matcher = htm.Matcher(16,
                       ra=ref_stars[truth_config['ra_key']],
                       dec=ref_stars[truth_config['dec_key']]
                       )

        catind, starind, dist = star_matcher.match(
                                ra=imcat[run_config['input_catalog']['ra_key']],
                                dec=imcat[run_config['input_catalog']['dec_key']],
                                radius=0.5/3600., maxmatch=1
                                )

        og_len = len(imcat); imcat = imcat[catind]
        wg_stars = (imcat['SNR_WIN'] > filt_params['min_snr'])

        print(f'{len(dist)}/{og_len} objects ' +
              'matched to reference (truth) star catalog \n' +
               f'{len(imcat[wg_stars])} stars passed MIN_SNR threshold')
    else:
        # Do more standard stellar locus selection
        print("Selecting stars on CLASS_STAR, SIZE and MAG...")
        wg_stars = \
            (imcat['CLASS_STAR'] > star_config['class_star_thresh']) & \
            (imcat[star_config['size_key']] > filt_params['min_size']) & \
            (imcat[star_config['size_key']] < filt_params['max_size']) & \
            (imcat[star_config['size_key']] < filt_params['min_mag']) & \
            (imcat[star_config['size_key']] > filt_params['max_mag'])

    # Return good stars
    return imcat[wg_stars]

def _exclude_satpix(starcat, ext='DQ_VIGNET', sentinel=[1, 2]):
    """ Little utility script to exclude stars with saturated pixels.
    There are a several extensions, e.g., DQ and VAR_READNOISE that can be used.
    These values might be good for the star config? Also, subselect so that
    only the center of the star gets used.
    TO DO: this should probably be integrated into BoxCutter. Also, remove loop
    """

    if len(starcat) == 0:
        print("exclude_satpix: no matched stars, returning empty")
        return starcat
    else:
        stars = starcat[ext]

    # We only care about center of star, say 31 x 31 pixels
    bs = size = stars[0].shape[0]//2
    n = bs - 15
    substars = stars[:, n:-n, n:-n]

    all_good = np.full(len(stars), True)
    for i,substar in enumerate(substars):
        is_sat = np.size(np.intersect1d(substar, sentinel)) == 0
        all_good[i] *= is_sat
    print(f'Removed {len(all_good)-np.count_nonzero(all_good)} entries',
          'for having value in {sentinel}')

    return starcat[all_good]

def split_starcat_into_training_validation(starcat, training_fraction=0.9):
    """ For PSF validation analysis, split star catalog into training and
    validation samples. For real analysis, you will want to use all stars
    when PSF fitting.
    Input
        starcat: array-like catalog of "good" stars
        training fraction: what fraction of stars to reserve for training
    Returns
        training_stars: subset of input star catalog to use for PSF fitting
        validation_stars: subset of input star catalog to use for PSF validation
    """

    # Determine the number of rows for the training and validation catalogs
    total_rows = len(starcat)
    training_rows = int(total_rows * training_fraction)
    validation_rows = total_rows - training_rows

    # Generate random indices to shuffle the extension 2 rows
    random_indices = np.random.permutation(total_rows)

    # Split the extension 2 data into training and validation data
    training_data = starcat[random_indices[:training_rows]]
    validation_data = starcat[random_indices[training_rows:]]

    return training_data, validation_data


def make_starcat(image_file, cat_file, star_config, run_config):
    """
    Create a star catalog from the SExtractor image catalog using cuts on
    the SExtractor CLASS_STAR parameter and the supplied table of star properties.
    Also plot size-magnitude diagram for image and star catalogs.

    Inputs
        image_file: image file
        star_config: table with stellar locus parameters
        run_config: the overall config file
        cat_file: input sextractor catalog

    Returns
        star_cat_file: name of the star catalog (saved to file)
    """

    img_basename = os.path.basename(image_file)
    outdir = run_config['outdir']
    filter_name = fits.getval(image_file, 'FILTER', ext=0)

    starcat_name = os.path.join(
                   outdir, img_basename.replace('.fits','_starcat.fits')
                   )
    plot_name = os.path.join(
                outdir, img_basename.replace('.fits','_sizemag.pdf')
                )

    if not os.path.exists(cat_file):
        raise Exception(f'\n\ncould not find im_cat file {cat_file}\n\n')
    else:
        imcat_fits = fits.open(cat_file)
        imcat = imcat_fits['LDAC_OBJECTS'].data

    # Get a star catalog!
    selected_stars = _select_stars_for_psf(
                     imcat=imcat,
                     star_config=star_config,
                     run_config=run_config,
                     filter_name = filter_name
                     )
    # Filter out saturated stars
    selected_stars = _exclude_satpix(selected_stars,
                                     ext='DQ_VIGNET',
                                     sentinel=[1, 2]
                                     )
    # Make size-mag plot
    size_mag_plots(imcat, selected_stars, plot_name, filter_name)

    if run_config['split_stars_validation_training'] == True:

        #train_name = starcat_name.replace('starcat', 'train_starcat')
        valid_name = starcat_name.replace('starcat', 'valid_starcat')
        full_ind = range(len(selected_stars))

        rng = np.random.default_rng()
        train_ind = rng.choice(full_ind,
                               int(0.9*len(full_ind)),
                               replace=False
                               )
        valid_ind = np.setdiff1d(full_ind, train_ind)
        print(f"train ind is {train_ind}")
        print(f"valid ind is {valid_ind}")

        # Save training
        imcat_fits['LDAC_OBJECTS'].data = selected_stars[train_ind]
        print(f"saving training imcat to {starcat_name}")
        imcat_fits.writeto(starcat_name, overwrite=True)

        # Save validation
        validation_fits = copy.deepcopy(imcat_fits)
        validation_fits['LDAC_OBJECTS'].data = selected_stars[valid_ind]
        validation_fits.writeto(valid_name, overwrite=True)

        # Return the name of training catalog
        print(f'Returning training catalog')
        return starcat_name

    else:
        # Save to file
        imcat_fits['LDAC_OBJECTS'].data = selected_stars
        imcat_fits.writeto(starcat_name, overwrite=True)

        # Return default star catalaog name
        print(f'Returning full star catalog')
        return starcat_name

def add_cutouts(image_file, cat_file, boxcut, run_config, ext='ERR'):
    '''
    Wrapper to call BoxCutter.grab_boxes and add an extra ERR stamp (or other!)
    for chi2 calculations. Adding the extra column to the FITS HDU List is
    roundabout, but I couldn't figure out a better way to do it.

    Inputs
        image_file: the error file to add to catalog
        cat_file: catalog that's going to get an extra vignet/stamp
        boxcut: should be a BoxCutter instance
        ext: what extension are we reading in?
    '''

    # Read in the fits file so that we can add the column ourselves
    sc_fits = fitsio.FITS(cat_file, 'rw')
    cat_hdu = run_config['input_catalog']['hdu']
    imcat = sc_fits[cat_hdu].read()

    # This is not my favorite but here is a loop
    cutout_list = run_config['cutout_list']

    for hdu, ext in zip(cutout_list['hdu'], cutout_list['extname']):
        # Call to grab_boxes method
        print(f"WORKING ON {hdu}, {ext}")
        boxes = boxcut.grab_boxes(image_file=image_file,
                                  cat_file=imcat,
                                  ext=ext, hdu=hdu)
        data = np.zeros(len(boxes), dtype=[(f'{ext}_VIGNET', \
                        'f4', (boxcut.box_size, boxcut.box_size))])
        # Reformat column
        for i, box in enumerate(boxes):
            data[f'{ext}_VIGNET'][i] = box
        # Add column to catalog
        sc_fits[cat_hdu].insert_column(f'{ext}_VIGNET',
                                       data[f'{ext}_VIGNET'])
    sc_fits.close()

def run_psfex(image_file, starcat_file, run_config, star_config):
    """ Run PSFEx, creating an output directory if one doesn't already exist.
    Default is to create one directory per exposure. """

    astro_config_dir = run_config['astro_config_dir']
    base_name = os.path.basename(image_file)
    outcat_name  = starcat_file.replace('cat.fits', 'psfex_cat.fits')
    psfex_outdir = os.path.join(run_config['outdir'], \
        'psfex-output', base_name.split('.')[0])
    if not os.path.isdir(psfex_outdir):
        os.system(f'mkdir -p {psfex_outdir}')
    if star_config['truthcat']['use'] == True:
        autoselect_arg = '-SAMPLE_AUTOSELECT N'
    else:
        autoselect_arg = '-SAMPLE_AUTOSELECT Y'

    psfex_config_arg = '-c ' + os.path.join(astro_config_dir,'psfex.config')
    outcat_arg = f'-OUTCAT_NAME {outcat_name}'
    outdir_arg = f'-PSF_DIR {psfex_outdir}'

    cmd = ' '.join(['psfex',
                starcat_file, psfex_config_arg, outcat_arg,
                outdir_arg, autoselect_arg]
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


def main(args):

    i2d_images = args.images

    run_config = read_yaml(args.config)

    # Adds an outdir parameter to config if it was missing
    make_outdir(run_config)

    # Get astromatic configs
    astro_config_dir = run_config['astro_config_dir']

    # Stellar locus parameters
    star_config = read_yaml(run_config['star_param_file'])

    # Create a BoxCutter instance
    boxcut = BoxCutter(config_file=run_config)

    # Process exposures
    for j, i2d in enumerate(i2d_images):

        print(f'Working on file {i2d}...\n\n')

        image_file = i2d
        weight_file = i2d

        cat_file = run_sextractor(image_file=image_file,
                   weight_file=weight_file, star_config=star_config,
                   run_config=run_config
                   )

        add_cutouts(image_file=image_file, cat_file=cat_file,
                   boxcut=boxcut, run_config=run_config)

        cat_file = os.path.join(run_config['outdir'],
                   os.path.basename(image_file).replace('.fits', '.cat.fits'))

        starcat_file = make_starcat(image_file=image_file, cat_file=cat_file,
                       star_config=star_config, run_config=run_config)

        try:
            run_psfex(image_file=image_file, starcat_file=starcat_file,
                      run_config=run_config, star_config=star_config)

            # N.B. this adds MIRAGE, PIFF, WebbPSF... too, not just PSFEx models
            renderer = PSFRenderer(image_file=image_file, cat_file=starcat_file,
                                   run_config=run_config)
            renderer.render()

        except:
            # For the sake of paper analyses, if the PSFEx fitting failed,
            # exclude it from star/galaxy catalog
            print("\nWarning:")
            print(f"PSFEx and/or rendering failed for {image_file},",
                  "skipping analysis of this star...\n")
            pass

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
