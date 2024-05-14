###
### This is a quick and dirty script to get PSF models
### for JWST data in "i2d" format and output an augmented catalog with PSF
### renderings and ERR cutouts for each galaxy
###

import numpy as np
import os, re
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
from astropy.wcs import WCS

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('images', nargs='+',
                        help = 'Images to process (wildcards OK)')
    parser.add_argument('-config', default=None,
                        help = 'Configuration file for this script')
    parser.add_argument('-outdir', default=None,
                        help = 'Where to save output')
    parser.add_argument('--vb', action='store_true', default=True,
                        help = 'Print detailed messages [does nothing for now]')

    return parser.parse_args()

def _load_fits(cat_file, run_config):
    """ Utility function to load in FITS catalog """
    imcat_fits = fitsio.FITS(cat_file, 'rw')
    cat_hdu = run_config['input_catalog']['hdu']
    imcat = imcat_fits[cat_hdu].read()
    return imcat_fits, imcat

def run_sextractor(image_file, star_config, run_config):
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

    # Try accessing the filter name from the header...
    try:
        filter_name = fits.getval(image_file, 'FILTER', ext=0)
    # Otherwise, guess filter from the file name
    except:
        filter_name = re.search(r"f(\d){3}w", image_file).group().upper()

    # Establish filter params
    filter_params = star_config['filter_params']
    star_fwhm = filter_params[filter_name]['fwhm']

    # Slightly hacky way of ascertaining whether this is an "i2d" or not
    if 'sci' in image_file:
        weight_file = image_file.replace('sci', 'wht')
    else:
        weight_file = image_file

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

        cat_matcher = htm.Matcher(17,
                       ra=imcat['ra_corr'],
                       dec=imcat['dec_corr']
                       )

        starind, catind, dist = cat_matcher.match(
                                ra=ref_stars[truth_config['ra_key']],
                                dec=ref_stars[truth_config['dec_key']],
                                radius=1./3600., maxmatch=1
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
            (imcat[star_config['mag_key']] > filt_params['min_mag']) & \
            (imcat[star_config['mag_key']] < filt_params['max_mag'])

    # Return good stars
    return imcat[wg_stars]

def _exclude_satpix(starcat, ext='DQ_VIGNET', sentinel=[1, 2]):
    """ Little utility script to exclude stars with saturated pixels.
    There are a several extensions, e.g., DQ and VAR_READNOISE that can be used.
    These values might be good for the star config? Also, subselect so that
    only the center of the star gets used.
    TO DO: this should probably be integrated into BoxCutter. Also, remove loop

    Inputs
        starcat: star catalog
        ext: vignet in starcat that acts as bad pixel mask
        sentinel: value(s) of bad pixel mask that mark pixels as bad
    """

    if len(starcat) == 0:
        print("exclude_satpix: no matched stars, returning empty")
        return starcat
    else:
        stars = starcat[ext]

    # We only care about center of star, say 25 x 25 pixels
    bs = size = stars[0].shape[0]//2
    n = bs - 11
    substars = stars[:, n:-n, n:-n]

    all_good = np.full(len(stars), True)
    for i,substar in enumerate(substars):
        # This line picks out "good" star vignets whose (unraveled) intersection
        # with the sentinel values is empty, so the length of the list is 0
        is_sat = np.size(np.intersect1d(substar, sentinel)) == 0
        all_good[i] *= is_sat

    print(f'Removed {len(all_good)-np.count_nonzero(all_good)} entries',
          f'for having value in {sentinel}')

    return starcat[all_good]


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
        print(f'\nMaking star catalog from {cat_file}\n')

    imcat_fits = fitsio.FITS(cat_file, 'rw')
    cat_hdu = run_config['input_catalog']['hdu']
    imcat = imcat_fits[cat_hdu].read()

    # Get a star catalog!
    selected_stars = _select_stars_for_psf(
                     imcat=imcat,
                     star_config=star_config,
                     run_config=run_config,
                     filter_name = filter_name
                     )

    # Filter out saturated stars
    badstar_ext = 'ERR_VIGNET'; sentinel = 0
    #badstar_ext = 'DQ_VIGNET'; sentinel = [1, 2]

    selected_stars = _exclude_satpix(selected_stars,
                                     ext=badstar_ext,
                                     sentinel=sentinel
                                     )
    if len(selected_stars) == 0:
        raise ValueError("make_starcat: No good stars found!")

    # Make size-mag plot
    size_mag_plots(imcat, selected_stars, plot_name, filter_name)

    if run_config['split_stars_validation_training'] == True:

        valid_name = starcat_name.replace('starcat', 'valid_starcat')
        full_ind = np.arange(len(selected_stars))

        rng = np.random.default_rng()
        train_ind = rng.choice(full_ind,
                               int(0.9*len(full_ind)),
                               replace=False
                               )
        valid_ind = np.setdiff1d(full_ind, train_ind)
        print(f'train ind is {train_ind}')
        print(f'valid ind is {valid_ind}')

        # Save training
        with fitsio.FITS(starcat_name, 'rw', clobber=True) as fc:
            #fc.write_table(data=imcat_fits[0].read())
            fc.write_table(data=imcat_fits[1].read(),
                             extname='LDAC_IMHEAD')
            fc.write_table(data=selected_stars[train_ind],
                             extname='LDAC_OBJECTS')
        fc.close()

        with fitsio.FITS(valid_name, 'rw', clobber=True) as fc:
            #fc.write_table(data=imcat_fits[0].read())
            fc.write_table(data=imcat_fits[1].read(),
                             extname='LDAC_IMHEAD')
            fc.write_table(data=selected_stars[valid_ind],
                             extname='LDAC_OBJECTS')
        fc.close()

        # Save a votable for the sake of it
        Table(selected_stars[train_ind]).write(
            starcat_name.replace('.fits', '.vot'),
            format='votable',
            overwrite=True
            )

        # Return the name of training catalog
        print(f'Returning training catalog')
        return starcat_name

    else:
        # Save to file
        imcat_fits['LDAC_OBJECTS'].data = selected_stars
        imcat_fits.close()

        # Return default star catalaog name
        print(f'Returning full star catalog')
        return starcat_name

def convert_pix2wcs(imfile, imcat_fits, run_config):
    """
    SExtractor's astrometry doesn't seem to handle JWST cal file
    astrometry, so use astropy.wcs to convert X, Y to RA, Dec.
    Inputs
        imfile: image file to use for WCS
        imcat_fits: fitsio.FITS instance of catalog
        run_config: main run configuration dict (as opposed to star config)
    """

    # Get header & build a WCS from it
    hdr = fits.getheader(imfile, ext=run_config['sci_image']['hdu'])
    w = WCS(hdr)

    # Load in catalog
    cat_hdu = run_config['input_catalog']['hdu']
    imcat = imcat_fits[cat_hdu].read()

    # identify X & Y keys and convert to astrometric coordinates
    x_col = run_config['input_catalog']['psf_x_key']
    y_col = run_config['input_catalog']['psf_y_key']

    coords = w.all_pix2world(
             np.array([imcat[x_col], imcat[y_col]]).T, 1
             )

    # Add corrected RA & Dec to catalog
    print(f'Adding corrected RA/Dec columns to catalog')
    imcat_fits[cat_hdu].insert_column('ra_corr', coords[:, 0])
    imcat_fits[cat_hdu].insert_column('dec_corr', coords[:, 1])

    # Return imcat_fits, now with corrected RA, Dec column
    return imcat_fits

def add_cutouts(image_file, cat_file, boxcut, run_config,
                ext='ERR', add_pix2wcs=True):
    '''
    Wrapper to call BoxCutter.grab_boxes and add an extra ERR stamp (or other!)
    for chi2 calculations. Adding the extra column to the FITS HDU List is
    roundabout, but I couldn't figure out a better way to do it.

    Inputs
        image_file: the error file to add to catalog
        cat_file: catalog that's going to get an extra vignet/stamp
        boxcut: should be a BoxCutter instance
        ext: what extension are we reading in?
        add_pix2wcs: also add corrected RA/Dec? Default=True
    '''

    # Read in the fits file so that we can add the column ourselves
    sc_fits = fitsio.FITS(cat_file, 'rw')
    cat_hdu = run_config['input_catalog']['hdu']
    imcat = sc_fits[cat_hdu].read()

    # add a corrected RA/Dec column; doing it here to minimize disk I/O
    if add_pix2wcs == True:
        print("\nAdding corrected RA, Dec to image catalog...\n")
        sc_fits = convert_pix2wcs(image_file, sc_fits, run_config)

    # This is not my favorite but run it in a loop
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
                    outdir_arg, autoselect_arg])
    print(f'psfex cmd is {cmd}')
    os.system(cmd)

    # Again, if there's a better way to do this, I don't know it.
    cleanup_cmd = ' '.join([f'mv *_{os.path.basename(starcat_file)}',
                    f'*_{os.path.basename(starcat_file).replace(".fits", ".pdf")}',
                    '*.xml', psfex_outdir])
    os.system(cleanup_cmd)

    # And for convenience... save a PSFEx stars-only file
    pexcat = Table.read(outcat_name, hdu=2)
    starcat = Table.read(starcat_file, hdu=2)

    pexstar_name = os.path.join(run_config['outdir'],
                    base_name.replace('.fits', '_pex_stars.fits'))

    pex_stars = pexcat['FLAGS_PSF']==0
    starcat[pex_stars].write(pexstar_name, format='fits', overwrite=True)


    return

def run_piffy(image_file, starcat_file, run_config, echo=True):
    '''
    Run PIFF using supplied im_file and star catalog!

    Inputs:
        im_file : the exposure to characterize
        star_cat_file : catalog of stars for PSF fitting
        run_config : the run config
    '''
    sci = run_config['sci_image']['hdu']

    astro_config_dir = run_config['astro_config_dir']

    # Get the image filename root to name outputs
    base_name   = os.path.basename(image_file)

    output_name = base_name.replace('.fits','.piff')

    piff_outdir = os.path.join(run_config['outdir'], \
                    'piff-output', base_name.split('.')[0])

    # PIFF wants central RA in hours, not degrees
    ra = fits.getval(image_file, 'CRVAL1', ext=sci)/15.0
    dec = fits.getval(image_file, 'CRVAL2', ext=sci)

    # Load PIFF configuration file
    piff_config_arg = os.path.join(astro_config_dir, 'piff.config')

    # Now run PIFF on that image and accompanying catalog
    image_arg   = f'input.image_file_name={image_file}'
    weight_arg  = f'input.weight_hdu={run_config["weight_image"]["hdu"]}'
    coord_arg   = f'input.ra={ra} input.dec={dec}'
    psfcat_arg  = f'input.cat_file_name={starcat_file}'
    output_arg  = f'output.file_name={output_name} output.dir={piff_outdir}'

    cmd = ' '.join(['piffify', piff_config_arg,
                    image_arg, weight_arg,
                    coord_arg, psfcat_arg, output_arg
                    ])
    print('piff cmd is ' + cmd)

    if echo == False:
        os.system(cmd)

    return

def main(args):

    images = args.images

    run_config = read_yaml(args.config)

    # Adds an outdir parameter to config if it was missing
    run_config = make_outdir(run_config, arg=args.outdir)

    # Get astromatic configs
    astro_config_dir = run_config['astro_config_dir']

    # Stellar locus parameters
    star_config = read_yaml(run_config['star_param_file'])

    # Create a BoxCutter instance
    boxcut = BoxCutter(config_file=run_config)

    # Process exposures
    for j, image_file in enumerate(images):

        print(f'Working on file {image_file}...\n\n')

        try:

            # Make SExtractor catalog
            cat_file = run_sextractor(image_file=image_file,
                                      star_config=star_config,
                                      run_config=run_config)

            # Add cutouts to catalogs
            add_cutouts(image_file=image_file, cat_file=cat_file,
                        boxcut=boxcut, run_config=run_config, add_pix2wcs=True)

            # Make star catalog
            starcat_file = make_starcat(image_file=image_file,
                           cat_file=cat_file, star_config=star_config,
                           run_config=run_config)

            # Run PSFEx to get PSF model
            run_psfex(image_file=image_file, starcat_file=starcat_file,
                      run_config=run_config, star_config=star_config)

            # Add MIRAGE, PIFF, WebbPSF, PSFEx, ... models to star catalog
            renderer = PSFRenderer(image_file=image_file,
                                   cat_file=starcat_file,
                                   run_config=run_config)
            renderer.render()

            # And also to galaxy catalog, if that's what's up
            if run_config['augment_galcat'] == True:
                renderer = PSFRenderer(image_file=image_file,
                                       cat_file=cat_file,
                                       run_config=run_config)
                renderer.render()

            # Add MIRAGE, PIFF, WebbPSF, PSFEx, ... models to validation cat
            if run_config['split_stars_validation_training'] == True:
                validcat_file = starcat_file.replace('starcat', 'valid_starcat')
                val_renderer = PSFRenderer(image_file=image_file,
                                           cat_file=validcat_file,
                                           run_config=run_config)
                val_renderer.render()


        except:
            print("\nWarning:")
            print(f"PSFEx and/or rendering failed for {image_file},",
                  "skipping analysis of this star...\n")

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
