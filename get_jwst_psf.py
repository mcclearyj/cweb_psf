###
### This is a quick and dirty script to get PSF models
### for JWST data in "i2d" format.
###
import numpy as np
import os, re
from astropy.io import fits
import pdb
from astropy.table import Table, vstack, hstack
import glob
from esutil import htm
from argparse import ArgumentParser
import ipdb
from diagnostics.plotter import size_mag_plots

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('images', nargs='+',
                        help = 'Images to process (wildcards OK)')
    parser.add_argument('-outdir', default=None,
                        help = 'Where to save files')
    parser.add_argument('-configdir', default=None,
                        help = 'Location of SEx/PIFF config files')
    parser.add_argument('-truthstars', default=None,
                        help = 'Star catalog to use for PSF fits')
    parser.add_argument('--overwrite', action='store_true', default=False,
                        help = 'Overwrite sci/weight files')
    parser.add_argument('--vb', action='store_true', default=True,
                        help = 'Print detailed messages [does nothing for now]')

    return parser.parse_args()


def get_star_params(config=None):
    '''
    Return a look-up table containing approximate stellar
    FWHMs for each of the filters based on DEC2022 COSMOS-Web simulations
    This is horrendous coding and should be fixed ASAP -- try a yaml?

    Here, selecting on FLUX_RADIUS parameter in pixels. The longer wavelength
    images have a flux_rad vs. flux_auto
    '''

    filter_names = ['f115w','f150w','f277w', 'f444w']
    fwhms = [0.058, 0.0628, 0.125, 0.165]
    min_size = [0.88, 1, 1.83, 2.33]
    max_size = [1.7, 1.5, 2.75, 3.2]
    max_mag = [27.4, 28, 27, 27.1]
    min_mag = [19, 19, 19, 19]

    star_params = Table([filter_names, fwhms, min_size,\
                        max_size, max_mag, min_mag],
                        names=['filter_name', 'star_fwhm', 'min_size',\
                        'max_size', 'max_mag', 'min_mag']
                        )

    return star_params


def extract_sci_wht(i2d, outdir, overwrite=False):
    '''
    Extract and save to file the "sci" and "WHT" extensions from the input
    i2d-format image because SExtractor can't use a weight in a MEF, apparently?
    Return the science image and weight file names to be passed to SExtractor.

    Inputs:
        data_dir : directory containing images
        overwrite : if true, overwrite any sci/WHT files saved to disk
    '''

    i2d_name = os.path.basename(i2d)
    sci_name = os.path.join(outdir, i2d_name.replace('i2d', 'sci'))
    weight_name = os.path.join(outdir, i2d_name.replace('i2d', 'weight'))

    f = fits.open(i2d)
    sci_im = f['SCI']
    weight_im = f['WHT']

    try:
        sci_im.writeto(sci_name, overwrite=overwrite)
        weight_im.writeto(weight_name, overwrite=overwrite)
        print(f'saved {sci_name} and {weight_name} to file...\n\n')

    except OSError as e:
        print(e)

    return sci_name, weight_name

def crop(image_file, weight_file, wg):

    im = fits.read(image_file)
    wt = fits.read(weight_file)

    pass
    #im_trimed = im.data[]

    return image_file, weight_file


def run_sextractor(image_file, weight_file,
                    configdir, outdir, star_params):
    '''
    Run Source Extractor, including the appropriate star FWHM for
    preliminary star identification.

    Inputs
        image_file: image for SExtractor
        weight_file: weight file for SExtractor
        configdir: directory of SEx configuration files
        outdir: directory for outputs
        star_params: table with appx. size of stars in image
    '''


    img_basename = os.path.basename(image_file)

    cat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','cat.fits')
                )
    bkg_sub  = os.path.join(
                    outdir, img_basename.replace('sci.fits','sub.fits')
                    )
    sgm_name = os.path.join(
                outdir, img_basename.replace('sci.fits','sgm.fits')
                )
    filter_name = re.search(r"f(\d){3}w", img_basename).group()

    wg = np.isin(star_params['filter_name'], filter_name)
    star_fwhm = star_params[wg]['star_fwhm']
    star_fwhm = np.float64(star_fwhm)

    seeing_arg = f'-SEEING_FWHM {star_fwhm}'
    weight_arg = f'-WEIGHT_IMAGE {weight_file} -WEIGHT_TYPE MAP_WEIGHT'
    name_arg   = f'-CATALOG_NAME {cat_name}'
    check_arg  = f'-CHECKIMAGE_NAME  {bkg_sub},{sgm_name}'
    param_arg  = '-PARAMETERS_NAME ' + os.path.join(configdir, 'sextractor.param')
    nnw_arg    = '-STARNNW_NAME ' + os.path.join(configdir,'default.nnw')
    filter_arg = '-FILTER_NAME ' +  os.path.join(configdir,'default.conv')
    config_arg = '-c ' + os.path.join(configdir, 'sextractor.mock.config')

    cmd = ' '.join([
            'sex', image_file, weight_arg, name_arg, check_arg, param_arg,
                nnw_arg, filter_arg, seeing_arg, config_arg
                ])
    print("sex cmd is " + cmd)
    os.system(cmd)

    return

def make_starcat(image_file, err_file, outdir,
                truthstars=None, star_params=None, thresh=0.92):
    '''
    Create a star catalog from the SExtractor image catalog using cuts on
    the SExtractor CLASS_STAR parameter and the supplied table of star properties.
    Alternatively, can match against a reference star catalog (truthstars).
    Also plot size-magnitude diagram for image and star catalogs.

    Inputs:
        image_file : image file
        out_dir : where to save star catalog (also assumed to be location of image catalog)
        truthstars : reference star catalog
        star_params : table with stellar locus parameters
        thresh : CLASS_STAR threshold for accepting as star

    Outputs:
        star_cat_file: name of the star catalog (saved to file)

    TO DO:
        - make make_starcat robust to input catalogs that are not LDACs
        - make starparams a config file. Alternatively, improve auto-star selection
    '''

    # Make sure that one of either truthstars or star_params is supplied
    if (truthstars is None) and (star_params is None):
        print('Reference star catalog and star param table both NoneType, exiting')

    img_basename = os.path.basename(image_file)

    try:
        filter_name = re.search(r"f(\d){3}w", image_file).group() \
                        +'_'+ re.search(r"(\d){2}mas",image_file).group()
    except AttributeError:
        filter_name = re.search(r"f(\d){3}w", image_file).group()

    imcat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','cat.fits')
                )
    starcat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','starcat.fits')
                )
    plot_name = os.path.join(
                outdir, img_basename.replace('sci.fits','sizemag.pdf')
                )

    if not os.path.exists(imcat_name):
        print(f'could not find image im_cat file {imcat_name}')
    else:
        im_cat_fits = fits.open(imcat_name)
        im_cat = im_cat_fits['LDAC_OBJECTS'].data

    if truthstars is not None:

        truth_star_tab = Table.read(truthstars, format='ascii')
        truthmatch = htm.Matcher(16, ra=truth_star_tab['x_or_RA'],
                                    dec=truth_star_tab['y_or_Dec'])
        cat_ind, truth_ind, dist = truthmatch.match(
                                    ra=im_cat['ALPHAWIN_J2000'],
                                    dec=im_cat['DELTAWIN_J2000'],
                                    maxmatch=1, radius = 0.5/3600)

    else:

        filter_name = re.search(r"f(\d){3}w", image_file).group()
        wg = np.isin(star_params['filter_name'], filter_name)
        min_size = star_params[wg]['min_size']
        min_size = np.float64(min_size)
        max_size = star_params[wg]['max_size']
        max_size = np.float64(max_size)
        max_mag = star_params[wg]['max_mag']
        max_mag = np.float64(max_mag)
        min_mag = star_params[wg]['min_mag']
        min_mag = np.float64(min_mag)

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


def cookie_cutter(filename, catalog):
    '''
    Make sure that
    '''


def run_piffy(im_file, star_cat_file, configdir, outdir):
    '''
    Run PIFF using supplied im_file and star catalog!

    Inputs:
        im_file : the exposure to characterize
        star_cat_file : catalog of stars for PSF fitting
        configdir : path to PIFF config file
        outdir : where to save PIFF results
    '''

    # Get the image filename root to name outputs
    base_name   = os.path.basename(im_file)

    output_name = base_name.replace('.fits','.piff')

    piff_outdir = os.path.join(outdir, \
                    'piff-output', base_name.split('.')[0])

    bkg_sub = os.path.join(outdir,
                base_name.replace('sci.fits', 'sub.fits'))

    # PIFF wants RA in hours, not degrees
    ra = fits.getval(im_file, 'CRVAL1')/15.0
    dec = fits.getval(im_file, 'CRVAL2')

    # Load PIFF configuration file
    piff_config_arg = os.path.join(configdir, 'piff.config')

    # Now run PIFF on that image and accompanying catalog
    image_arg   = f'input.image_file_name={bkg_sub}'
    coord_arg   = f'input.ra={ra} input.dec={dec}'
    psfcat_arg  = f'input.cat_file_name={star_cat_file}'
    output_arg  = f'output.file_name={output_name} output.dir={piff_outdir}'

    cmd = ' '.join([
             'piffify', piff_config_arg, image_arg, coord_arg, \
                psfcat_arg, output_arg
                ])
    print('piff cmd is ' + cmd)
    os.system(cmd)

    return

def run_psfex(image_file, configdir, outdir):

    psfex_config_arg = '-c ' + os.path.join(configdir,'psfex.mock.config')

    img_basename = os.path.basename(image_file)
    cat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','cat.fits')
                )
    starcat_name = cat_name.replace('cat.fits', 'starcat.fits')
    outcat_name  = cat_name.replace('cat.fits', 'psfex_cat.fits')
    outcat_arg   = f'-OUTCAT_NAME {outcat_name}'

    psfex_outdir = outdir

    psfcat_name = starcat_name #cat_name

    cmd = ' '.join(
            ['psfex', psfcat_name, psfex_config_arg, outcat_arg]
            )

    print(f'psfex cmd is {cmd}')
    os.system(cmd)

    cleanup_cmd = ' '.join(
                    ['mv chi* resi* samp* snap* proto* *.xml', psfex_outdir]
                    )
    cleanup_cmd2 = ' '.join(
                    ['mv count*pdf ellipticity*pdf fwhm*pdf', psfex_outdir]
                    )
    os.system(cleanup_cmd)
    os.system(cleanup_cmd2)

    return


def main(args):

    i2d_images = args.images
    outdir = args.outdir
    configdir = args.configdir
    truthstars = args.truthstars
    overwrite = args.overwrite
    vb = args.vb
    crop = False

    # Set default output directory values if none provided
    if configdir is None:
        configdir = '/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config/'

    if outdir is None:
        basedir = os.path.commonpath(images)
        outdir = os.path.join(basedir,'working')

    # Make output directory
    if not os.path.isdir(outdir):
        cmd = 'mkdir -p {outdir}'.format(outdir=outdir)
        os.system(cmd)
        print(f'Made output directory {outdir}')

    else:
        print(f'Output directory {outdir} exists, continuing...')


    # Placeholder: table of stellar locus parameters, though probably
    # not necessary long-term -- need to develop better auto-starselect
    star_params = get_star_params()

    # Process exposures
    for i2d in i2d_images:

        print(f'Working on file {i2d}...\n\n')

        #image_file, weight_file = extract_sci_wht(i2d, outdir,
        #                                            overwrite=overwrite)

        image_file = i2d
        weight_file = i2d.replace('sci', 'wht')
        err_file = i2d.replace('sci', 'err')

        run_sextractor(image_file, weight_file,
                                configdir, outdir, star_params)

        starcat_name = make_starcat(image_file, err_file, outdir,
                                        truthstars=truthstars,
                                        star_params=star_params
                                        )

        run_piffy(image_file, starcat_name,
                    configdir=configdir, outdir=outdir)


        run_psfex(image_file, configdir, outdir)

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
