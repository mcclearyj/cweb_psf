###
### This is a quick and dirty script to get PSF models
### for JWST data in "i2d" format.
###
import numpy as np
import os, re
from astropy.io import fits
from matplotlib import rc,rcParams
rc('font',**{'family':'serif'})
rc('text', usetex=True)
import matplotlib.pyplot as plt
import pdb
from astropy.table import Table, vstack, hstack
import glob
from argparse import ArgumentParser

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('-basedir', default=None,
                        help = 'Location of simulated images')
    parser.add_argument('-outdir', default=None,
                        help = 'Where to save files')
    parser.add_argument('-configdir', default=None,
                        help = 'Location of SEx/PIFF config files')
    parser.add_argument('-star_fwhm_tab', default=None,
                        help = 'Lookup table with stellar FWHMs for each filter')
    parser.add_argument('--overwrite', action='store_true', default=False,
                        help = 'Overwrite sci/weight files')

    return parser.parse_args()


def make_fwhm_tab():
    '''
    Return a look-up table containing approximate stellar
    FWHMs for each of the filters based on DEC2022 COSMOS-Web simulations
    This is horrendous coding and should be fixed ASAP -- try a yaml?
    '''

    filter_names = ['f115w','f150w','f277w', 'f444w']
    fwhms = [0.058, 0.0628, 0.146, 0.175]
    min_fwhm_im = [1.5, 1.6, 1.6, 2.26]
    max_fwhm_im = [2.4, 2.8, 2.8, 3.5]
    star_fwhms = Table([filter_names, fwhms, min_fwhm_im, max_fwhm_im],
        names=['filter_name', 'star_fwhm', 'min_fwhm_im', 'max_fwhm_im'])

    return star_fwhms


def extract_sci_wht(i2d, data_dir=None, overwrite=False):
    '''
    Extract and save to file the "SCI" and "WHT" extensions from the input
    i2d-format imagebecause SExtractor can't SExtract a weight in a MEF, apparently.
    Return the science image and weight file names to be passed to SExtractor.

    Inputs:
        data_dir : directory containing images
        overwrite : if true, overwrite any SCI/WHT files saved to disk
    '''

    sci_name = i2d.replace('i2d', 'sci')
    weight_name = i2d.replace('i2d', 'weight')

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


def run_sextractor(image_file, weight_file,
                    configdir, outdir, star_fwhms):
    '''
    Run Source Extractor, including the appropriate star FWHM for
    preliminary star identification.

    Inputs:
        star_fwhms : table with appx. size of stars in image
        image_file : image for SExtractor
        weight_file : weight file for SExtractor
        configdir : directory of SEx configuration files
    '''

    cat_name = os.path.join(
                outdir, image_file.replace('sci.fits','cat.fits')
                )
    aperture_name = os.path.join(
                    outdir, image_file.replace('sci.fits','apertures.fits')
                    )
    sgm_name = os.path.join(
                outdir, image_file.replace('sci.fits','sgm.fits')
                )

    # Assuming that the JWST filter name appears in the image file name
    filter_name = re.search(r"f(\d){3}w", image_file).group()
    wg = np.isin(star_fwhms['filter_name'], filter_name)
    star_fwhm = star_fwhms[wg]['star_fwhm']
    star_fwhm = np.float64(star_fwhm)

    seeing_arg = f'-SEEING_FWHM {star_fwhm}'
    weight_arg = f'-WEIGHT_IMAGE {weight_file} -WEIGHT_TYPE MAP_WEIGHT'
    name_arg   = f'-CATALOG_NAME {cat_name}'
    check_arg  = f'-CHECKIMAGE_NAME  {aperture_name},{sgm_name}'
    param_arg  = '-PARAMETERS_NAME ' + os.path.join(configdir, 'sextractor.param')
    nnw_arg    = '-STARNNW_NAME ' + os.path.join(configdir,'default.nnw')
    filter_arg = '-FILTER_NAME ' +  os.path.join(configdir,'default.conv')
    config_arg = '-c ' + os.path.join(configdir, 'sextractor.mock.config')

    cmd = ' '.join([
        'sex', image_file, weight_arg, name_arg,  check_arg,\
        param_arg, nnw_arg, filter_arg, seeing_arg, config_arg\
        ])

    print("sex cmd is " + cmd)
    os.system(cmd)

    return

def get_starcat(image_file, outdir, star_fwhms=None, thresh=0.92):
    '''
    Create a star catalog from the SExtractor image catalog using cuts on
    the SExtractor CLASS_STAR parameter and the supplied table of star properties.

    Inputs:
        image_file : image file
        out_dir : where to save star catalog
        thresh : CLASS_STAR threshold for accepting as star

    Outputs:
        star_cat_file: name of the star catalog (saved to file)
    '''

    cat_name = image_file.replace('sci.fits','cat.fits')
    star_cat_name = image_file.replace('sci.fits','starcat.fits')
    cat_file = os.path.join(outdir, cat_name)
    star_cat_file = os.path.join(outdir, star_cat_name)

    if star_fwhms is not None:
        filter_name = re.search(r"f(\d){3}w", image_file).group()
        wg = np.isin(star_fwhms['filter_name'], filter_name)
        min_fwhm_im = star_fwhms[wg]['min_fwhm_im']
        min_fwhm_im = np.float64(min_fwhm_im)
        max_fwhm_im = star_fwhms[wg]['max_fwhm_im']
        max_fwhm_im = np.float64(max_fwhm_im)
    else:
        max_fwhm_im = -9999.
        max_fwhm_im = 9999.

    if os.path.exists(cat_file):
        catalog = Table.read(cat_file)
        star_selec = (catalog['CLASS_STAR'] >= thresh) \
                        & (catalog['MAG_AUTO'] < 35) \
                        & (catalog['FWHM_IMAGE'] > min_fwhm_im) \
                        & (catalog['FWHM_IMAGE'] < max_fwhm_im)
        catalog[star_selec].write(star_cat_file, format='fits', overwrite=True)

    else:
        print(f'could not find image catalog file {cat_file}')
        pdb.set_trace()

    return star_cat_file


def run_piffy(im_file, star_cat_file, configdir, outdir):
    '''
    Run PIFF using supplied im_file and star catalog!

    Inputs:
        im_file : the exposure to characterize
        star_cat_file : catalog of stars for PSF fitting
        configdir : path to PIFF config file
        outdir : where to save PIFF results
    '''

    # Load PIFF configuration file
    run_piff_config = os.path.join(configdir, 'piff.config')
    outdir = os.path.join(outdir,'piff-output')

    # Now run PIFF on that image and accompanying catalog
    image_arg   = f'input.image_file_name={im_file}'
    psfcat_arg  = f'input.cat_file_name={star_cat_file}'
    output_name = im_file.split('/')[-1].replace('.fits', '.piff')
    full_output_name = os.path.join(outdir, output_name)
    output_arg  = f'output.file_name={output_name} output.dir={outdir}'

    cmd = ' '.join([
             'piffify', run_piff_config, image_arg, psfcat_arg, output_arg
             ])

    print('piff cmd is ' + cmd)
    os.system(cmd)

    return

def main(args):

    basedir = args.basedir
    outdir = args.outdir
    configdir = args.configdir
    overwrite = args.overwrite

    # Set default parameter values if none provided
    if basedir is None:
        basedir = '/Users/j.mccleary/Research/jwst_cosmos/mock_data/DEC2022/mosaics'
    if outdir is None:
        outdir = os.path.join(basedir,'working')
    if configdir is None:
        configdir = '/Users/j.mccleary/Research/jwst_cosmos/astro_config/'

    # Make output directory
    if not os.path.isdir(outdir):
        cmd = 'mkdir -p {outdir}'.format(outdir=outdir)
        os.system(cmd)
        print('Made output directory {outdir}').format(outdir=outdir)
    else:
        print(f'Output directory {outdir} exists, continuing...')

    # Load in images
    i2d_path = os.path.join(basedir,'*_i2d.fits')
    all_i2ds = glob.glob(i2d_path)
    all_i2ds.sort()
    print(f'Found {all_i2ds}')

    # Make reference stellar FWHM table
    star_fwhms = make_fwhm_tab()

    # Process exposures
    for i2d in all_i2ds[1:2]:
        print(f'Working on file {i2d}...\n\n')
        image_file, weight_file = extract_sci_wht(i2d, overwrite=overwrite)
        run_sextractor(image_file, weight_file, configdir, outdir, star_fwhms)
        starcat_file = get_starcat(image_file, outdir, star_fwhms)
        plotter(image_file, starcat_file, outdir)
        run_piffy(image_file, starcat_file, configdir=configdir, outdir=outdir)

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
