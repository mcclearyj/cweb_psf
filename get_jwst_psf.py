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
from diagnostics.plotter import size_mag_plot

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


def make_fwhm_tab():
    '''
    Return a look-up table containing approximate stellar
    FWHMs for each of the filters based on DEC2022 COSMOS-Web simulations
    This is horrendous coding and should be fixed ASAP -- try a yaml?
    '''

    filter_names = ['f115w','f150w','f277w', 'f444w']
    fwhms = [0.058, 0.0628, 0.125, 0.165]
    min_fwhm = [0.04, 0.04, 0.10, 0.15]
    max_fwhm = [0.085, 0.085, 0.14, 0.18]

    star_fwhms = Table([filter_names, fwhms, min_fwhm, max_fwhm],
        names=['filter_name', 'star_fwhm', 'min_fwhm', 'max_fwhm'])

    return star_fwhms


def extract_sci_wht(i2d, outdir, overwrite=False):
    '''
    Extract and save to file the "SCI" and "WHT" extensions from the input
    i2d-format image because SExtractor can't use a weight in a MEF, apparently?
    Return the science image and weight file names to be passed to SExtractor.

    Inputs:
        data_dir : directory containing images
        overwrite : if true, overwrite any SCI/WHT files saved to disk
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
                    configdir, outdir, star_fwhms):
    '''
    Run Source Extractor, including the appropriate star FWHM for
    preliminary star identification.

    Inputs
        image_file: image for SExtractor
        weight_file: weight file for SExtractor
        configdir: directory of SEx configuration files
        outdir: directory for outputs
        star_fwhms: table with appx. size of stars in image
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

    wg = np.isin(star_fwhms['filter_name'], filter_name)
    star_fwhm = star_fwhms[wg]['star_fwhm']
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
            'sex', image_file, weight_arg, name_arg, check_arg,param_arg, nnw_arg,
            filter_arg, seeing_arg, config_arg
            ])
    print("sex cmd is " + cmd)
    os.system(cmd)

    return

def make_starcat(image_file, outdir,
                truthstars=None, star_fwhms=None, thresh=0.75):
    '''
    Create a star catalog from the SExtractor image catalog using cuts on
    the SExtractor CLASS_STAR parameter and the supplied table of star properties.
    Alternatively, can match against a reference star catalog (truthstars).
    Also plot size-magnitude diagram for image and star catalogs.

    Inputs:
        image_file : image file
        out_dir : where to save star catalog (also assumed to be location of image catalog)
        truthstars : reference star catalog
        star_fwhms : table with stellar locus parameters
        thresh : CLASS_STAR threshold for accepting as star

    Outputs:
        star_cat_file: name of the star catalog (saved to file)
    '''

    # Make sure that one of either truthstars or star_fwhms is supplied
    if (truthstars is None) and (star_fwhms is None):
        print('Reference star catalog and star param table both NoneType, exiting')

    img_basename = os.path.basename(image_file)
    filter_name = re.search(r"f(\d){3}w", image_file).group() \
                    +'_'+ re.search(r"(\d){2}mas",image_file).group()

    imcat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','cat.fits')
                )
    starcat_name = os.path.join(
                outdir, img_basename.replace('sci.fits','starcat.fits')
                )
    plot_name = os.path.join(
                outdir, img_basename.replace('sci.fits','sizemag.png')
                )

    if not os.path.exists(imcat_name):
        print(f'could not find image im_cat file {cat_name}')

    else:
        im_cat = Table.read(imcat_name)

    if truthstars is not None:
        truth_star_tab = Table.read(truthstars, format='ascii')
        truthmatch = htm.Matcher(16, ra=truth_star_tab['x_or_RA'],
                                    dec=truth_star_tab['y_or_Dec'])

        cat_ind, truth_ind, dist = truthmatch.match(
                                    ra=im_cat['ALPHAWIN_J2000'],
                                    dec=im_cat['DELTAWIN_J2000'],
                                    maxmatch=1, radius = 0.5/3600)
        star_cat = cat_name[cat_ind]

    else:
        filter_name = re.search(r"f(\d){3}w", image_file).group()
        wg = np.isin(star_fwhms['filter_name'], filter_name)
        min_fwhm = star_fwhms[wg]['min_fwhm']
        min_fwhm = np.float64(min_fwhm)
        max_fwhm = star_fwhms[wg]['max_fwhm']
        max_fwhm = np.float64(max_fwhm)

        star_selec = (im_cat['CLASS_STAR'] >= thresh) \
                        & (im_cat['MAG_AUTO'] < 30) \
                        & (im_cat['FWHM_WORLD']*3600 > min_fwhm) \
                        & (im_cat['FWHM_WORLD']*3600 < max_fwhm) \
                        & (im_cat['SNR_WIN'] > 15)

        star_cat = im_cat[star_selec]

    # Save star catalog to file
    star_cat.write(starcat_name, format='fits', overwrite=True)

    # Make size-mag plot
    size_mag_plot(im_cat, star_cat, plot_name, filter_name)

    return starcat_name


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

    ra = fits.getval(im_file, 'CRVAL1')/15.0 # PIFF wants it in hours
    dec = fits.getval(im_file, 'CRVAL2')

    # Load PIFF configuration file
    run_piff_config = os.path.join(configdir, 'piff.config')

    # Now run PIFF on that image and accompanying catalog
    image_arg   = f'input.image_file_name={bkg_sub}'
    coord_arg   = f'input.ra={ra} input.dec={dec}'
    psfcat_arg  = f'input.cat_file_name={star_cat_file}'
    output_arg  = f'output.file_name={output_name} output.dir={piff_outdir}'
    cmd = ' '.join([
             'piffify', run_piff_config, image_arg, coord_arg, \
                psfcat_arg, output_arg
                ])

    print('piff cmd is ' + cmd)
    os.system(cmd)

    return

def run_psfex():

    pass

    psfcat_name = im_cat

    # Now run PSFEx on that image and accompanying catalog

    psfex_config_arg = '-c '+config_path+'psfex.mock.config'
    outcat_name = im_cat.replace('.fits','.psfex.star')
    cmd = ' '.join(
            ['psfex', psfcat_name,psfex_config_arg,'-OUTCAT_NAME', outcat_name]
            )
    self.logprint("psfex cmd is " + cmd)
    os.system(cmd)


    cleanup_cmd = ' '.join(
                    ['mv chi* resi* samp* snap* proto* *.xml', psfex_plotdir]
                    )
    cleanup_cmd2 = ' '.join(
                    ['mv count*pdf ellipticity*pdf fwhm*pdf', psfex_plotdir]
                    )
    os.system(cleanup_cmd)
    os.system(cleanup_cmd2)
# utils.run_command(cmd, logprint=self.logprint)


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
    # not necessary long-term
    star_fwhms = make_fwhm_tab()

    # Process exposures
    for i2d in i2d_images:

        print(f'Working on file {i2d}...\n\n')

        #image_file, weight_file = extract_sci_wht(i2d, outdir,
        #                                            overwrite=overwrite)

        image_file = i2d
        weight_file = i2d.replace('sci.fits', 'wht.fits')

        if crop == True:
            image_file, weight_file = crop(image_file, weight_file)

        run_sextractor(image_file, weight_file,
                                configdir, outdir, star_fwhms)

        starcat_file = make_starcat(image_file, outdir,
                                        truthstars=truthstars,
                                        star_fwhms=star_fwhms)

        run_piffy(image_file, starcat_file,
                    configdir=configdir, outdir=outdir)

    return 0


if __name__ == '__main__':

    args = parse_args()
    rc = main(args)

    if rc == 0:
        print('get_jwst_psfs.py has completed succesfully')
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
