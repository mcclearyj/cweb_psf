import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import os, re
import glob
import yaml
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy import units as u
import jwst
from jwst.datamodels import ImageModel


def read_yaml(yaml_file):
    """
    Load a YAML-format configuration file 
    current package has a problem reading scientific notation as
    floats; see
    https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    """

    loader = yaml.SafeLoader
    loader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
        [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(yaml_file, 'r') as stream:
        # return yaml.safe_load(stream) # see above issue
        return yaml.load(stream, Loader=loader)

def grab_coords(imcat, config):
    """
    Load in coordinates from catalog, return SkyCoord object.
    Inputs
        imcat: FITS instance of catalog
        config: main run configuration dict
    """
    
    # Load in catalog
    cat_hdu = config['input_catalog']['hdu']
    imcat = imcat_fits[cat_hdu].read()

    # Grab RA, Dec
    ra_key = config['input_catalog']['ra_key']
    dec_key = config['input_catalog']['dec_key']
    imcat_ra = imcat[ra_key]; imcat_dec = imcat[dec_key]
    
    return (imcat_ra, imcat_dec)

def convert_wcs2pix(imfile, coords, config):
    """
    Use astropy.wcs to pick out objects that lie in the footprint of a given single-visit
    mosaic, then use its WCS to convert RA, Dec to pixel X, Y. 
    Inputs
        imfile: image file to use for WCS
        imcat_fits: fitsio.FITS instance of catalog
        run_config: main run configuration dict (as opposed to star config)
    """

    # Get header & build a WCS from it
    with ImageModel(imfile) as model:
        radec_to_xy = model.meta.wcs.get_transform('world', 'detector')
        x, y = radec_to_xy(ra, dec)

    # Not entirely sure how one might do this; polygons? Assuming images
    # are all the same dimensions.

    # CRPIX are 1 in this WCS convention
    x_cen = 9163.4587; y_cen = 3783.8069
    x_min = x_cen - (17717.755/2); x_max = x_cen + (17717.755/2)
    y_min = y_cen - (6958.1233/2); y_max = y_cen + (6958.1233/2) 
    
    wg = (x > x_min) & (x<x_max) & (y > y_min) & (y < y_max)
    
    
    # Add corrected RA & Dec to catalog
    print(f'Adding corrected RA/Dec columns to catalog')
    imcat_fits[cat_hdu].insert_column('x_visit', coords[:, 0])
    imcat_fits[cat_hdu].insert_column('y_visit', coords[:, 1])

    # Return imcat_fits, now with corrected RA, Dec column

    
    return imcat_fits


def create_catalog(imcat_fits):
    """
    Create a dummy FITS catalog, which will be populated
    with RA, Dec, size measures from real catalog and X/Y from visit
    """
    
    tab = Table()
    table.add_columns([
        imcat['Id'],
        imcat['RA_DETEC'], imcat['DEC_DETEC'],
        imcat['RA_MODEL'], imcat['DEC_MODEL'],
        imcat['AREA'], imcat['RADIUS'], imcat['RADIUS_err'],
        imcat['AXRATIO'], imcat['AXRATIO_err'],
        imcat['E1'], imcat['E1_err'],
        imcat['E2'], imcat['E2_err']
        ]
    )

    return tab

def main()

    # Load in config file
    config = read_yaml('./xy_from_visit_config.yaml')

    # Change to config value ASAP
    catalog_path = '/Users/j.mccleary/Research/jwst_cosmos/real_data/catalogs/'
    catalog_name = 'COSMOSWeb_master_v1.6.0-sersic+BD-em_cgs_LePhare_nodupl.fits'
    catalog = fits.open(os.path.join(catalog_path, catalog_name))
    
    image_path = '/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2024/20mas_resamp'
    visit_mosaics = glob.glob(
        os.path.join(image_path, '*fits')
    )

    # Load in catalog, retrieve coords
    coords = grab_coords(image_path, config)

    # Loop over visit mosaics, pick out overlapping objects, 

    

    

