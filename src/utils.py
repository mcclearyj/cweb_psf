import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
from esutil import htm
import astropy.wcs as wcs
import yaml
import re
import os
import glob
from astropy.table import Table, vstack
import pdb
import functools
import numpy as np

class AttrDict(dict):
    '''
    More convenient to access dict keys with dict.key than dict['key'],
    so cast the input dict into a class!
    '''

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def set_rc_params(fontsize=None):
    '''
    Set figure parameters
    '''

    if fontsize is None:
        fontsize=16
    else:
        fontsize=int(fontsize)

    rc('font',**{'family':'serif'})
    rc('text', usetex=True)

    #plt.rcParams.update({'figure.facecolor':'w'})
    plt.rcParams.update({'axes.linewidth': 1.3})
    plt.rcParams.update({'xtick.labelsize': fontsize})
    plt.rcParams.update({'ytick.labelsize': fontsize})
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
    plt.rcParams.update({'axes.labelsize': fontsize})
    plt.rcParams.update({'axes.titlesize': fontsize})
    plt.rcParams.update({'legend.fontsize': int(fontsize-2)})

    return

def match_coords(cat1, cat2, ratag1=None, dectag1=None, ratag2=None, dectag2=None, radius=0.5):
    '''
    Utility function to match cat1 to cat 2 using celestial coordinates
    '''

    # Either 'ra/dec' or 'ALPHAWIN_J2000/DELTAWIN_J2000'!

    try:
        if (ratag1 is not None) and (dectag1 is not None):
            cat1_ra = cat1[ratag1]
            cat1_dec =  cat1[dectag1]
        elif 'ra' in cat1.colnames:
            cat1_ra = cat1['ra']
            cat1_dec =  cat1['dec']
        elif 'ALPHAWIN_J2000' in cat1.colnames:
            cat1_ra = cat1['ALPHAWIN_J2000']
            cat1_dec =  cat1['DELTAWIN_J2000']
        else:
            raise KeyError('cat1: no "ra,dec" or "{ALPHA,DELTA}WIN_J2000" columns')
    except:
        print("Couldn't load catalog 1 RA & Dec")

    try:
        if (ratag2 is not None) and (dectag2 is not None):
            cat2_ra = cat2[ratag1]
            cat2_dec =  cat2[dectag1]
        elif 'ra' in cat2.colnames:
            cat2_ra = cat2['ra']
            cat2_dec =  cat2['dec']
        elif 'ALPHAWIN_J2000' in cat2.colnames:
            cat2_ra = cat2['ALPHAWIN_J2000']
            cat2_dec =  cat2['DELTAWIN_J2000']
        else:
            raise KeyError('cat2: no "ra,dec" or "{ALPHA,DELTA}WIN_J2000" columns')
    except:
        print("Couldn't load catalog 2 RA & Dec")

    cat1_matcher = htm.Matcher(16, ra=cat1_ra, dec=cat1_dec)

    cat2_ind, cat1_ind, dist = cat1_matcher.match(
        ra=cat2_ra, dec=cat2_dec, maxmatch=1, radius=radius/3600.
    )

    print(f'{len(dist)}/{len(cat1)} gals matched to truth')

    return cat1[cat1_ind], cat2[cat2_ind]

def read_yaml(yaml_file):
    '''
    current package has a problem reading scientific notation as
    floats; see
    https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    '''

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


def make_outdir(config, arg=None, path='./', cval='outdir'):
    '''
    Make an output directory. Once we finally move to OOP, we won't have
    to return the config like it's 2003.

    Inputs:
        arg: the full outdir path. If it doesn't exist, one will be made.
        path: the base path in which outdir should be saved
              [default: './']
        cval: the config parameter specifying output directory
              [default: 'outdir']
    Outputs:
        config: A parameter named 'outdir' is added if it was missing
    '''

    if arg == None:
        outdir = config[cval]
    else:
        outdir = arg
        config[cval] = outdir
    # Set a sensible default if outdir is still none
    if outdir == None:
        #basedir = os.path.commonpath(images)
        basedir = path
        outdir = os.path.join(basedir,'working')
        config['outdir'] = os.path.join(basedir,'working')

    if not os.path.isdir(outdir):
        cmd = 'mkdir -p {outdir}'.format(outdir=outdir)
        os.system(cmd)
        print(f'Made output directory {outdir}')
    else:
        print(f'Output directory {outdir} exists, continuing...')

    return config

def concatenate_catalogs(catnames, outname=None, hdu=1, pad_size=None):
    """
    This I want to run on all the teeny single exposure catalogs.
    I am going to operate like I have the run config available, too.
    It might make sense to trim the tables to the value in the run_config...
    but since the catalog is being saved to file, maybe not?
    Another option might be to specify an output stamp size and pad or not?
    """

    if outname == None:
        outname='combined_catalog.fits'

    if type(catnames) == str:
        fits_files = glob.glob(catnames)
    elif type(catnames) == list:
        fits_files = catnames
    else:
        raise ValueError("concatenate_catalogs: catnames must be string or list")

    print(f"utils.concatenate_catalogs: joining tables {catnames}")

    # List to hold tables
    table_list = []

    # Loop over all your FITS files
    for fits_file in fits_files:
        table = Table.read(fits_file, hdu=hdu)
        # Pad stamps if requested...
        if pad_size:
            table = pad_stamps(table, pad_size)
        table_list.append(table)

    # Vertically stack tables
    combined_table = vstack(table_list)

    # Save the combined table to a new FITS file
    combined_table.write(outname, format='fits', overwrite=True)

    return combined_table

def pad_stamps(table, pad_size):
    """
    Function to pad stamps to the same size. Eventually, it would be good to
    be able to specify extensions, refer to run_config, or something like that.
    Maybe even combine it with a stamp trimming operation? So standardize all
    stamps to be the same size specified in run_config? It would probably be
    the simplest, especially if I combine catalogs from different runs.

    Inputs
      table: astropy Table instance
      pad_size: desired output size of stamps

    Returns: updated Table instance
    """

    # Access columns with "VIGNET" in the name that need to be padded.
    # There has GOT to be a way to clean up...
    colnames = table.colnames
    pad_colnames = []
    for col in colnames:
        if "VIGNET" in col:
            if table[col][0].shape[0] < pad_size:
                pad_colnames.append(col)

    # Pad columns that we need to pad
    for col in pad_colnames:
        unpadded = table[col].data
        padded = list(
            map(functools.partial(_pad_stamps, pad_size), table[col].data)
        )
        table[col] = np.array(padded)

    return table

def _pad_stamps(pad_size, stamp):
    """Function that gets called to actually pad stamps based on pad_size"""
    npix = int((pad_size - stamp.shape[0])/2)

    # Make sure stamp is actually smaller than desired padding size
    if npix < 0:
        raise ValueError(
            f"pad_stamps: pad size {pad_size} must be larger than stamp size {stamp.shape[0]}!!"
        )

    # Do the actual padding
    padded_stamp = np.pad(stamp, pad_width=npix)

    # Make sure output stamp is the right size!
    if padded_stamp.shape[0] != pad_size:
        raise ValueError(
            f"Output padded stamp size {padded_stamp.shape[0]} not desired stamp size of {pad_size}!"
        )
    return padded_stamp
