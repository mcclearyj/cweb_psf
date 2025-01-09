import numpy as np
import os, re
from astropy.io import fits
import pdb
from astropy.table import Table, vstack, hstack
import glob
from esutil import htm
from argparse import ArgumentParser
import ipdb
from src.plotter import size_mag_plots
from src.utils import read_yaml
from src.box_cutter import BoxCutter
import fitsio

class PSFRunner:
    '''
    Do things like hold configuration file, list of exposures, etc.
    '''

    def __init__(self, config_file):
        self.exposure_list = []
        self.config_file = config_file
        self.config = {}
        self.outdir = None

        self._load_config()

    def _load_config(self):
        '''
        Load up a supplied config file. TO DO: add in error-checking;
        if configdir, outdir aren't defined in config, add them here
        '''
        self.config = read_yaml(self.config_file)
        self.configdir = config['configdir']
        self.outdir = config['outdir']

        # Load in a truth star catalog or don't, idc
        try:
            truthstars = config['truthstars']
            if ['none', 'None'] in truthstars:
                truthstars = None
        except KeyError:
            truthstars = None

        # Set default output directory values if none provided
        if configdir is None:
            configdir = '/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config/'

        if outdir is None:
            basedir = os.path.commonpath(images)
            outdir = os.path.join(basedir,'working')

        if not os.path.isdir(outdir):
            cmd = 'mkdir -p {outdir}'.format(outdir=outdir)
            os.system(cmd)
            print(f'Made output directory {outdir}')

        else:
            print(f'Output directory {outdir} exists, continuing...')


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
