import numpy as np
from numpy.lib import recfunctions as rf
from pathlib import Path
from glob import glob
import fitsio
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import galsim as gs
import astropy.units as u
from astropy.table import Table, vstack
import glob
from time import time
import copy

from .config import CookieCutterConfig
import ipdb

class CookieCutter():

    def __init__(self, config=None):

        if config is not None:
            utils.check_type('config', config, (str, Path, dict))

        if not isinstance(config, dict):
            config = utils.read_yaml(str(config))

            # also sets defaults for optional params
            self.config = CookieCutterConfig(config)
        else:
            self.config = None


    def initialize(self):
        '''
        Initalize CookieCutter from either a config or file. Actually
        loads the catalog or stamp data into memory
        '''

        config = self.config

        # Parse the config file.
        # Decide if we're initializing a cookie cutter object from an existing
        # cookie cutout file or if we're building one from scratch.

        # Initialize from config
        input_dir = config['input']['dir']

        catalog_file = Path(config['input']['catalog'])

            if input_dir is not None:
                catalog_file = input_dir / catalog_file

            self.input_dir = input_dir

            ext = config['input']['catalog_ext']
            catalog = fitsio.read(str(catalog_file), ext=ext)

            self.catalog_file = catalog_file
            self.catalog = catalog

            self.wcs_type = config['input']['wcs_type']


        else:
            self._fits = fitsio.FITS(cc_file, 'rw')
            return

        ret

    def compute_image_pos(self, image, obj, wcs=None):
        '''
        Compute the position of an object in pixel coordinates
        given an image with a registered WCS
        NOTE: This method is quite slow, but uses the typical high-level
        astropy WCS interface
        image: fitsio.hdu.image.ImageHDU
            A FITS Image HDU object
        obj: astropy.table.Row
            A row of an astropy table
        wcs: astropy.WCS
            An astropy WCS instance, if you want to use one other than
            the SCI image (such as the coadd segmap)
        '''

        id_tag = self.config['input']['id_tag']
        ra_tag = self.config['input']['ra_tag']
        dec_tag = self.config['input']['dec_tag']

        ra_unit = u.Unit(self.config['input']['ra_unit'])
        dec_unit = u.Unit(self.config['input']['dec_unit'])


        coord = SkyCoord(
            ra=obj[ra_tag]*ra_unit, dec=obj[dec_tag]*dec_unit
            )

        if wcs is None:
            # default behavior is to grab it from the passed image header
            wcs = WCS(fits.getheader(imfile))

        image_shape = image.get_info()['dims']

        # x, y = wcs.world_to_pixel(coord)
        x, y = wcs.all_world2pix(
            coord.ra.value, coord.dec.value, 0
            )

        y, x = wcs.all_world2pix(
                ra, dec, 0, ra_dec_order=True
                )
        object_pos_in_image = [x.item(), y.item()]

        # NOTE: reversed as numpy arrays have opposite convention!
        object_pos_in_image_array = object_pos_in_image[-1::-1]

        return object_pos_in_image_array
