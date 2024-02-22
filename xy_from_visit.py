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
import pdb

# Local imports
from max_polygons import which_visit

class GrabXYCoords:
    """
    Utility to grab resampled single-visit mosaic XY coordinates from RA, Dec.

    TO DO
        - Incorporate a loop over desired bandpasses
        - Find a way to avoid looping over coordfile names 380,000 times
    """

    def __init__(self, catalog, visit_mosaics, config):
        self.catalog = catalog
        self.imcat = catalog[config['input_catalog']['hdu']].data
        self.ras = self.imcat[config['input_catalog']['ra_colname']]
        self.decs = self.imcat[config['input_catalog']['dec_colname']]
        self.visit_mosaics = visit_mosaics
        self.config = config
        self.output_cat = None

    def grab_coords(self, imcat, config):
        """
        Load in coordinates from catalog, return coord tuple. This is in here
        just in case it becomes useful later
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

    def convert_wcs2pix(self, imfile, ras, decs):
        """
        Use JWST ImageModel to convert the RA, Dec of galaxies that lie in the
        footprint of the given single-visit mosaic, to pixel X, Y. Function assumes
        that

        Inputs
            imfile: image file to use for WCS
            ras, decs: astrometric coordinates

        Returns
            x, y: pixel coordinates
        """

        # Get header & build a WCS from it
        with ImageModel(imfile) as model:
            radec_to_xy = model.meta.wcs.get_transform('world', 'detector')
            x, y = radec_to_xy(ra, dec)

        # Return imcat_fits, now with corrected RA, Dec column
        return x, y

    def create_output_catalog(self, visits):
        """
        Create a dummy FITS catalog, which will be populated
        with RA, Dec, size measures from real catalog and X/Y from visit
        """
        imcat = self.imcat[10000:12000]
        ras = self.ras[10000:12000]
        decs = self.decs[10000:12000]
        number = self.imcat['NUMBER'][10000:12000]
        visits = np.array(visits, dtype=int)

        # These two will hold the pixel X, Y
        dummy_x = np.ones(len(imcat))
        dummy_y = np.ones(len(imcat))

        tab = Table()
        tab.add_columns(
            [number, visits, ras, decs, dummy_x, dummy_y],
            names=['number', 'visit_num', 'ra', 'dec', 'visit_X', 'visit_Y']
        )
        '''
        # Add these back in when I have real catalog
        tab.add_columns([
            imcat['Id'],
            imcat['RA_DETEC'], imcat['DEC_DETEC'],
            imcat['RA_MODEL'], imcat['DEC_MODEL'],
            imcat['AREA'], imcat['RADIUS'], imcat['RADIUS_err'],
            imcat['AXRATIO'], imcat['AXRATIO_err'],
            imcat['E1'], imcat['E1_err'],
            imcat['E2'], imcat['E2_err']
            ]
        )
        '''

        self.output_cat = tab

    def assign_visit_number(self, coords_file='coords_f150w_CW_JAN2024.txt'):
        """
        Invoke which_visit() to assign a particular visit to every RA, Dec

        TO DO: loop over bandpasses!
        could rename visits as: bandpass_visit = [i + '_f277w' for i in visits]
        """
        ras = self.ras
        decs = self.decs

        return which_visit(coords_file, ras, decs)

    def get_pixcoord_from_visit(self):
        """
        Couldn't think of a better name, but this will do the looping
        Saving this snippet in case I need it
        #unique_visits = [int(i) for i in np.unique(output_cat['visit_num'])]
        """

        # For convenience
        output_cat = self.output_cat

        # Loop over visits, reading in image WCS one at a time
        for visit_mosaic in self.visit_mosaics:

            # Friendly alert
            print(f"Working on mosaic {visit_mosaic}")

            # Find entries in catalog in this visit
            this_visit = int(os.path.basename(mosaic_name).split('_')[1])
            wg = output_cat['visit_num'] == this_visit

            # Call to transform coords
            x, y = self.convert_wcs2pix(
                visit_mosaic, self.ras[wg], self.decs[wg]
            )

            # Replace coordinates
            output_cat['visit_X'][wg] = x
            output_cat['visit_Y'][wg] = y

        return

    def save_catalog(self):
        """ Save augmented catalog to file """

        outcatpath = os.path.join(
            self.config['output_catalog']['path'],
            self.config['output_catalog']['name']
        )

        self.output_cat.write(outcatpath, format='fits', overwrite=True)

    def run(self):
        """ Perform all functions """

        # Get corresponding visit numbers for every catalog entry
        visits = self.assign_visit_number()

        # Create dummy output catalog
        self.create_output_catalog(visits)

        # Save to file as a pseudocheckpoint
        self.save_catalog()

        # Loop over visit visits, select matching gals, get X, Y
        self.get_pixcoord_from_visit()

        # Save to file again
        self.save_catalog()
