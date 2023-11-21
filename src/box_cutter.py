import numpy as np
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import fitsio
import pdb

class BoxCutter:
    """ TO DO: ADD A CONFIG_FILE CHECKER """

    def __init__(self, config_file, image_file=None, box_size=None,
                 x=None, y=None, vignet_size=None):
        """Busted cookie-cutter."""
        self.config = config_file
        self.image_file = image_file
        self.x = x
        self.y = y
        self.vignet_size = vignet_size
        self.imcat = []
        self.image = []

        # We refer to box size a lot so read it in
        self.box_size = np.int32(self.config['box_size'])

    def _grab_wcs_box(self, obj):
        """Placeholder for what might be a good RA/Dec matcher."""
        ra_key = self.config['input_catalog']['ra_key']
        dec_key = self.config['input_catalog']['dec_key']
        ra_unit = self.config['input_catalog']['ra_unit']
        dec_unit = self.config['input_catalog']['dec_unit']

        coord = SkyCoord(ra=obj[ra_key] * ra_unit,
                         dec=obj[dec_key] * dec_unit
                         )

        wcs = WCS(fits.getheader(self.image_file))
        x, y = wcs.all_world2pix(coord.ra.value, coord.dec.value, 0)
        object_pos_in_image = [x.item(), y.item()]
        return

    def _grab_box(self, x, y):
        """Grab a region of the image.

        Parameters:
            x, y : float
                location of star from catalog
            im : array-like
                image data
            box_size : int
                vignet size to cut around star

        Returns:
            box
                Cutout box from the image.
        """
        bs = np.int(self.box_size)
        bb = self.box_size/2
        im = self.image * 1.0
        j1 = int(np.floor(x-bb))
        j2 = int(np.floor(x+bb))
        k1 = int(np.floor(y-bb))
        k2 = int(np.floor(y+bb))

        box = im[k1:k2, j1:j2]

        if np.shape(box) != (bs, bs):
            box = np.zeros([bs, bs])

        return box

    def _check_boxsize(self):
        """Check box size and update if needed."""
        if self.vignet_size is None:
            self.vignet_size = self.imcat['VIGNET'][0].shape[0]

        if self.box_size != self.vignet_size:
            print(f'\nrun_config box_size={self.box_size} and catalog',
                  f'vignet size={self.vignet_size} differ')
            print(f'Overriding run_config box_size to {self.vignet_size}\n')
            self.box_size = self.vignet_size

    def grab_boxes(self, image_file, cat_file, ext='ERR', hdu=None):
        """Load image files, call box grabber."""
        x_key = self.config['input_catalog']['x_key']
        y_key = self.config['input_catalog']['y_key']

        if isinstance(image_file, str):
            # FITS extension in which the image to cut out is located
            if hdu == None:
                hdu = self.config[f'{ext.lower()}_image']['hdu']
            else:
                print("reading supplied hdu") # Probably using an HDU supplied as an argument
            print(f'reading image extension {ext} at hdu {hdu}')
            imf = fitsio.FITS(image_file, 'r')[hdu]
            self.image = imf.read()
        else:
            self.image = image_file

        if isinstance(cat_file, str):
            imcat_fits = fitsio.FITS(cat_file, 'rw')
            cat_hdu = self.config['input_catalog']['hdu']
            self.imcat = imcat_fits[cat_hdu].read()
        else:
            self.imcat = cat_file

        self._check_boxsize()

        x = self.imcat[x_key]
        y = self.imcat[y_key]
        new_boxes = list(map(self._grab_box, x, y))

        return new_boxes

    def add_boxes_to_imcat(self, imcat, new_boxes, ext):
        """Convenience function that adds new_boxes to imcat as
        a column called {ext}_VIGNET. Should finish this later"""
        pass
