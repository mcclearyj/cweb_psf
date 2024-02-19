import numpy as np
import os
from astropy.io import fits
import pdb
from astropy.table import Table, Column
import glob
import psfex
import fitsio
from src.box_cutter import BoxCutter
import re
import piff

class PSFRenderer:
    def __init__(self, image_file, cat_file, run_config):
        self.image_file = image_file # Just there to set PSF file names
        self.run_config = run_config

        # Load the catalog file
        self.gc_fits = fitsio.FITS(cat_file, 'rw')
        self.hdu = run_config['input_catalog']['hdu']
        self.gal_catalog = self.gc_fits[self.hdu].read()

        # Set x, y values of objects
        if run_config['input_catalog'].get('psf_x_key') != None:
            self.x = self.gal_catalog[run_config['input_catalog']['psf_x_key']]
            self.y = self.gal_catalog[run_config['input_catalog']['psf_y_key']]
        else:
            self.x = self.gal_catalog[run_config['input_catalog']['x_key']]
            self.y = self.gal_catalog[run_config['input_catalog']['y_key']]

        # Set the baseline image/vignet size
        self.vignet_size = self.gal_catalog['VIGNET'][0].shape[0]

        # Map models to their respective methods
        self.model_map = {
            'psfex': self._render_psfex,
            'webbpsf': self._render_webbpsf,
            'single': self._render_single,
            'mirage': self._grab_mirage,
            'piff': self._render_piff
        }

    def _trim_box(self, psf_im):
        if psf_im.shape[0] > self.vignet_size:
            n = int((psf_im.shape[0]-self.vignet_size)/2)
            k = n; j = -n
        else:
            k = 0; j = None

        psf_im = psf_im[k:j,k:j]

    def _add_to_cat(self, psfims, ext):
        """Convenience function to add PSF renderings to imcat.
        Inputs
            boxes: list of PSF renderings
            ext: name of PSF model (e.g., PSFEX, WEBBPSF, ...)
        """
        # Box size
        bs = psfims[0].shape[0]

        # Dummy column to hold PSF data
        data = np.zeros(len(psfims), dtype=[(f'{ext}_VIGNET',
                        'f4', (bs, bs))])

        for i, psfim in enumerate(psfims):
            data[f'{ext}_VIGNET'][i] = psfim

        # Add column & print comfort message
        self.gc_fits[self.hdu].insert_column(f'{ext}_VIGNET',
                                             data[f'{ext}_VIGNET'])
        print(f'\nAdded {ext} vignets to {self.gc_fits._filename}\n')

    def _render_psfex(self):
        """ Render PSFEx image. """
        print("Rendering PSFEx")

        # Preliminaries: find PSF, define output directory
        basename = os.path.basename(self.image_file)
        psf_name = basename.replace('.fits', '_starcat.psf')
        psfex_dir = os.path.join(
                    self.run_config['outdir'],
                    'psfex-output', basename.split('.')[0]
                    )
        psf_path = os.path.join(psfex_dir, psf_name)

        # Load in the PSF
        psf = psfex.PSFEx(psf_path)

        # Create empty list that will hold the PSF renderings
        psfex_images = []

        # Loop prevents memory overflow errors for very large catalogs
        for xi, yi in zip(self.x, self.y):
            psfex_images.append(psf.get_rec(yi, xi))

        # Add psf images to catalog
        self._add_to_cat(psfex_images, ext="PSFEX")

    def _render_webbpsf(self):
        # Code for rendering WebbPSF goes here
        pass

    def _render_single(self):
        # Render Marko-format PSF (no spatial variation)

        psfdir = os.path.dirname(self.image_file)
        bandpass = re.search(r"f(\d){3}w", self.image_file).group()
        tile = re.search(r"[A,B](\d){1}", self.image_file).group()
        #psfname = f'JWST-PSFEx_out_{bandpass}_{tile}_psf_v4.0.psf'
        psfname = '/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2024/marko_psfex/JWST-PSFEx_out_f277w_B4_psf_v4.0.psf'

        # Open FITS-format PSF, grab extensions
        single_psf = fits.open(psfname)
        extension_names = [ext.name for ext in single_psf]

        if 'DET_DIST' in extension_names:
            ext = 'DET_DIST'
        elif 'PSF_DATA' in extension_names:
            ext = 'PSF_DATA'
        else:
            raise "Extension not found in WebbPSF/SinglePSF model"

        try:
            # Marko's PSFEx format
            psf_im = single_psf[ext].data['PSF_MASK'][0, 0, :,:]
        except:
            # Probably WebbPSF model
            psf_im = single_psf[ext].data

        # Since it is the same for all catalog entries, create a list as long
        # as the catalog entries and populate it with images.
        single_psf_images = []

        # Loop prevents memory overflow errors for very large catalogs
        for i in range(len(self.x)):
            single_psf_images.append(psf_im)

        # Add psf images to catalog
        self._add_to_cat(single_psf_images, ext="SINGLE")


    def _locate_mirage_file(self):
        """
        Utility method to extract the path to a MIRAGE point source image for
        simulated data. Returns the absolute path of the point source image.
        """

        basename = os.path.basename(self.image_file)
        dirname = os.path.dirname(self.image_file)
        try:
            # This saves a few seconds, but will break for different file orgs
            filter_name = re.search(r"f(\d){3}w", dirname).group()
        except:
            filter_name=''

        mirage_path = os.path.join(self.run_config.get('mirage_dir'),
                      filter_name,
                      basename.replace('cal.fits', \
                      '*CLEAR_ptsrc_seed_image.fits')
                      )

        return mirage_path

    def _grab_mirage(self, ext='MIRAGE'):
        """ This is specific to single-exposure sims, as noiseless,
        perfect point source images aren't usually available for real data!
        Note: MIRAGE PSF models (aka point source seed images) are assumed to
        live in the image directory"""

        if self.run_config.get('mirage_dir') != None:
            mirage_path = self._locate_mirage_file()
        else:
            mirage_path = self.image_file.replace('cal.fits',
                          '*CLEAR_ptsrc_seed_image.fits')

        mirage_file = glob.glob(mirage_path)[0]

        # Create "image_file" for boxcutter
        try:
            mirf = fits.open(mirage_file)
            mirage_im = mirf['DATA'].data
        except FileNotFoundError as fnf:
            print(
            f'Could not find a MIRAGE point source image at {mirage_name}', fnf
            )

        # Start by create a BoxCutter instance
        boxcut = BoxCutter(config_file=self.run_config,
                           vignet_size=self.vignet_size)

        # Call to grab_boxes method
        psf_images = boxcut.grab_boxes(image_file=mirage_im,
                                       cat_file=self.gal_catalog)

        # Add psf images to catalog
        self._add_to_cat(psf_images, ext="MIRAGE")

    def _render_piff(self):
        """ Render PIFF image. """

        # Preliminaries: find PSF, define output directory
        print("Rendering PIFF")
        basename = os.path.basename(self.image_file)
        psf_name = basename.replace('.fits', '.piff')
        piff_dir = os.path.join(self.run_config['outdir'],
                                'piff-output',
                                basename.split('.')[0])
        psf_path = os.path.join(piff_dir, psf_name)

        # Load in the PSF
        piff_psf = piff.read(psf_path)

        # Create empty list that will hold the PSF renderings
        piff_images = []

        # Loop prevents memory overflow errors for very large catalogs
        for xi, yi in zip(self.x, self.y):
            piff_im = piff_psf.draw(x=xi, y=yi, stamp_size=self.vignet_size)
            piff_images.append(piff_im.array)

        # Add psf images to catalog
        self._add_to_cat(piff_images, ext="PIFF")

    def render(self):
        for model, render_method in self.model_map.items():
            if self.run_config['psf_models'].get(model, False):
                try:
                    render_method()
                except:
                    pdb.set_trace()
        self.gc_fits.close()

# Usage
#renderer = PSFRenderer(image_file, cat_file, run_config)
#renderer.render()
