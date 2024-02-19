import psfex
import galsim
import treecorr
import piff
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from astropy.table import Table, vstack
import pdb
from argparse import ArgumentParser
import yaml
#from astropy.io import fits
import fitsio
from src.plotter import compare_rho_stats
from src.star_psf_holder import StarPSFHolder
from src.plotter import make_resid_plot, plot_rho_stats
from src.quiverplot import QuiverPlot
from src.residplots import ResidPlots


class PSFDiagnosticsRunner:
    """
    Load in the catalog, calculate star stamp sky bkg stats, and
    and for each of the PSF types supplied in the config, create a bunch of
    diagnostic images. You know what? I was right the first time.
    """

    def __init__(self, catalog_file, config, psf_type=None, vignet_size=None):
        self.catalog_file = catalog_file
        self.config = config # Config should have been read in with a runner script
        self.vignet_size = vignet_size
        self.stars = []
        self.psfs = []
        self.sky_level = []
        self.sky_std = []
        self.star_fluxes = []

        # Load FITS catalog with stars and everything
        try:
            self.catfitsObj = fitsio.FITS(catalog_file, 'rw')
            hdu = config['input_catalog']['hdu']
            self.catalog = self.catfitsObj[hdu].read()
        except FileNotFoundError as fnf:
            print('PSFDiagnostics: no catalog found at supplied filepath:')
            print(fnf)

        # Allowed PSF models
        self.model_map = [
            'psfex',
            'webbpsf',
            'single',
            'mirage',
            'piff',
        ]

    def load_star_cat(self, star_holder):
        """
        Set stars attributes: VIGNET, x, y, RA, dec. Maybe this should be
        moved to StarPSFHolder?
        """
        # For convenience
        starcat_params = self.config['input_catalog']

        star_holder.stamps = self.catalog['VIGNET']
        star_holder.x = self.catalog[starcat_params['psf_x_key']]
        star_holder.y = self.catalog[starcat_params['psf_y_key']]
        star_holder.ra = self.catalog[starcat_params['ra_key']]
        star_holder.dec = self.catalog[starcat_params['dec_key']]
        star_holder.star_fluxes = self.catalog[starcat_params['flux_key']]
        if 'err_image_key' in starcat_params.keys():
            star_holder.err_stamps = self.catalog[
                starcat_params['err_image_key']
            ]
        # Non-catalog params
        star_holder.pixel_scale = self.config['pixel_scale']
        star_holder.vignet_size = self.config['box_size']
        print(f'PSFmaker vignet size is {self.vignet_size}')

    def _run_all_diagnostics(self, stars, psfs):
        """ Wrapper for diagnostics plots making """

        # Make output quiverplot
        quiv_name = os.path.join(
            self.config['outdir'],
            '_'.join([psfs.psf_type,'quiverplot.png'])
        )

        # Go.
        quiverplot = QuiverPlot(starmaker=stars, psfmaker=psfs)
        quiverplot.run(scale=1, outname=quiv_name)

        # Define resid plot names, incl. chi2
        resid_name = os.path.join(
            self.config['outdir'], '_'.join([psfs.psf_type,'flux_resid.png'])
        )
        chi2_name = os.path.join(
            self.config['outdir'], '_'.join([psfs.psf_type,'chi2.png'])
        )

        # Go.
        resid_plot = ResidPlots(starmaker=stars, psfmaker=psfs)
        resid_plot.run(resid_name=resid_name, chi2_name=chi2_name)


    def _run_psf(self, psf_model, stars):

        psfs = StarPSFHolder(psf_type=psf_model)
        psfs.copy_from_holder(stars)

        # Load star stamps
        psfs.load_psf_stamps(catalog=self.catalog)

        # Trim stamps. Maybe one day we can import as a function?
        psfs.trim_stamps(stars)

        # Add flux!
        psfs.add_flux()

        # Do HSM fits
        psfs.do_hsm_fit()

        return psfs

    def run_all(self, vb=False):
        """
        Loop over all PSF types and get diagnostics, creating an instance of
        CalcPSFDiagnostics for each type. I wonder whether it makes more sense
        to do the looping over all PSF types here or in the runner script.
        """
        stars = StarPSFHolder(psf_type='star', vb=vb)

        # If someone didn't already call it...
        self.load_star_cat(stars)

        # Calculate sky background
        stars.calc_sky_bkg()

        # Get star HSM fits
        stars.do_hsm_fit()

        # Loop over allowed model types
        for model in self.model_map:

            # If model is "True" in config, go!
            if self.config['psf_models'].get(model, False):
                # Friendly alert
                print(f'\n Running fits for PSF model type: {model.upper()}\n')

                # Initialize PSF object
                psfs = self._run_psf(model, stars)

                # Run diagnostics
                self._run_all_diagnostics(stars, psfs)
