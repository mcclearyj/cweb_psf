import json
import numpy as np
import os
import sys
from astropy.table import Table, vstack
import pdb
from argparse import ArgumentParser
import fitsio

# Local imports
from src.plotter import compare_rho_stats
from src.star_psf_holder import StarPSFHolder
from src.plotter import make_resid_plot, plot_rho_stats
from src.quiverplot import QuiverPlot
from src.residplots import ResidPlots
from src.chisqplots import ChiSqPlots

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
        self.outdir = config['outdir']

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

    def run_psf(self, psf_model, stars):

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

    def run_all_diagnostics(self, stars, psfs):
        """ Wrapper for diagnostics plots making """

        # Main dictionary to hold all results
        all_results = {}

        # Define plot names, incl. quiver, resid and chi2
        quiv_name = os.path.join(
            self.outdir,
            '_'.join([psfs.psf_type, 'quiverplot.png'])
        )
        resid_name = os.path.join(
            self.outdir, '_'.join([psfs.psf_type, 'flux_resid.png'])
        )
        chisq_name = os.path.join(
            self.outdir, '_'.join([psfs.psf_type, 'chi2.png'])
        )

        # Make quiverplots
        quiverplot = QuiverPlot(starmaker=stars, psfmaker=psfs)
        quiverplot.run(scale=1, outname=quiv_name)

        # Save quiverplot results
        quiver_resid_vals = quiverplot.return_hsm_resid_vals()

        # Go for residplots, incl. SSIM
        resid_plot = ResidPlots(starmaker=stars, psfmaker=psfs)
        resid_plot.run(resid_name=resid_name)

        # Save resid_plot results
        flux_resid_vals = resid_plot.return_flux_resid_vals()
        ssim_resid_vals = resid_plot.return_ssim_vals()

        # Handle polydeg based on psf_type
        if psfs.psf_type == 'single':
            polydeg = 0
        else:
            polydeg = 1

        # Go for chi2plots
        chisq_plot = ChiSqPlots(starmaker=stars, psfmaker=psfs)
        chisq_plot.run(chisq_name=chisq_name, polydeg=polydeg)

        # Save chisq_plot results
        chisq_resid_vals = chisq_plot.return_chi2_resid_vals()

        # Organize results by model type
        all_results[psfs.psf_type] = {
            'HSM_resids': quiver_resid_vals,
            'ssim_resids': ssim_resid_vals,
            'flux_resids': flux_resid_vals,
            'chisq_resids': chisq_resid_vals
        }

        # Define the filename for the JSON file
        results_filename = os.path.join(
            self.outdir, f"{psfs.psf_type}_diagnostics_results.json"
        )

        # Write results to a JSON file
        with open(results_filename, 'w') as outfile:
            json.dump(all_results, outfile, indent=4)

        print(
            f"Saved diagnostics results for {psfs.psf_type}" +
            f" model to {results_filename}"
        )

    def make_ouput_subdir(self):
        """
        Bit of a kludge, but make a directory with the name/tile of
        image inside of 'outdir' for easier organization
        """
        basename = os.path.basename(self.catalog_file).split('.fits')[0]
        outdir = self.config['outdir']
        subdir = os.path.join(outdir, basename)

        # Now replace outdir with subdir!
        self.outdir = subdir

        if not os.path.isdir(subdir):
            cmd = f'mkdir -p {subdir}'
            os.system(cmd)
            print(f'Made output directory {subdir}')
        else:
            print(f'Output directory {subdir} exists, continuing...')

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
                psfs = self.run_psf(model, stars)

                # Run diagnostics
                self.run_all_diagnostics(stars, psfs)
