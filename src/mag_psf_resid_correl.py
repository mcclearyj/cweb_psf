import numpy as np
import re, os
import pdb
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
import yaml
import time

from .utils import set_rc_params, read_yaml

class MagPsfResidCorrel:
    """
    This class is intended to compute PSF size and ellipticity residuals in
    bins of stellar magnitude and then save to file. Code is structured such
    that it can be run either inside some larger cweb_psf diagnostic run, i.e.,
    run_psf_diagnostics.py, OR independently on a cweb_psf output HSM fit table.

    The inputs are (1) a catalog with star positions and fluxes plus HSM moment
    fits for stars and PSF models and (2) a standard cweb_psf run config

    After reading in data and config, next step is to bin data by magnitude
    and compute mean and standard deviation of star-PSF size and ellipticity
    residuals. Finally, save results as a plot and (maybe) a file for later
    comparison.

    The code assumes that the input table is a cweb_psf HSM output file. Error
    checking should be added at some point.

    TO DO:
        - Add error checking for input file.
        - Save to a table?
        - Do a comparison figure
        - Improve documentation

    """

    def __init__(self, hsm_cat, run_config, model, zeropoint=28.086519392, nbins=12):

        t_start = time.perf_counter()
        # Catalog!
        self.hsm_cat = None
        # Cweb_psf-format run configuration file
        self.run_config = run_config
        # PSF model to use
        self.model = model
        # Zeropoint for magnitude calculation
        self.zeropoint = zeropoint
        # How many bins for mag histogram
        self.nbins = nbins
        # Holds magnitude histogram
        self.mag_bin_hist = []
        # Holds magnitude bin index for each star
        self.mag_bin_numbers = []

        # First, set configuration file atttribute.
        if type(run_config) == str:
            self.run_config = read_yaml(run_config)
        else:
            self.run_config = run_config

        # Read in catalog
        if type(hsm_cat) == astropy.table.table.Table:
            self.hsm_cat = hsm_cat
        else:
            try:
                self.hsm_cat = Table.read(
                    hsm_cat, hdu=run_config["input_catalog"]["hdu"]
                )
            except FileNotFoundError as fnf:
                #print(f"File {hsm_cat} not found")
                raise(fnf)

        t_stop = time.perf_counter()
        print(f"init elapsed time: {t_stop - t_start}")

    def bin_by_magnitude(self):
        """
        Convert fluxes to magnitudes, create magnitude bins, digitize
        accordingly.
        """

        # First convert fluxes to mags
        mags = self.zeropoint -2.5 * np.log10(self.hsm_cat['flux_auto'])

        # Create magnitude bins & including number in each, then digitize.
        # Pad rightmost edge so that all points get included in digitization
        mag_bin_hist = np.histogram(
            mags, bins=self.nbins, range=[np.min(mags), np.max(mags)+0.1]
        )
        bin_numbers = np.digitize(mags, mag_bin_hist[1])

        # Set class attributes and go
        self.mag_bin_hist = mag_bin_hist
        self.mag_bin_numbers = bin_numbers

    def calculate_statistics(self):
        """
        Calculate mean and standard deviation of e1, e2, T and FWHM residuals
        Do this for each PSF model in configuration file!
        """

        # Initialize dict that will hold the residuals for each model
        #full_resid_dict = {}

        # Set up star stats
        star_fwhm = self.hsm_cat["star_fwhm"]
        star_g1 = self.hsm_cat["star_hsm_g1"]
        star_g2 = self.hsm_cat["star_hsm_g2"]
        star_T  = 2.0*(self.hsm_cat["star_hsm_sig"]**2)

        # Now PSF model stats; this will eventually be a loop over
        # all models defined in config but let's start simple
        psf_fwhm = self.hsm_cat[f"{self.model}_fwhm"]
        psf_g1 = self.hsm_cat[f"{self.model}_hsm_g1"]
        psf_g2 = self.hsm_cat[f"{self.model}_hsm_g2"]
        psf_T  = 2.0*(self.hsm_cat[f"{self.model}_hsm_sig"]**2)

        # Compute residuals
        dfwhm = star_fwhm - psf_fwhm
        dg1 = star_g1 - psf_g1
        dg2 = star_g2 - psf_g2
        dT = star_T - psf_T; dTT = dT/star_T

        # Initialize dict that will hold this model's residuals
        mag_bins = np.unique(self.mag_bin_numbers)

        model_resid = {
            "mag_bin_number": [], "bin_count": [],
            "mean_dg1": [], "std_dg1": [], "mean_dg2": [], "std_dg2": [],
            "mean_dT": [], "std_dT": [], "mean_dfwhm": [], "std_dfwhm": [],
        }

        # Now do loop over bin numbers; there is probably a clever way to
        # vectorize this, but we'll do a loop for now.
        for mag_bin in mag_bins:
            wg = self.mag_bin_numbers == mag_bin
            bin_count = np.ma.count(dg1[wg]) # Counts only unmasked values
            model_resid["mag_bin_number"].append(mag_bin)
            model_resid["bin_count"].append(bin_count)
            model_resid["mean_dg1"].append(np.mean(dg1[wg]))
            model_resid["std_dg1"].append(np.std(dg1[wg])/np.sqrt(bin_count))
            model_resid["mean_dg2"].append(np.mean(dg2[wg]))
            model_resid["std_dg2"].append(np.std(dg2[wg])/np.sqrt(bin_count))
            model_resid["mean_dT"].append(np.mean(dT[wg]))
            model_resid["std_dT"].append(np.std(dT[wg])/np.sqrt(bin_count))
            model_resid["mean_dfwhm"].append(np.mean(dfwhm[wg]))
            model_resid["std_dfwhm"].append(np.std(dfwhm[wg])/np.sqrt(bin_count))

        # Append this model's residuals to the resid_dict
        #full_resid_dict[model] = model_resid

        """ Not sure whether it makes more sense to store full_resid_dict as a
        lazy attribute or just return it to plotter... going to make it an
        attribute so that run statement is a little tidier """
        self.resid_dict = model_resid


    def make_plot(self):
        """ Plot all the things """

        # Set RC params
        fontsize = 15; set_rc_params(fontsize)

        # Grab resid dict
        #model = 'psfex'; model_resid = self.full_resid_dict[model]
        model_resid = self.resid_dict

        # Calculate bin midpoints, list concatenation for the win
        mag_bins = [
            (self.mag_bin_hist[1][i+1] + self.mag_bin_hist[1][i])/2.0
            for i in range(self.nbins)
        ]

        # OK, make plots!
        fig, axs = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

        # Size resid
        axs[0].errorbar(
            mag_bins, model_resid["mean_dT"], yerr=model_resid["std_dT"],
            marker='o', lw=1.5, capsize=5, color='green',
            label=r'$\delta T$',
        )

        axs[0].legend(fontsize=fontsize)
        axs[0].set_ylabel('HSM $T_\mathrm{star}-T_\mathrm{psf}$ (arcsec$^2$)')

        # g1 resid
        axs[1].errorbar(
            mag_bins, model_resid["mean_dg1"], yerr=model_resid["std_dg1"],
            marker='o', lw=1.5, capsize=5, color='red',
            label=r'$\delta g_1$',
        )

        # g2 resid
        axs[1].errorbar(
            mag_bins, model_resid["mean_dg2"], yerr=model_resid["std_dg2"],
            marker='o', lw=1.5, capsize=5, color='blue',
            label=r'$\delta g_2$',
        )

        axs[1].set_ylabel('HSM $g_\mathrm{star}-g_\mathrm{psf}$')
        axs[1].set_xlabel('Magnitude')
        axs[1].legend(fontsize=fontsize)

        plt.tight_layout()

        # Save to file
        fig.savefig(
            os.path.join(
                self.run_config["outdir"], f"{self.model}_mag_psf_resids.pdf"
            )
        )

    def run(self):

        # Create magnitude bins and digitize HSM catalog
        t1_start = time.perf_counter()
        self.bin_by_magnitude()
        t1_stop = time.perf_counter()
        print(f"bin_by_magnitude elapsed time: {t1_stop - t1_start}")

        # Compute residuals and statistics thereof in each mag bin
        t2_start = time.perf_counter()
        self.calculate_statistics()
        t2_stop = time.perf_counter()
        print(f"calculate_statistics elapsed time: {t2_stop - t2_start}")

        # Plot and exit
        t3_start = time.perf_counter()
        self.make_plot()
        t3_stop = time.perf_counter()
        print(f"make_plot elapsed time: {t3_stop - t3_start}")
