import glob
from astropy.io import fits
from astropy.table import Table, vstack
from catalogaugmenter import catalog, psf
from catalogaugmenter import webb_psf, epsfex, shopt, piff_psf 
#from catalogplotter import ResidPlots
import os
import re 
from datetime import datetime, timedelta
import catplot as ctp 
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#Make augmented catalogs with columns for each psf fitter than use the plotter with these new catalogs

#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/test_two_single_exposures/jw0*cal.fits')

#f115w_cat_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/working/mosaic_nircam_f115w_COSMOS-Web_i2d_valid_starcat.fits'
f115w_cat_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/mock_data_outputs_and_cats/mosaic_nircam_f115w_COSMOS-Web_i2d_valid_starcat.fits'
f115w_catalog = catalog(f115w_cat_name)
#f150w_catalog = catalog('/home/eddieberman/research/mcclearygroup/cweb_psf/working/mosaic_nircam_f150w_COSMOS-Web_i2d_valid_starcat.fits')
#f277w_catalog = catalog('/home/eddieberman/research/mcclearygroup/cweb_psf/working/mosaic_nircam_f277w_COSMOS-Web_i2d_valid_starcat.fits')
#f444w_catalog = catalog('/home/eddieberman/research/mcclearygroup/cweb_psf/working/mosaic_nircam_f444w_COSMOS-Web_i2d_valid_starcat.fits')

f115w_epsfex_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/mock_data_outputs_and_cats/mosaic_nircam_f115w_COSMOS-Web_i2d_train_starcat.psf'
epsfex_f115w_object = epsfex(f115w_epsfex_name)
#epsfex_f150w_object = epsfex('/home/eddieberman/research/mcclearygroup/cweb_psf/working/psfex-output/mosaic_nircam_f150w_COSMOS-Web_i2d/mosaic_nircam_f150w_COSMOS-Web_i2d_train_starcat.psf')
#epsfex_f277w_object = epsfex('/home/eddieberman/research/mcclearygroup/cweb_psf/working/psfex-output/mosaic_nircam_f277w_COSMOS-Web_i2d/mosaic_nircam_f277w_COSMOS-Web_i2d_train_starcat.psf')
#epsfex_f444w_object = epsfex('/home/eddieberman/research/mcclearygroup/cweb_psf/working/psfex-output/mosaic_nircam_f444w_COSMOS-Web_i2d/mosaic_nircam_f444w_COSMOS-Web_i2d_train_starcat.psf')

f115w_shopt_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/mock_data_outputs_and_cats/summary.shopt'
shopt_object_f115w = shopt(f115w_shopt_name) 
#shopt_object_f150w = shopt('/home/eddieberman/research/mcclearygroup/shopt/outdir/f150w_mosaic_mock/summary.shopt')
#shopt_object_f277w = shopt('/home/eddieberman/research/mcclearygroup/shopt/outdir/f277w_mosaic_mock/summary.shopt')
#shopt_object_f444w = shopt('/home/eddieberman/research/mcclearygroup/shopt/outdir/f444w_mosaic_mock/summary.shopt')

f115w_piff_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/mock_data_outputs_and_cats/mosaic_nircam_f115w_COSMOS-Web_i2d.piff'
piff_object_f115w = piff_psf(f115w_piff_name)

f115w_catalog.augment(epsfex_f115w_object)
f115w_catalog.augment(shopt_object_f115w)
f115w_catalog.augment(piff_object_f115w, size=75)
f115w_catalog.save_new(outname='f115w_catalog_mosaics_mock_data.fits')
f115w_cat_name = 'f115w_catalog_mosaics_mock_data.fits'


def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

chi_square_visit_psfex = []
chi_square_visit_shopt = []
chi_square_visit_piff = []
 #columns = ['Detector', 'Filter', 'PSFex: SSIM ', 'PSFex: Chi-Square', 'PSFex: Absolute Error', 'PSFex: Relative Error','WebbPSF: SSIM ', 'WebbPSF: Chi-Square', 'WebbPSF: Absolute Error', 'WebbPSF: Relative Error']
columns = ['Detector', 'Filter', 'PSFex: Absolute Error', 'WebbPSF: Absolute Error', 'PSFex: Error', 'WebbPSF: Error','PSFex: Chi-Square', 'WebbPSF: Chi-Square','PSFex: Relative Error', 'WebbPSF: Relative Error', 'PSFex: SSIM ', 'WebbPSF: SSIM ', 'PSFex: FWHM Residual', 'WebbPSF: FWHM Residual', 'PSFex: Reduced Chi Square', 'WebbPSF: Reduced Chi Square']
data_table = Table(names=columns, dtype=['S10', 'S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'])
chi2_name = "f115w_mosaic_mock_chi2_psfex.png"
mre_name = "f115w_mosaic_mock_mre_psfex.png"
abs_name = "f115w_mosaic_mock_abs_psfex.png"

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), epsfex(f115w_epsfex_name))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
mean_absolute_error_plot_psfex.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_psf())
mean_absolute_error_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Normalized Absolute Error = {sum_residuals_abs_psfex}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)
psfex_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_psfex.return_fwhm_star()) - np.array(mean_absolute_error_plot_psfex.return_fwhm_psf()))

abs_name = "f115w_mosaic_mock_abs_shopt.png"
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), shopt(f115w_shopt_name))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_absolute_error_plot_shopt.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_shopt.return_fwhm_psf())
mean_absolute_error_plot_shopt.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Normalized Absolute Error = {sum_residuals_abs_shopt}'])
mean_absolute_error_plot_shopt.save_figure(outname=abs_name)
shopt_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_shopt.return_fwhm_star()) - np.array(mean_absolute_error_plot_shopt.return_fwhm_psf()))

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(catalog(f115w_cat_name), epsfex(f115w_epsfex_name))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
mean_relative_error_plot_psfex.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_relative_error_plot_psfex.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_relative_error_plot_psfex.return_fwhm_psf())
mean_relative_error_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error = {sum_residuals_mre_psfex}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)
psfex_fwhm_residual = np.nanmean(np.array(mean_relative_error_plot_psfex.return_fwhm_star()) - np.array(mean_relative_error_plot_psfex.return_fwhm_psf()))

mre_name = "f115w_mosaic_mock_mre_shopt.png"
mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f115w_cat_name), shopt(f115w_shopt_name))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
mean_relative_error_plot_shopt.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_relative_error_plot_shopt.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_relative_error_plot_shopt.return_fwhm_psf())
mean_relative_error_plot_shopt.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error = {sum_residuals_mre_shopt}'])
mean_relative_error_plot_shopt.save_figure(outname=mre_name)
shopt_fwhm_residual = np.nanmean(np.array(mean_relative_error_plot_shopt.return_fwhm_star()) - np.array(mean_relative_error_plot_shopt.return_fwhm_psf()))

mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f115w_cat_name), epsfex(f115w_epsfex_name))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
mean_chi2_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_psfex}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

chi2_name = "f115w_mosaic_mock_chi2_shopt.png"
mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f115w_cat_name), shopt(f115w_shopt_name))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
mean_chi2_plot_shopt.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_shopt}'])
mean_chi2_plot_shopt.save_figure(outname=chi2_name)
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

mean_chi2_plot_piff = ctp.chi_2_error_plot(catalog(f115w_cat_name), piff_psf(f115w_piff_name))
mean_chi2_plot_piff.preprocessing(hsm_fit=True)
mean_chi2_plot_piff.set_residuals()
sum_residuals_chi2_piff = mean_chi2_plot_piff.return_residuals_sum()
mean_chi2_plot_piff.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_piff}'])
mean_chi2_plot_piff.save_figure(outname=chi2_name)
chi_square_visit_piff += [element for element in mean_chi2_plot_piff.return_chi2_vals()]


newfig, newaxs = plt.subplots(1, 1, figsize=(10, 10))
bins = np.logspace(np.log10(0.1), np.log10(25), 50) # Generates 50 bins between 0.1 and 25 on a log scale.
newaxs.hist(chi_square_visit_psfex, bins=bins, label='PSFex', alpha=0.16)
newaxs.hist(chi_square_visit_shopt, bins=bins, label='shopt', alpha=0.16)
newaxs.hist(chi_square_visit_piff, bins=bins, label='PIFF', alpha=0.16)
newaxs.set_xscale('log')
# Set the limits for the x-axis
newaxs.set_xlim(0.5, 25)
newaxs.set_xlabel('Reduced Chi-Square')
newaxs.set_ylabel('Number of Occurrences')
greater_than_25_psfex = len([value for value in chi_square_visit_psfex if value > 25])
greater_than_25_shopt = len([value for value in chi_square_visit_shopt if value > 25])
greater_than_25_piff = len([value for value in chi_square_visit_piff if value > 25])
newaxs.set_title(f'Mock Data Mosaic F115W Reduced Chi-Square Distribution\n 75 x 75 Pixels')
newaxs.legend()
newfig.savefig('mosaic_reduced_chi_square_distribution.png')
print(f'Not shown are {greater_than_25_psfex} PSFex Reduced chi square values greater than 25\nand {greater_than_25_shopt} Shopt Reduced chi square values greater than 25\n and {greater_than_25_piff} PIFF Reduced chi square values greater than 25')

def sem_with_nans(data):
    """
    Calculate the standard error of the mean (SEM) for a dataset that might contain NaNs.
    
    Args:
    - data (list or np.array): The dataset
    
    Returns:
    - float: The standard error of the mean
    """
    # Remove NaNs
    filtered_data = np.array(data)[~np.isnan(data)]
    
    # Calculate standard deviation
    sd = np.std(filtered_data, ddof=1)  # Using ddof=1 for sample standard deviation
    
    # Calculate SEM
    sem = sd / np.sqrt(len(filtered_data))
    
    return sem

values = [np.nanmean(chi_square_visit_psfex), np.nanmean(chi_square_visit_shopt), np.nanmean(chi_square_visit_piff)]
errors = [sem_with_nans(chi_square_visit_psfex), sem_with_nans(chi_square_visit_shopt), sem_with_nans(chi_square_visit_piff)]

dataset_names = ['PSFex', 'shopt', 'PIFF']


def sem_with_nans(data):
    """
    Calculate the standard error of the mean (SEM) for a dataset that might contain NaNs.
    
    Args:
    - data (list or np.array): The dataset
    
    Returns:
    - float: The standard error of the mean
    """
    filtered_data = np.array(data)[~np.isnan(data)]
    sd = np.std(filtered_data, ddof=1)
    sem = sd / np.sqrt(len(filtered_data))
    return sem

# Data
datasets = ['Comparison Point']
data_samples = [
    chi_square_visit_psfex,
    chi_square_visit_shopt,
    chi_square_visit_piff
]

errors = [sem_with_nans(data) for data in data_samples]

# Colors for different datasets
colors = ['blue', 'red', 'green']

# Plotting
fig, ax = plt.subplots()

for i, (value, error, color) in enumerate(zip(values, errors, colors)):
    ax.errorbar(datasets, value, yerr=error, fmt='o', color=color, ecolor='black', capsize=5, label=f'Dataset {i + 1}')

ax.set_ylabel('Value')
ax.set_title('Comparison of Three Datasets at One Datapoint')
ax.yaxis.grid(True)
ax.legend()

# Save and show the plot
plt.tight_layout()
plt.savefig('overlay_error_bars.png')





'''
426 
427         mean_chi2_plot_psfex = chi_2_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
428         mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
429         mean_chi2_plot_psfex.set_residuals()
430         sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
431         mean_chi2_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Sum Chi2 Error = {sum_residuals_chi2_psfex}'])
432         mean_chi2_plot_psfex.save_figure(outname=chi2_name)
433         reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
434         chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]
435         print(mean_chi2_plot_psfex.return_chi2_vals())
436 
437         mean_mre_plot_psfex = mean_relative_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
438         mean_mre_plot_psfex.preprocessing()
439         mean_mre_plot_psfex.set_residuals()
440         sum_residuals_mre_psfex = mean_mre_plot_psfex.return_residuals_sum()
441         mean_mre_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average MRE Error \n Sum MRE Error = {sum_residuals_mre_psfex}'])
442         mean_mre_plot_psfex.save_figure(outname=mre_name)
443 
'''
'''
f150w_catalog_augmented = f150w_catalog.augment(epsfex_f150w_object)
f150w_catalog_augmented = f150w_catalog_augmented.augment(shopt_object_f150w)
f150w_catalog_augmented.save_new(outname='f150w_catalog_mosaics_mock_data.fits')

f277w_catalog_augmented = f277w_catalog.augment(epsfex_f277w_object)
f277w_catalog_augmented = f277w_catalog_augmented.augment(shopt_object_f277w)
f277w_catalog_augmented.save_new(outname='f277w_catalog_mosaics_mock_data.fits')

f444w_catalog_augmented = f444w_catalog.augment(epsfex_f444w_object)
f444w_catalog_augmented = f444w_catalog_augmented.augment(shopt_object_f444w)
f444w_catalog_augmented.save_new(outname='f444w_catalog_mosaics_mock_data.fits')
'''
