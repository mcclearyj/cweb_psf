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
import pandas as pd
import pdb
#Make augmented catalogs with columns for each psf fitter than use the plotter with these new catalogs
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None) # or a large number like 1000

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
#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/test_two_single_exposures/jw0*cal.fits')

#f115w_cat_name = '/home/eddieberman/research/mcclearygroup/cweb_psf/working/mosaic_nircam_f115w_COSMOS-Web_i2d_valid_starcat.fits'
f115w_cat_name = 'new_f115_sse_combined_catalog.fits'
f115w_catalog = catalog(f115w_cat_name)

f150w_cat_name = 'new_f150_sse_combined_catalog.fits'
f150w_catalog = catalog(f150w_cat_name)

f277w_cat_name = 'new_f277_sse_combined_catalog.fits'
f277w_catalog = catalog(f277w_cat_name)

f444w_cat_name = 'new_f444_sse_combined_catalog.fits'
f444w_catalog = catalog(f444w_cat_name)


def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

chi_square_visit_psfex = []
chi_square_visit_shopt = []
#chi_square_visit_piff = []
#columns = ['Filter', 'PIFF_Avg_Reduced_Chisq', 'PIFF_STD_Reduced_Chisq', 'PIFF_Median_Reduced_Chisq', 'PIFF_MAE', 'PIFF_MAE_STD', 'PIFF_MRE' , 'PIFF_MRE_STD','PSFex_Avg_Reduced_Chisq', 'PSFex_STD_Reduced_Chisq', 'PSFex_Median_Reduced_Chisq', 'PSFex_MAE', 'PSFex_MAE_STD', 'PSFex_MRE' , 'PSFex_MRE_STD','ShOpt_Avg_Reduced_Chisq', 'ShOpt_STD_Reduced_Chisq', 'ShOpt_Median_Reduced_Chisq', 'ShOpt_MAE', 'ShOpt_MAE_STD', 'ShOpt_MRE' , 'ShOpt_MRE_STD']

# Columns without PIFF
columns = ['Filter', 'PSFex_Avg_Reduced_Chisq', 'PSFex_STD_Reduced_Chisq', 'PSFex_Median_Reduced_Chisq', 'PSFex_MAE', 'PSFex_MAE_STD', 'PSFex_MRE' , 'PSFex_MRE_STD','ShOpt_Avg_Reduced_Chisq', 'ShOpt_STD_Reduced_Chisq', 'ShOpt_Median_Reduced_Chisq', 'ShOpt_MAE', 'ShOpt_MAE_STD', 'ShOpt_MRE' , 'ShOpt_MRE_STD']

df = pd.DataFrame(columns=columns)


chi2_name = "f115w_sse_chi2_psfex.pdf"
mre_name = "f115w_sse_mre_psfex.pdf"
abs_name = "f115w_sse_abs_psfex.pdf"

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_sem()
mean_absolute_error_plot_psfex.set_titles(['Average Star F115W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)

'''
abs_name='f115w_sse_abs_piff.png'
mean_absolute_error_plot_piff = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), piff_psf(''))
mean_absolute_error_plot_piff.preprocessing()
mean_absolute_error_plot_piff.set_residuals()
sum_residuals_abs_piff = mean_absolute_error_plot_piff.return_residuals_sum()
std_abs_piff = mean_absolute_error_plot_piff.return_sem()
mean_absolute_error_plot_piff.set_titles(['Average Star F115W', 'Average PIFF PSF', f'MAE = {sum_residuals_abs_piff} +/- {round(std_abs_piff,2)}'])
mean_absolute_error_plot_piff.save_figure(outname=abs_name)
print(f'Average Normalized Absolute Error = {round(sum_residuals_abs_piff,2)} +/- {round(std_abs_piff,2)}')
'''

abs_name = "f115w_sse_abs_shopt.png"
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_sem()
mean_absolute_error_plot_shopt.set_titles([f'Average Star F115W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname=abs_name)


mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(catalog(f115w_cat_name), epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F115W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

'''
mre_name = "f115w_sse_mre_piff.png"
mean_relative_error_plot_piff = ctp.mean_relative_error_plot(catalog(f115w_cat_name), piff_psf(''))
mean_relative_error_plot_piff.preprocessing()
mean_relative_error_plot_piff.set_residuals()
sum_residuals_mre_piff = mean_relative_error_plot_piff.return_residuals_sum()
std_mr_piff = mean_relative_error_plot_piff.return_sem()
mean_relative_error_plot_piff.set_titles([f'Average Star F115W', f'Average PIFF PSF', f'MRE = {round(sum_residuals_mre_piff,2)} +/- {round(std_mr_piff,2)}'])
mean_relative_error_plot_piff.save_figure(outname=mre_name)
'''

mre_name = "f115w_sse_mre_shopt.png"
mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f115w_cat_name), shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname=mre_name)

#+- {round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}
chi2_name =  "f115w_sse_chi2_psfex.png"
mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f115w_cat_name), epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F115W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]


chi2_name = "f115w_sse_chi2_shopt.png"
mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f115w_cat_name), shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname=chi2_name)
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

'''
chi2_name = "f115w_sse_chi2_piff.png"
mean_chi2_plot_piff = ctp.chi_2_error_plot(catalog(f115w_cat_name), piff_psf(''))
mean_chi2_plot_piff.preprocessing(hsm_fit=True)
mean_chi2_plot_piff.set_residuals()
sum_residuals_chi2_piff = mean_chi2_plot_piff.return_residuals_sum()
reduced_chi2_piff = np.nanmean(mean_chi2_plot_piff.return_chi2_vals())
median_reduced_chi2_piff = np.nanmedian(mean_chi2_plot_piff.return_chi2_vals())
std_chi_piff = sem_with_nans(mean_chi2_plot_piff.return_chi2_vals())
mean_chi2_plot_piff.set_titles([f'Average Star F115W', f'Average PIFF PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_piff, 2)} +/- {round(std_chi_piff,2)}'])
mean_chi2_plot_piff.save_figure(outname=chi2_name)
chi_square_visit_piff += [element for element in mean_chi2_plot_piff.return_chi2_vals()]
'''

#{round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}']

print('F115W')
print(f' PSFex Avg Reduced Chi2 = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex, 2)}, PSFex Median Reduced Chi2 = {round(median_reduced_chi2_psfex, 2)}, PSFex MAE = {round(sum_residuals_abs_psfex, 2)} +/- {round(std_abs_psfex, 2)}, PSFex MRE = {round(sum_residuals_mre_psfex, 2)} +/- {round(std_mr_psfex, 2)}, ShOpt Avg Reduced Chi2 = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt, 2)}, ShOpt Median Reduced Chi2 = {round(median_reduced_chi2_shopt, 2)}, ShOpt MAE = {round(sum_residuals_abs_shopt, 2)} +/- {round(std_abs_shopt, 2)}, ShOpt MRE = {round(sum_residuals_mre_shopt, 2)} +/- {round(std_mr_shopt, 2)}')


new_row = {
    'Filter': 'F115W',
    'PSFex_Avg_Reduced_Chisq': round(reduced_chi2_psfex, 2),
    'PSFex_STD_Reduced_Chisq': round(std_chi_psfex, 2),
    'PSFex_Median_Reduced_Chisq': round(median_reduced_chi2_psfex, 2),
    'PSFex_MAE': round(sum_residuals_abs_psfex, 2),
    'PSFex_MAE_STD': round(std_abs_psfex, 2),
    'PSFex_MRE': round(sum_residuals_mre_psfex, 2),
    'PSFex_MRE_STD': round(std_mr_psfex, 2),
    'ShOpt_Avg_Reduced_Chisq': round(reduced_chi2_shopt, 2),
    'ShOpt_STD_Reduced_Chisq': round(std_chi_shopt, 2),
    'ShOpt_Median_Reduced_Chisq': round(median_reduced_chi2_shopt, 2),
    'ShOpt_MAE': round(sum_residuals_abs_shopt, 2),
    'ShOpt_MAE_STD': round(std_abs_shopt, 2),
    'ShOpt_MRE': round(sum_residuals_mre_shopt, 2),
    'ShOpt_MRE_STD': round(std_mr_shopt, 2)
}

df = df.append(new_row, ignore_index=True)
print(df)
print('done F115W')


chi2_name = "f150w_sse_chi2_psfex.png"
mre_name = "f150w_sse_mre_psfex.png"
abs_name = "f150w_sse_abs_psfex.png"

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f150w_cat_name), epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_sem()
mean_absolute_error_plot_psfex.set_titles(['Average Star F150W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)

'''
abs_name='f150w_sse_abs_piff.png'
mean_absolute_error_plot_piff = ctp.mean_absolute_error_plot(catalog(f150w_cat_name), piff_psf(''))
mean_absolute_error_plot_piff.preprocessing()
mean_absolute_error_plot_piff.set_residuals()
sum_residuals_abs_piff = mean_absolute_error_plot_piff.return_residuals_sum()
std_abs_piff = mean_absolute_error_plot_piff.return_sem()
mean_absolute_error_plot_piff.set_titles(['Average Star F150W', 'Average PIFF PSF', f'MAE = {sum_residuals_abs_piff} +/- {round(std_abs_piff,2)}'])
mean_absolute_error_plot_piff.save_figure(outname=abs_name)
print(f'Average Normalized Absolute Error = {round(sum_residuals_abs_piff,2)} +/- {round(std_abs_piff,2)}')
'''

abs_name = "f150w_sse_abs_shopt.png"
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f150w_cat_name), shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_sem()
mean_absolute_error_plot_shopt.set_titles([f'Average Star F150W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname=abs_name)


mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(catalog(f150w_cat_name), epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F150W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

'''
mre_name = "f150w_sse_mre_piff.png"
mean_relative_error_plot_piff = ctp.mean_relative_error_plot(catalog(f150w_cat_name), piff_psf(''))
mean_relative_error_plot_piff.preprocessing()
mean_relative_error_plot_piff.set_residuals()
sum_residuals_mre_piff = mean_relative_error_plot_piff.return_residuals_sum()
std_mr_piff = mean_relative_error_plot_piff.return_sem()
mean_relative_error_plot_piff.set_titles([f'Average Star F150W', f'Average PIFF PSF', f'MRE = {round(sum_residuals_mre_piff,2)} +/- {round(std_mr_piff,2)}'])
mean_relative_error_plot_piff.save_figure(outname=mre_name)
'''

mre_name = "f150w_sse_mre_shopt.png"
mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f150w_cat_name), shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F150W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname=mre_name)

#+- {round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}
chi2_name =  "f150w_sse_chi2_psfex.png"
mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f150w_cat_name), epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = np.nanstd(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F150W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

chi2_name = "f150w_sse_chi2_shopt.png"
mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f150w_cat_name), shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = np.nanstd(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F150W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname=chi2_name)
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

'''
chi2_name = "f150w_sse_chi2_piff.png"
mean_chi2_plot_piff = ctp.chi_2_error_plot(catalog(f150w_cat_name), piff_psf(''))
mean_chi2_plot_piff.preprocessing(hsm_fit=True)
mean_chi2_plot_piff.set_residuals()
sum_residuals_chi2_piff = mean_chi2_plot_piff.return_residuals_sum()
reduced_chi2_piff = np.nanmean(mean_chi2_plot_piff.return_chi2_vals())
median_reduced_chi2_piff = np.nanmedian(mean_chi2_plot_piff.return_chi2_vals())
std_chi_piff = np.nanstd(mean_chi2_plot_piff.return_chi2_vals())
mean_chi2_plot_piff.set_titles([f'Average Star F150W', f'Average PIFF PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_piff, 2)} +/- {round(std_chi_piff,2)}'])
mean_chi2_plot_piff.save_figure(outname=chi2_name)
chi_square_visit_piff += [element for element in mean_chi2_plot_piff.return_chi2_vals()]
'''

print("F150W")
print(f' PSFex Avg Reduced Chi2 = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex, 2)}, PSFex Median Reduced Chi2 = {round(median_reduced_chi2_psfex, 2)}, PSFex MAE = {round(sum_residuals_abs_psfex, 2)} +/- {round(std_abs_psfex, 2)}, PSFex MRE = {round(sum_residuals_mre_psfex, 2)} +/- {round(std_mr_psfex, 2)}, ShOpt Avg Reduced Chi2 = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt, 2)}, ShOpt Median Reduced Chi2 = {round(median_reduced_chi2_shopt, 2)}, ShOpt MAE = {round(sum_residuals_abs_shopt, 2)} +/- {round(std_abs_shopt, 2)}, ShOpt MRE = {round(sum_residuals_mre_shopt, 2)} +/- {round(std_mr_shopt, 2)}')

#{round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}']
new_row = {
    'Filter': 'F150W',
    #'PIFF_Avg_Reduced_Chisq': round(reduced_chi2_piff, 2),
    #'PIFF_STD_Reduced_Chisq': round(std_chi_piff, 2),
    #'PIFF_Median_Reduced_Chisq': round(median_reduced_chi2_piff, 2),
    #'PIFF_MAE': round(sum_residuals_abs_piff, 2),
    #'PIFF_MAE_STD': round(std_abs_piff, 2),
    #'PIFF_MRE': round(sum_residuals_mre_piff, 2),
    #'PIFF_MRE_STD': round(std_mr_piff, 2),
    'PSFex_Avg_Reduced_Chisq': round(reduced_chi2_psfex, 2),
    'PSFex_STD_Reduced_Chisq': round(std_chi_psfex, 2),
    'PSFex_Median_Reduced_Chisq': round(median_reduced_chi2_psfex, 2),
    'PSFex_MAE': round(sum_residuals_abs_psfex, 2),
    'PSFex_MAE_STD': round(std_abs_psfex, 2),
    'PSFex_MRE': round(sum_residuals_mre_psfex, 2),
    'PSFex_MRE_STD': round(std_mr_psfex, 2),
    'ShOpt_Avg_Reduced_Chisq': round(reduced_chi2_shopt, 2),
    'ShOpt_STD_Reduced_Chisq': round(std_chi_shopt, 2),
    'ShOpt_Median_Reduced_Chisq': round(median_reduced_chi2_shopt, 2),
    'ShOpt_MAE': round(sum_residuals_abs_shopt, 2),
    'ShOpt_MAE_STD': round(std_abs_shopt, 2),
    'ShOpt_MRE': round(sum_residuals_mre_shopt, 2),
    'ShOpt_MRE_STD': round(std_mr_shopt, 2)
}

df = df.append(new_row, ignore_index=True)
print(df)
print('done F150W')


chi2_name = "f277w_sse_chi2_psfex.png"
mre_name = "f277w_sse_mre_psfex.png"
abs_name = "f277w_sse_abs_psfex.png"

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f277w_cat_name), epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_sem()
mean_absolute_error_plot_psfex.set_titles(['Average Star F277W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)

'''
abs_name='f277w_sse_abs_piff.png'
mean_absolute_error_plot_piff = ctp.mean_absolute_error_plot(catalog(f277w_cat_name), piff_psf(''))
mean_absolute_error_plot_piff.preprocessing()
mean_absolute_error_plot_piff.set_residuals()
sum_residuals_abs_piff = mean_absolute_error_plot_piff.return_residuals_sum()
std_abs_piff = mean_absolute_error_plot_piff.return_sem()
mean_absolute_error_plot_piff.set_titles(['Average Star F277W', 'Average PIFF PSF', f'MAE = {sum_residuals_abs_piff} +/- {round(std_abs_piff,2)}'])
mean_absolute_error_plot_piff.save_figure(outname=abs_name)
print(f'Average Normalized Absolute Error = {round(sum_residuals_abs_piff,2)} +/- {round(std_abs_piff,2)}')
'''

abs_name = "f277w_sse_abs_shopt.png"
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f277w_cat_name), shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_sem()
mean_absolute_error_plot_shopt.set_titles([f'Average Star F150W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname=abs_name)

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(catalog(f277w_cat_name), epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F277W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

'''
mre_name = "f277w_sse_mre_piff.png"
mean_relative_error_plot_piff = ctp.mean_relative_error_plot(catalog(f277w_cat_name), piff_psf(''))
mean_relative_error_plot_piff.preprocessing()
mean_relative_error_plot_piff.set_residuals()
sum_residuals_mre_piff = mean_relative_error_plot_piff.return_residuals_sum()
std_mr_piff = mean_relative_error_plot_piff.return_sem()
mean_relative_error_plot_piff.set_titles([f'Average Star F277W', f'Average PIFF PSF', f'MRE = {round(sum_residuals_mre_piff,2)} +/- {round(std_mr_piff,2)}'])
mean_relative_error_plot_piff.save_figure(outname=mre_name)
'''

mre_name = "f277w_sse_mre_shopt.png"
mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f277w_cat_name), shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F277W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname=mre_name)

#+- {round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}
chi2_name =  "f277w_sse_chi2_psfex.png"
mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f277w_cat_name), epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = np.nanstd(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F277W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

chi2_name = "f277w_sse_chi2_shopt.png"
mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f277w_cat_name), shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = np.nanstd(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F277W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname=chi2_name)
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

'''
chi2_name = "f277w_sse_chi2_piff.png"
mean_chi2_plot_piff = ctp.chi_2_error_plot(catalog(f277w_cat_name), piff_psf(''))
mean_chi2_plot_piff.preprocessing(hsm_fit=True)
mean_chi2_plot_piff.set_residuals()
sum_residuals_chi2_piff = mean_chi2_plot_piff.return_residuals_sum()
reduced_chi2_piff = np.nanmean(mean_chi2_plot_piff.return_chi2_vals())
median_reduced_chi2_piff = np.nanmedian(mean_chi2_plot_piff.return_chi2_vals())
std_chi_piff = np.nanstd(mean_chi2_plot_piff.return_chi2_vals())
mean_chi2_plot_piff.set_titles([f'Average Star F277W', f'Average PIFF PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_piff, 2)} +/- {round(std_chi_piff,2)}'])
mean_chi2_plot_piff.save_figure(outname=chi2_name)
chi_square_visit_piff += [element for element in mean_chi2_plot_piff.return_chi2_vals()]
'''

print("F277W")
print(f' PSFex Avg Reduced Chi2 = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex, 2)}, PSFex Median Reduced Chi2 = {round(median_reduced_chi2_psfex, 2)}, PSFex MAE = {round(sum_residuals_abs_psfex, 2)} +/- {round(std_abs_psfex, 2)}, PSFex MRE = {round(sum_residuals_mre_psfex, 2)} +/- {round(std_mr_psfex, 2)}, ShOpt Avg Reduced Chi2 = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt, 2)}, ShOpt Median Reduced Chi2 = {round(median_reduced_chi2_shopt, 2)}, ShOpt MAE = {round(sum_residuals_abs_shopt, 2)} +/- {round(std_abs_shopt, 2)}, ShOpt MRE = {round(sum_residuals_mre_shopt, 2)} +/- {round(std_mr_shopt, 2)}')



#{round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}']
new_row = {
    'Filter': 'F277W',
    #'PIFF_Avg_Reduced_Chisq': round(reduced_chi2_piff, 2),
    #'PIFF_STD_Reduced_Chisq': round(std_chi_piff, 2),
    #'PIFF_Median_Reduced_Chisq': round(median_reduced_chi2_piff, 2),
    #'PIFF_MAE': round(sum_residuals_abs_piff, 2),
    #'PIFF_MAE_STD': round(std_abs_piff, 2),
    #'PIFF_MRE': round(sum_residuals_mre_piff, 2),
    #'PIFF_MRE_STD': round(std_mr_piff, 2),
    'PSFex_Avg_Reduced_Chisq': round(reduced_chi2_psfex, 2),
    'PSFex_STD_Reduced_Chisq': round(std_chi_psfex, 2),
    'PSFex_Median_Reduced_Chisq': round(median_reduced_chi2_psfex, 2),
    'PSFex_MAE': round(sum_residuals_abs_psfex, 2),
    'PSFex_MAE_STD': round(std_abs_psfex, 2),
    'PSFex_MRE': round(sum_residuals_mre_psfex, 2),
    'PSFex_MRE_STD': round(std_mr_psfex, 2),
    'ShOpt_Avg_Reduced_Chisq': round(reduced_chi2_shopt, 2),
    'ShOpt_STD_Reduced_Chisq': round(std_chi_shopt, 2),
    'ShOpt_Median_Reduced_Chisq': round(median_reduced_chi2_shopt, 2),
    'ShOpt_MAE': round(sum_residuals_abs_shopt, 2),
    'ShOpt_MAE_STD': round(std_abs_shopt, 2),
    'ShOpt_MRE': round(sum_residuals_mre_shopt, 2),
    'ShOpt_MRE_STD': round(std_mr_shopt, 2)
}

df = df.append(new_row, ignore_index=True)
print(df)
print('done F277W')


chi2_name = "f444w_sse_chi2_psfex.png"
mre_name = "f444w_sse_mre_psfex.png"
abs_name = "f444w_sse_abs_psfex.png"

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f444w_cat_name), epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_sem()
mean_absolute_error_plot_psfex.set_titles(['Average Star F444W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)

'''
abs_name='f444w_sse_abs_piff.png'
mean_absolute_error_plot_piff = ctp.mean_absolute_error_plot(catalog(f444w_cat_name), piff_psf(''))
mean_absolute_error_plot_piff.preprocessing()
mean_absolute_error_plot_piff.set_residuals()
sum_residuals_abs_piff = mean_absolute_error_plot_piff.return_residuals_sum()
std_abs_piff = mean_absolute_error_plot_piff.return_sem()
mean_absolute_error_plot_piff.set_titles(['Average Star F444W', 'Average PIFF PSF', f'MAE = {sum_residuals_abs_piff} +/- {round(std_abs_piff,2)}'])
mean_absolute_error_plot_piff.save_figure(outname=abs_name)
print(f'Average Normalized Absolute Error = {round(sum_residuals_abs_piff,2)} +/- {round(std_abs_piff,2)}')
'''

abs_name = "f444w_sse_abs_shopt.png"
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f444w_cat_name), shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_sem()
mean_absolute_error_plot_shopt.set_titles([f'Average Star F444W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname=abs_name)

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(catalog(f444w_cat_name), epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F444W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

'''
mre_name = "f444w_sse_mre_piff.png"
mean_relative_error_plot_piff = ctp.mean_relative_error_plot(catalog(f444w_cat_name), piff_psf(''))
mean_relative_error_plot_piff.preprocessing()
mean_relative_error_plot_piff.set_residuals()
sum_residuals_mre_piff = mean_relative_error_plot_piff.return_residuals_sum()
std_mr_piff = mean_relative_error_plot_piff.return_sem()
mean_relative_error_plot_piff.set_titles([f'Average Star F444W', f'Average PIFF PSF', f'MRE = {round(sum_residuals_mre_piff,2)} +/- {round(std_mr_piff,2)}'])
mean_relative_error_plot_piff.save_figure(outname=mre_name)
'''

mre_name = "f444w_sse_mre_shopt.png"
mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f444w_cat_name), shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F444W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname=mre_name)

#+- {round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}
chi2_name =  "f444w_sse_chi2_psfex.png"
mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f444w_cat_name), epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = np.nanstd(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F444W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

chi2_name = "f444w_sse_chi2_shopt.png"
mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f444w_cat_name), shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = np.nanstd(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F444W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname=chi2_name)
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

'''
chi2_name = "f444w_sse_chi2_piff.png"
mean_chi2_plot_piff = ctp.chi_2_error_plot(catalog(f444w_cat_name), piff_psf(''))
mean_chi2_plot_piff.preprocessing(hsm_fit=True)
mean_chi2_plot_piff.set_residuals()
sum_residuals_chi2_piff = mean_chi2_plot_piff.return_residuals_sum()
reduced_chi2_piff = np.nanmean(mean_chi2_plot_piff.return_chi2_vals())
median_reduced_chi2_piff = np.nanmedian(mean_chi2_plot_piff.return_chi2_vals())
std_chi_piff = np.nanstd(mean_chi2_plot_piff.return_chi2_vals())
mean_chi2_plot_piff.set_titles([f'Average Star F444W', f'Average PIFF PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_piff, 2)} +/- {round(std_chi_piff,2)}'])
mean_chi2_plot_piff.save_figure(outname=chi2_name)
chi_square_visit_piff += [element for element in mean_chi2_plot_piff.return_chi2_vals()]
'''

print("F444W")
print(f' PSFex Avg Reduced Chi2 = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex, 2)}, PSFex Median Reduced Chi2 = {round(median_reduced_chi2_psfex, 2)}, PSFex MAE = {round(sum_residuals_abs_psfex, 2)} +/- {round(std_abs_psfex, 2)}, PSFex MRE = {round(sum_residuals_mre_psfex, 2)} +/- {round(std_mr_psfex, 2)}, ShOpt Avg Reduced Chi2 = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt, 2)}, ShOpt Median Reduced Chi2 = {round(median_reduced_chi2_shopt, 2)}, ShOpt MAE = {round(sum_residuals_abs_shopt, 2)} +/- {round(std_abs_shopt, 2)}, ShOpt MRE = {round(sum_residuals_mre_shopt, 2)} +/- {round(std_mr_shopt, 2)}')


#{round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}']
new_row = {
    'Filter': 'F444W',
    #'PIFF_Avg_Reduced_Chisq': round(reduced_chi2_piff, 2),
    #'PIFF_STD_Reduced_Chisq': round(std_chi_piff, 2),
    #'PIFF_Median_Reduced_Chisq': round(median_reduced_chi2_piff, 2),
    #'PIFF_MAE': round(sum_residuals_abs_piff, 2),
    #'PIFF_MAE_STD': round(std_abs_piff, 2),
    #'PIFF_MRE': round(sum_residuals_mre_piff, 2),
    #'PIFF_MRE_STD': round(std_mr_piff, 2),
    'PSFex_Avg_Reduced_Chisq': round(reduced_chi2_psfex, 2),
    'PSFex_STD_Reduced_Chisq': round(std_chi_psfex, 2),
    'PSFex_Median_Reduced_Chisq': round(median_reduced_chi2_psfex, 2),
    'PSFex_MAE': round(sum_residuals_abs_psfex, 2),
    'PSFex_MAE_STD': round(std_abs_psfex, 2),
    'PSFex_MRE': round(sum_residuals_mre_psfex, 2),
    'PSFex_MRE_STD': round(std_mr_psfex, 2),
    'ShOpt_Avg_Reduced_Chisq': round(reduced_chi2_shopt, 2),
    'ShOpt_STD_Reduced_Chisq': round(std_chi_shopt, 2),
    'ShOpt_Median_Reduced_Chisq': round(median_reduced_chi2_shopt, 2),
    'ShOpt_MAE': round(sum_residuals_abs_shopt, 2),
    'ShOpt_MAE_STD': round(std_abs_shopt, 2),
    'ShOpt_MRE': round(sum_residuals_mre_shopt, 2),
    'ShOpt_MRE_STD': round(std_mr_shopt, 2)
}

df = df.append(new_row, ignore_index=True)
print(df)
print('done F444W')

df.to_csv('sse_results.csv', index=False)

tab = Table(
      [chi_square_visit_psfex, chi_square_visit_shopt],
      names=['psfex_chi2', 'shopt_chi2']
      )
tab.write('sse_chi2_vals.csv', format='ascii.csv', overwrite=True)

# tab = Table.read('simulated_mosaics_chi2_vals')
fig, axs = plt.subplots(1, 1, figsize=(8, 7), tight_layout=True)
plt.rcParams.update({'xtick.labelsize': 18})
plt.rcParams.update({'ytick.labelsize': 18})
bins = np.logspace(np.log10(1), np.log10(100), 20) # Generates 50 bins between 0.1 and 50 on a log scale.
#bins = np.arange(0.1, 100)
axs.hist(chi_square_visit_psfex, bins=bins, label='PSFex',
            alpha=0.5, range=[1, 100], histtype='bar')
axs.hist(chi_square_visit_shopt, bins=bins, label=r'\texttt{ShOpt}',
            alpha=0.5, range=[1, 100], histtype='bar')
#axs.hist(chi_square_visit_piff, bins=bins, label='PIFF',
#            alpha=0.5, range=[0, 100], histtype='bar')
axs.set_xscale('log')
# Set the limits for the x-axis
axs.set_xlim(1.5, 110)
axs.set_xlabel(r'$\chi^2$', fontsize=18)
axs.set_ylabel('Number of Occurrences', fontsize=18)
greater_than_100_psfex = len([value for value in chi_square_visit_psfex if value > 100])
greater_than_100_shopt = len([value for value in chi_square_visit_shopt if value > 100])
#greater_than_100_piff = len([value for value in chi_square_visit_piff if value > 100])
axs.set_title(f'Simulated Single Exposure ' + r'$\chi^2$' + ' Distribution', fontsize=18)
axs.legend()
fig.savefig('sse_reduced_chi_square_distribution.pdf')

ShOpt = r"\texttt{ShOpt}"
print(f'Not shown are {greater_than_100_psfex} PSFex reduced $\chi^2$ values greater than 100 and {greater_than_100_shopt} {ShOpt} reduced $\chi^2$ values greater than 100')



plt.rcParams.update({'xtick.major.size': 8})
plt.rcParams.update({'xtick.major.width': 1.3})
plt.rcParams.update({'xtick.minor.visible': True})
plt.rcParams.update({'xtick.minor.width': 1.1})
plt.rcParams.update({'xtick.minor.size': 6})
plt.rcParams.update({'xtick.direction': 'out'})
plt.rcParams.update({'ytick.major.width': 1.3})
plt.rcParams.update({'ytick.major.size': 8})
plt.rcParams.update({'ytick.minor.visible': True})
plt.rcParams.update({'ytick.minor.width': 1.1})
plt.rcParams.update({'ytick.minor.size':6})
plt.rcParams.update({'xtick.labelsize': 22})
plt.rcParams.update({'ytick.labelsize': 22})

fig, ax = plt.subplots(figsize=(15, 7), linewidth=2)
positions = range(len(df['Filter']))
ax.set_yscale('log')

# Plotting PIFF data
#ax.errorbar(positions, df['PIFF_Avg_Reduced_Chisq'], yerr=df['PIFF_STD_Reduced_Chisq'],
#            fmt='o', label='PIFF', capsize=5, capthick=1, color='blue')

# Plotting PSFex data (with a slight x shift for better visibility)
ax.errorbar(positions, df['PSFex_Avg_Reduced_Chisq'], yerr=df['PSFex_STD_Reduced_Chisq'],
            fmt='o', label='PSFex', capsize=5, ms=8, capthick=1.1, color='red', lw=2)

# Plotting ShOpt data (with another slight x shift for better visibility)
ax.errorbar([p + 0.2 for p in positions], df['ShOpt_Avg_Reduced_Chisq'], yerr=df['ShOpt_STD_Reduced_Chisq'],
            fmt='o', label='ShOpt', capsize=5, ms=8, capthick=1.1, color='green', lw=2)

# Customization
ax.set_title('Comparison of PSFex and ShOpt ' + r'$\overline{\chi^2}$', fontsize=25)
ax.set_xlabel('Filter', fontsize=25)
ax.set_ylabel('Average ' + r'$\chi^2$', fontsize=25)
ax.set_xticks([p + 0.2 for p in positions])
ax.set_xticklabels(df['Filter'], rotation=45, ha="right")

ax.legend(fontsize=22)
plt.grid(True, which='both', linestyle='--', linewidth=1)
plt.tight_layout()
fig.savefig('sse_mean_comparison.png')
