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

f115w_catalog_name = 'new_f115_apr_mosaic_combined_catalog.fits'
f115w_catalog = catalog(f115w_catalog_name)

f150w_catalog_name = 'new_f150_apr_mosaic_combined_catalog.fits'
f150w_catalog = catalog(f150w_catalog_name)

f277w_catalog_name = 'new_f277_apr_mosaic_combined_catalog.fits'
f277w_catalog = catalog(f277w_catalog_name)

f444w_catalog_name = 'new_f444_apr_mosaic_combined_catalog.fits'
f444w_catalog = catalog(f444w_catalog_name)

def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

chi_square_visit_psfex = []
chi_square_visit_shopt = []
chi_square_visit_piff = []
columns = ['Filter', 'PSFex_Avg_Reduced_Chisq', 'PSFex_STD_Reduced_Chisq', 'PSFex_Median_Reduced_Chisq', 'PSFex_MAE', 'PSFex_MAE_STD', 'PSFex_MRE' , 'PSFex_MRE_STD','ShOpt_Avg_Reduced_Chisq', 'ShOpt_STD_Reduced_Chisq', 'ShOpt_Median_Reduced_Chisq', 'ShOpt_MAE', 'ShOpt_MAE_STD', 'ShOpt_MRE' , 'ShOpt_MRE_STD']
df = pd.DataFrame(columns=columns)

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(f115w_catalog, epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F115W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname="f115w_mosaic_real_mre_psfex.png")

mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(f115w_catalog, shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname="f115w_mosaic_real_mre_shopt.png")

mean_chi2_plot_psfex = ctp.chi_2_error_plot(f115w_catalog, epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F115W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname="f115w_mosaic_real_chi2_psfex.png")
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

mean_chi2_plot_shopt = ctp.chi_2_error_plot(f115w_catalog, shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname="f115w_mosaic_real_chi2_shopt.png")
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(f115w_catalog, epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
mean_absolute_error_plot_psfex.set_titles(['Average Star F115W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname="f115w_mosaic_real_abs_psfex.png")

mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(f115w_catalog, shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.set_titles(['Average Star F115W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname="f115w_mosaic_real_abs_shopt.png")

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

mean_relative_error_plot_psfs = ctp.mean_relative_error_plot(f150w_catalog, epsfex(''))
mean_relative_error_plot_psfs.preprocessing()
mean_relative_error_plot_psfs.set_residuals()
sum_residuals_mre_psfs = mean_relative_error_plot_psfs.return_residuals_sum()
std_mr_psfs = mean_relative_error_plot_psfs.return_sem()
mean_relative_error_plot_psfs.set_titles([f'Average Star F150W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfs,2)} +/- {round(std_mr_psfs,2)}'])
mean_relative_error_plot_psfs.save_figure(outname="f150w_mosaic_real_mre_psfex.png")

mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(f150w_catalog, shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F150W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname="f150w_mosaic_real_mre_shopt.png")

mean_chi2_plot_psfs = ctp.chi_2_error_plot(f150w_catalog, epsfex(''))
mean_chi2_plot_psfs.preprocessing(hsm_fit=True)
mean_chi2_plot_psfs.set_residuals()
sum_residuals_chi2_psfs = mean_chi2_plot_psfs.return_residuals_sum()
reduced_chi2_psfs = np.nanmean(mean_chi2_plot_psfs.return_chi2_vals())
median_reduced_chi2_psfs = np.nanmedian(mean_chi2_plot_psfs.return_chi2_vals())
std_chi_psfs = sem_with_nans(mean_chi2_plot_psfs.return_chi2_vals())
mean_chi2_plot_psfs.set_titles([f'Average Star F150W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfs, 2)} +/- {round(std_chi_psfs,2)}'])
mean_chi2_plot_psfs.save_figure(outname="f150w_mosaic_real_chi2_psfex.png")
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfs.return_chi2_vals()]

mean_chi2_plot_shopt = ctp.chi_2_error_plot(f150w_catalog, shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F150W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname="f150w_mosaic_real_chi2_shopt.png")
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

mean_absolute_error_plot_psfs = ctp.mean_absolute_error_plot(f150w_catalog, epsfex(''))
mean_absolute_error_plot_psfs.preprocessing()
mean_absolute_error_plot_psfs.set_residuals()
sum_residuals_abs_psfs = mean_absolute_error_plot_psfs.return_residuals_sum()
std_abs_psfs = mean_absolute_error_plot_psfs.return_residuals_sum()
mean_absolute_error_plot_psfs.set_titles(['Average Star F150W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfs,2)} +/- {round(std_abs_psfs,2)}'])
mean_absolute_error_plot_psfs.save_figure(outname="f150w_mosaic_real_abs_psfex.png")

mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(f150w_catalog, shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.set_titles(['Average Star F150W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname="f150w_mosaic_real_abs_shopt.png")

new_row = {
    'Filter': 'F150W',
    'PSFex_Avg_Reduced_Chisq': round(reduced_chi2_psfs, 2),
    'PSFex_STD_Reduced_Chisq': round(std_chi_psfs, 2),
    'PSFex_Median_Reduced_Chisq': round(median_reduced_chi2_psfs, 2),
    'PSFex_MAE': round(sum_residuals_abs_psfs, 2),
    'PSFex_MAE_STD': round(std_abs_psfs, 2),
    'PSFex_MRE': round(sum_residuals_mre_psfs, 2),
    'PSFex_MRE_STD': round(std_mr_psfs, 2),
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

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(f277w_catalog, epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F277W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname="f277w_mosaic_real_mre_psfex.png")

mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(f277w_catalog, shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F277W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname="f277w_mosaic_real_mre_shopt.png")

mean_chi2_plot_psfex = ctp.chi_2_error_plot(f277w_catalog, epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F277W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname="f277w_mosaic_real_chi2_psfex.png")
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

mean_chi2_plot_shopt = ctp.chi_2_error_plot(f277w_catalog, shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F277W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname="f277w_mosaic_real_chi2_shopt.png")
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(f277w_catalog, epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
mean_absolute_error_plot_psfex.set_titles(['Average Star F277W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname="f277w_mosaic_real_abs_psfex.png")

mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(f277w_catalog, shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.set_titles(['Average Star F277W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname="f277w_mosaic_real_abs_shopt.png")

new_row = {
    'Filter': 'F277W',
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

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(f444w_catalog, epsfex(''))
mean_relative_error_plot_psfex.preprocessing()
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Star F444W', f'Average PSFex PSF', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname="f444w_mosaic_real_mre_psfex.png")

mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(f444w_catalog, shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F444W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname="f444w_mosaic_real_mre_shopt.png")

mean_chi2_plot_psfex = ctp.chi_2_error_plot(f444w_catalog, epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Star F444W', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname="f444w_mosaic_real_chi2_psfex.png")
chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]

mean_chi2_plot_shopt = ctp.chi_2_error_plot(f444w_catalog, shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F444W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname="f444w_mosaic_real_chi2_shopt.png")
chi_square_visit_shopt += [element for element in mean_chi2_plot_shopt.return_chi2_vals()]

mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(f444w_catalog, epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
mean_absolute_error_plot_psfex.set_titles(['Average Star F444W', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname="f444w_mosaic_real_abs_psfex.png")

mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(f444w_catalog, shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.set_titles(['Average Star F444W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname="f444w_mosaic_real_abs_shopt.png")

new_row = {
    'Filter': 'F444W',
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

df.to_csv('real_mosaic_results.csv', index=False)

newfig, newaxs = plt.subplots(1, 1, figsize=(10, 10))
bins = np.logspace(np.log10(0.1), np.log10(100), 200) # Generates 50 bins between 0.1 and 25 on a log scale.
newaxs.hist(chi_square_visit_psfex, bins=bins, label='PSFex', alpha=0.33)
newaxs.hist(chi_square_visit_shopt, bins=bins, label='shopt', alpha=0.33)
newaxs.set_xscale('log')
# Set the limits for the x-axis
newaxs.set_xlim(0.5, 100)
newaxs.set_xlabel(r'$\chi^2$', fontsize=30)
newaxs.set_ylabel('Number of Occurrences', fontsize=30)
greater_than_100_psfex = len([value for value in chi_square_visit_psfex if value > 100])
greater_than_100_shopt = len([value for value in chi_square_visit_shopt if value > 100])
newaxs.set_title(f'April Mosaic ' + r'$\chi^2$' + ' Distribution', fontsize=35)
newaxs.legend()
newfig.savefig('apr_real_mosaic_reduced_chi_square_distribution.png')
print(f'Not shown are {greater_than_100_psfex} PSFex Reduced chi square values greater than 100\nand {greater_than_100_shopt} Shopt Reduced chi square values greater than 100')

fig, ax = plt.subplots(figsize=(15, 7))
positions = range(len(df['Filter']))
#ax.set_yscale('log')

# Plotting PIFF data
#ax.errorbar(positions, df['PIFF_Avg_Reduced_Chisq'], yerr=df['PIFF_STD_Reduced_Chisq'], 
 #           fmt='o', label='PIFF', capsize=5, capthick=1, color='blue')

# Plotting PSFex data (with a slight x shift for better visibility)
ax.errorbar([p + 0.2 for p in positions], df['PSFex_Avg_Reduced_Chisq'], yerr=df['PSFex_STD_Reduced_Chisq'], 
            fmt='o', label='PSFex', capsize=5, capthick=1, color='red')

# Plotting ShOpt data (with another slight x shift for better visibility)
ax.errorbar([p + 0.4 for p in positions], df['ShOpt_Avg_Reduced_Chisq'], yerr=df['ShOpt_STD_Reduced_Chisq'], 
            fmt='o', label='ShOpt', capsize=5, capthick=1, color='green')

# Customization      
ax.set_title('Comparison of PSFex, and ShOpt ' + r'$\overline{\chi^2}$', fontsize=25)
ax.set_xlabel('Filter', fontsize=30)
ax.set_ylabel('Average ' + r'$\chi^2$', fontsize=30)
ax.set_xticks([p + 0.2 for p in positions])  
ax.set_xticklabels(df['Filter'], rotation=45, ha="right")  
ax.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
fig.savefig('apr_real_mosaic_mean_comparison.png')

