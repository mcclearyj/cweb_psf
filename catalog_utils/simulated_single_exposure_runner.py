import glob
from astropy.io import fits
from astropy.table import Table, vstack
from catalogaugmenter import catalog, psf
from catalogaugmenter import webb_psf, epsfex, shopt, piff_psf 
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
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def sem_with_nans(data):        
    filtered_data = np.array(data)[~np.isnan(data)]        
    sd = np.std(filtered_data, ddof=1)          
    sem = sd / np.sqrt(len(filtered_data))        
    return sem       
#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/test_two_single_exposures/jw0*cal.fits')

def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

#mirage_catalog = catalog('f115w_sse_combined_catalog.fits')
mirage_catalog = catalog('new_f444_sse_combined_catalog.fits')

#mre_name = "f115w_sse_psfex.png"
mre_name = "f444w_sse_psfex.png"

'''
mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), epsfex(''))
mean_absolute_error_plot_psfex.preprocessing()
mean_absolute_error_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
std_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
mean_absolute_error_plot_psfex.set_titles(['Average Mirage Cutout', 'Average PSFex PSF', f'MAE = {round(sum_residuals_abs_psfex,2)} +/- {round(std_abs_psfex,2)}'])
mean_absolute_error_plot_psfex.save_figure(outname=abs_name)
'''

mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, epsfex(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average PSFex PSF F444w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mre_name = "f444w_sse_shopt.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, shopt(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average ShOpt PSF F444w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mirage_catalog = catalog('new_f115_sse_combined_catalog.fits')
mre_name = "f115w_sse_psfex.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, epsfex(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average PSFex PSF F115w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mre_name = "f115w_sse_shopt.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, shopt(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average ShOpt PSF F115w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mirage_catalog = catalog('new_f150_sse_combined_catalog.fits')
mre_name = "f150w_sse_psfex.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, epsfex(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average PSFex PSF F150w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mre_name = "f150w_sse_shopt.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, shopt(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average ShOpt PSF F150w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mirage_catalog = catalog('new_f277_sse_combined_catalog.fits')
mre_name = "f277w_sse_psfex.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, epsfex(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average PSFex PSF F277w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)

mre_name = "f277w_sse_shopt.png"
mean_relative_error_plot_psfex = ctp.mean_relative_error_plot(mirage_catalog, shopt(''), Mirage=True)
mean_relative_error_plot_psfex.preprocessing(hsm_fit=True)
mean_relative_error_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_relative_error_plot_psfex.return_residuals_sum()
std_mr_psfex = mean_relative_error_plot_psfex.return_sem()
mean_relative_error_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average ShOpt PSF F277w', f'MRE = {round(sum_residuals_mre_psfex,2)} +/- {round(std_mr_psfex,2)}'])
mean_relative_error_plot_psfex.save_figure(outname=mre_name)


'''
mean_chi2_plot_psfex = ctp.chi_2_error_plot(catalog(f115w_cat_name), epsfex(''))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True)
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
std_chi_psfex = sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals())
mean_chi2_plot_psfex.set_titles([f'Average Mirage Cutout', f'Average PSFex PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_psfex, 2)} +/- {round(std_chi_psfex,2)}'])
mean_chi2_plot_psfex.save_figure(outname=chi2_name)
'''

'''
mean_absolute_error_plot_shopt = ctp.mean_absolute_error_plot(catalog(f115w_cat_name), shopt(''))
mean_absolute_error_plot_shopt.preprocessing()
mean_absolute_error_plot_shopt.set_residuals()
sum_residuals_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
std_abs_shopt = mean_absolute_error_plot_shopt.return_residuals_sum()
mean_absolute_error_plot_shopt.set_titles(['Average Star F115W', 'Average ShOpt PSF', f'MAE = {round(sum_residuals_abs_shopt,2)} +/- {round(std_abs_shopt,2)}'])
mean_absolute_error_plot_shopt.save_figure(outname='f115w_mosaic_real_abs_shopt.png')

mean_relative_error_plot_shopt = ctp.mean_relative_error_plot(catalog(f115w_cat_name), shopt(''))
mean_relative_error_plot_shopt.preprocessing()
mean_relative_error_plot_shopt.set_residuals()
sum_residuals_mre_shopt = mean_relative_error_plot_shopt.return_residuals_sum()
std_mr_shopt = mean_relative_error_plot_shopt.return_sem()
mean_relative_error_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', f'MRE = {round(sum_residuals_mre_shopt,2)} +/- {round(std_mr_shopt,2)}'])
mean_relative_error_plot_shopt.save_figure(outname="f115w_mosaic_real_mre_shopt.png")

mean_chi2_plot_shopt = ctp.chi_2_error_plot(catalog(f115w_cat_name), shopt(''))
mean_chi2_plot_shopt.preprocessing(hsm_fit=True)
mean_chi2_plot_shopt.set_residuals()
sum_residuals_chi2_shopt = mean_chi2_plot_shopt.return_residuals_sum()
reduced_chi2_shopt = np.nanmean(mean_chi2_plot_shopt.return_chi2_vals())
median_reduced_chi2_shopt = np.nanmedian(mean_chi2_plot_shopt.return_chi2_vals())
std_chi_shopt = sem_with_nans(mean_chi2_plot_shopt.return_chi2_vals())
mean_chi2_plot_shopt.set_titles([f'Average Star F115W', f'Average ShOpt PSF', r'$\overline{\chi^2}$' + f' = {round(reduced_chi2_shopt, 2)} +/- {round(std_chi_shopt,2)}'])
mean_chi2_plot_shopt.save_figure(outname="f115w_mosaic_real_chi2_shopt.png")
'''

