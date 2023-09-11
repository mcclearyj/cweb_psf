import unittest
test = unittest.TestCase
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.pyplot as plt
import ipdb, pdb
from astropy.io import fits
from scipy.stats import chi2
import galsim
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error, normalized_root_mse
from catalogutils import AttrDict, set_rc_params
from catalog_hsm_fitter import do_hsm_fit
import catalogaugmenter as ca
import glob
import re
from astropy.table import Table, vstack, hstack, Column
import pandas as pd
import catplot as ctp

directory_path = "/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits"
file_list = glob.glob(directory_path + "/*.fits")

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

def extract_3_numbers(filename):
    pattern = r'\d{3}'
    matches = re.findall(pattern, filename)
    return matches

chi_square_visit_psfex = []
chi_square_visit_webbpsf = []
chi_square_filter_psfex_25_plus = []
chi_square_filter_webbpsf_25_plus = []

#columns = ['Detector', 'Filter', 'PSFex: SSIM ', 'PSFex: Chi-Square', 'PSFex: Absolute Error', 'PSFex: Relative Error','WebbPSF: SSIM ', 'WebbPSF: Chi-Square', 'WebbPSF: Absolute Error', 'WebbPSF: Relative Error']
columns = ['Detector', 'Filter', 'PSFex: Absolute Error', 'WebbPSF: Absolute Error', 'PSFex: Error', 'WebbPSF: Error','PSFex: Chi-Square', 'WebbPSF: Chi-Square','PSFex: Relative Error', 'WebbPSF: Relative Error', 'PSFex: FWHM Residual', 'WebbPSF: FWHM Residual', 'PSFex: Reduced Chi Square (avg)', 'WebbPSF: Reduced Chi Square(avg)', 'PSFex: Reduced Chi Square (median)', 'WebbPSF: Reduced Chi Square (median)']

data_table = Table(names=columns, dtype=['S10', 'S10', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'])

columns = ['Filter_Detector', 'WebbPSF_Mean', 'WebbPSF_SEM', 'PSFex_Mean', 'PSFex_SEM']
df = pd.DataFrame(columns=columns)
'''
mean_abs_plot_psfex = ctp.mean_absolute_error_plot(ca.catalog('validation_thirty_visits/f277_a5_thirty_visits_valid_psf_starcat.fits'), ca.epsfex('epsfex_object_name'))
mean_abs_plot_psfex.preprocessing()
mean_abs_plot_psfex.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_abs_plot_psfex.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_abs_plot_psfex.return_fwhm_psf())
mean_abs_plot_psfex.set_residuals()
sum_residuals_abs_psfex = mean_abs_plot_psfex.return_residuals_sum()
mean_abs_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n MAE = {round(sum_residuals_abs_psfex,3)}'])
mean_abs_plot_psfex.save_figure(outname='test_mae.png')

mean_chi2_plot_psfex = ctp.chi_2_error_plot(ca.catalog('validation_thirty_visits/f277_a5_thirty_visits_valid_psf_starcat.fits'), ca.epsfex('epsfex_object_name'))
mean_chi2_plot_psfex.preprocessing(hsm_fit=True) #False for now
mean_chi2_plot_psfex.set_residuals()
sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
mean_chi2_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi Square \n Chi Square = {round(sum_residuals_chi2_psfex,3)}'])
mean_chi2_plot_psfex.save_figure(outname='test_chi2.png')

mean_mre_plot_psfex = ctp.mean_relative_error_plot(ca.catalog('validation_thirty_visits/f277_a5_thirty_visits_valid_psf_starcat.fits'), ca.epsfex('epsfex_object_name'))
mean_mre_plot_psfex.preprocessing()
mean_mre_plot_psfex.set_residuals()
sum_residuals_mre_psfex = mean_mre_plot_psfex.return_residuals_sum()
mean_mre_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error \n MRE = {round(sum_residuals_mre_psfex,3)}'])
mean_mre_plot_psfex.save_figure(outname='test_mre.png')

mean_abs_plot_webbpsf = ctp.mean_absolute_error_plot(ca.catalog('validation_thirty_visits/f277_a5_thirty_visits_valid_psf_starcat.fits'), ca.webb_psf('WebbPSF'))
mean_abs_plot_webbpsf.preprocessing()
mean_abs_plot_webbpsf.calc_fwhm()
avg_star_fwhm = np.nanmean(mean_abs_plot_webbpsf.return_fwhm_star())
avg_psf_fwhm = np.nanmean(mean_abs_plot_webbpsf.return_fwhm_psf())
mean_abs_plot_webbpsf.set_residuals()
sum_residuals_abs_webbpsf = mean_abs_plot_webbpsf.return_residuals_sum()
mean_abs_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n MAE = {round(sum_residuals_abs_webbpsf,3)}'])
mean_abs_plot_webbpsf.save_figure(outname='test_webbpsf_abs.png')

mean_chi2_plot_webbpsf = ctp.chi_2_error_plot(ca.catalog('validation_thirty_visits/f277_a5_thirty_visits_valid_psf_starcat.fits'), ca.webb_psf('WebbPSF'))
mean_chi2_plot_webbpsf.preprocessing(hsm_fit=True) #False for now
mean_chi2_plot_webbpsf.set_residuals()
sum_residuals_chi2_webbpsf = mean_chi2_plot_webbpsf.return_residuals_sum()
reduced_chi2_webbpsf = np.nanmean(mean_chi2_plot_webbpsf.return_chi2_vals())
median_reduced_chi2_webbpsf = np.nanmedian(mean_chi2_plot_webbpsf.return_chi2_vals())
chi_square_visit_webbpsf += [element for element in mean_chi2_plot_webbpsf.return_chi2_vals()]
mean_chi2_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Mean Reduced Chi2 = {round(reduced_chi2_webbpsf,3)}'])
mean_chi2_plot_webbpsf.save_figure(outname='test_webbpsf_chi2.png')

mean_mre_plot_webbpsf = ctp.mean_relative_error_plot(ca.catalog('validation_thirty_visits/f444_a5_thirty_visits_valid_psf_starcat.fits'), ca.webb_psf('WebbPSF'))
mean_mre_plot_webbpsf.preprocessing()
mean_mre_plot_webbpsf.set_residuals()
sum_residuals_mre_webbpsf = mean_mre_plot_webbpsf.return_residuals_sum()
mean_mre_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error \n MRE = {round(sum_residuals_mre_webbpsf,3)}'])
mean_mre_plot_webbpsf.save_figure(outname='test_webbpsf_mre.png')
print("done")
'''
for file in file_list:
    try:
        if (extract_3_numbers(file)[0] != 444) and (extract_3_numbers(file)[0] != 277):
            
            catalog_object_name = file
            epsfex_object_name = 'epsfex' 
            #It can be any .psf file since we only care about the method that names the column correctly
            
            chi2_name = file.replace(".fits", "psfex_chi2.png")
            mre_name = file.replace(".fits", "psfex_mre_resid.png")
            abs_name = file.replace(".fits", "psfex_abs_resid.png")
            err_name = file.replace(".fits", "psfex_err.png")
            
            mean_absolute_error_plot_psfex = ctp.mean_absolute_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
            mean_absolute_error_plot_psfex.preprocessing()
            mean_absolute_error_plot_psfex.set_residuals()
            sum_residuals_abs_psfex = mean_absolute_error_plot_psfex.return_residuals_sum()
            mean_absolute_error_plot_psfex.calc_fwhm()
            avg_star_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_star())
            avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_psfex.return_fwhm_psf())
            print(f'\n\t Mean PSFex MAE: {np.round(np.nanmean(sum_residuals_abs_psfex),decimals=2)} +- {np.round(mean_absolute_error_plot_psfex.return_sem(),decimals=3)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #mean_absolute_error_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n MAE = {round(sum_residuals_abs_psfex,3)}'])
            #mean_absolute_error_plot_psfex.save_figure(outname=abs_name)
            #psfex_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_psfex.return_fwhm_star()) - np.array(mean_absolute_error_plot_psfex.return_fwhm_psf()))
            '''
           
            mean_chi2_plot_psfex = ctp.chi_2_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
            mean_chi2_plot_psfex.preprocessing(hsm_fit=True) #False for now
            mean_chi2_plot_psfex.set_residuals()
            sum_residuals_chi2_psfex = mean_chi2_plot_psfex.return_residuals_sum()
            reduced_chi2_psfex = np.nanmean(mean_chi2_plot_psfex.return_chi2_vals())
            median_reduced_chi2_psfex = np.nanmedian(mean_chi2_plot_psfex.return_chi2_vals())
            chi_square_visit_psfex += [element for element in mean_chi2_plot_psfex.return_chi2_vals()]
            for i, val in enumerate(mean_chi2_plot_psfex.return_chi2_vals()):
                if val > 25:
                    chi_square_filter_psfex_25_plus.append(re.findall(r'[A-Za-z]\d', file)[1])
                    chi_square_filter_psfex_25_plus.append('F'+extract_3_numbers(file)[0]+'W_')
                    chi_square_filter_psfex_25_plus.append(i)
            #mean_chi2_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Mean Reduced Chi2 = {round(reduced_chi2_psfex,3)}'])
            #mean_chi2_plot_psfex.save_figure(outname='chi2_name_cmap_test.png')
            
            #print(mean_chi2_plot_psfex.return_chi2_vals())
            #print(f'\n\t Median PSFex Chi2 Error: {round(median_reduced_chi2_psfex,2)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #print(f'\n\t Mean PSFex Chi2 Error: {round(reduced_chi2_psfex,2)} +- {round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            
            mean_mre_plot_psfex = ctp.mean_relative_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
            mean_mre_plot_psfex.preprocessing()
            mean_mre_plot_psfex.set_residuals()
            sum_residuals_mre_psfex = mean_mre_plot_psfex.return_residuals_sum()
            print(f'\n\t Mean PSFex MRE: {np.round(sum_residuals_mre_psfex,decimals=2)} +- {np.round(mean_mre_plot_psfex.return_sem(),decimals=2)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #mean_mre_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error \n MRE = {round(sum_residuals_mre_psfex,3)}'])
            #mean_mre_plot_psfex.save_figure(outname=mre_name)

            mean_error_plot_psfex = ctp.mean_absolute_error_plot(ca.catalog(catalog_object_name), ca.epsfex(epsfex_object_name))
            mean_error_plot_psfex.preprocessing()
            mean_error_plot_psfex.set_residuals()
            sum_residuals_err_psfex = mean_error_plot_psfex.return_residuals_sum()
            #mean_error_plot_psfex.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Error \n Sum Absolute Error = {sum_residuals_abs_psfex}'])
            #mean_error_plot_psfex.save_figure(outname=err_name)
            '''
            webbpsf_object_name = 'WebbPSF'
            
            chi2_name = file.replace(".fits", "webb_chi2.png")
            mre_name = file.replace(".fits", "webb_mre_resid.png")
            abs_name = file.replace(".fits", "webb_abs_resid.png")
            err_name = file.replace(".fits", "webb_err.png")
            
            mean_absolute_error_plot_webbpsf = ctp.mean_absolute_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
            mean_absolute_error_plot_webbpsf.preprocessing()
            mean_absolute_error_plot_webbpsf.set_residuals()
            sum_residuals_abs_webb = mean_absolute_error_plot_webbpsf.return_residuals_sum()
            mean_absolute_error_plot_webbpsf.calc_fwhm()
            avg_star_fwhm = np.nanmean(mean_absolute_error_plot_webbpsf.return_fwhm_star())
            avg_psf_fwhm = np.nanmean(mean_absolute_error_plot_webbpsf.return_fwhm_psf())
            print(f'\n\t Mean WebbPSF Absolute Error: {np.round(sum_residuals_abs_webb,decimals=2)} +- {np.round(mean_absolute_error_plot_webbpsf.return_sem(),decimals=3)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #mean_absolute_error_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Absolute Error \n MAE = {round(sum_residuals_abs_webb,3)}'])
            #mean_absolute_error_plot_webbpsf.save_figure(outname=abs_name)
            #webbpsf_fwhm_residual = np.nanmean(np.array(mean_absolute_error_plot_webbpsf.return_fwhm_star()) - np.array(mean_absolute_error_plot_webbpsf.return_fwhm_psf()))
            '''
            
            mean_chi2_plot_webbpsf = ctp.chi_2_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
            mean_chi2_plot_webbpsf.preprocessing(hsm_fit=True) #hsm_fit=False for now
            mean_chi2_plot_webbpsf.set_residuals()
            sum_residuals_chi2_webb = mean_chi2_plot_webbpsf.return_residuals_sum()
            reduced_chi2_webb = np.nanmean(mean_chi2_plot_webbpsf.return_chi2_vals())
            median_reduced_chi2_webb = np.nanmedian(mean_chi2_plot_webbpsf.return_chi2_vals())
            chi_square_visit_webbpsf += [element for element in mean_chi2_plot_webbpsf.return_chi2_vals()]
            for i, val in enumerate(mean_chi2_plot_webbpsf.return_chi2_vals()):
                if val > 25:
                    chi_square_filter_webbpsf_25_plus.append(re.findall(r'[A-Za-z]\d', file)[1])
                    chi_square_filter_webbpsf_25_plus.append('F'+extract_3_numbers(file)[0]+'W_')
                    chi_square_filter_webbpsf_25_plus.append(i)
            #mean_chi2_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Chi2 Error \n Mean Reduced Chi2 Error = {round(reduced_chi2_webb,3)}'])
            #mean_chi2_plot_webbpsf.save_figure(outname=chi2_name)
            
            #print(f'\tMedian WebbPSF Chi2 Error: {round(median_reduced_chi2_webb,3)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #print(f'\tMean WebbPSF Chi2 Error: {round(reduced_chi2_webb,2)}  +- {round(sem_with_nans(mean_chi2_plot_webbpsf.return_chi2_vals()),2)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            
            mean_mre_plot_webbpsf = ctp.mean_relative_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
            mean_mre_plot_webbpsf.preprocessing()
            mean_mre_plot_webbpsf.set_residuals()
            sum_residuals_mre_webb = mean_mre_plot_webbpsf.return_residuals_sum()
            print(f'\tMean WebbPSF MRE: {np.round(sum_residuals_mre_webb,decimals=2)}  +- {np.round(mean_mre_plot_webbpsf.return_sem(),decimals=2)}', file.replace('/home/eddieberman/research/mcclearygroup/cweb_psf/validation_thirty_visits/', ''))
            #mean_mre_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Relative Error \n MRE = {round(sum_residuals_mre_webb,3)}'])
            #mean_mre_plot_webbpsf.save_figure(outname=mre_name)
            '''

            mean_error_plot_webbpsf = ctp.mean_error_plot(ca.catalog(catalog_object_name), ca.webb_psf(webbpsf_object_name))
            mean_error_plot_webbpsf.preprocessing()
            mean_error_plot_webbpsf.set_residuals()
            sum_residuals_err_webb = mean_error_plot_webbpsf.return_residuals_sum()
            #mean_error_plot_webbpsf.set_titles([f'Average Star\n Avg FWHM = {avg_star_fwhm}', f'Average PSF\n Avg FWHM = {avg_psf_fwhm}', f'Average Error \n Sum Error = {sum_residuals_err_webb}'])
            #mean_error_plot_webbpsf.save_figure(outname=err_name)
            #new_row = {'Filter_Detector': f'F{extract_3_numbers(file)[0]}W_'+re.findall(r"[A-Za-z]\d", file)[1], 'WebbPSF_Mean': round(reduced_chi2_webb,2), 'WebbPSF_SEM': round(sem_with_nans(mean_chi2_plot_webbpsf.return_chi2_vals()),2), 'PSFex_Mean': round(reduced_chi2_psfex,2), 'PSFex_SEM': round(sem_with_nans(mean_chi2_plot_psfex.return_chi2_vals()),2)}
            #df = df.append(new_row, ignore_index=True)
            #print(df)
            #letter_number_codes = re.findall(r'[A-Za-z]\d', file)
            #print('\n', 'F'+extract_3_numbers(file)[0]+'W_', letter_number_codes[1], ' Sum Absolute Error PSFex = ', sum_residuals_abs_psfex, ' Sum Absolute Error WebbPSF = ', sum_residuals_abs_webb, 'Sum Error PSFex =', sum_residuals_err_psfex, 'Sum Error WebbPSF', sum_residuals_err_webb, ' Sum Chi2 Error PSFex = ', f'{sum_residuals_chi2_psfex:.1f}', ' Sum Chi2 Error WebbPSF = ', f'{sum_residuals_chi2_webb:.1f}', ' Sum MRE Error PSFex = ', sum_residuals_mre_psfex, ' Sum MRE Error WebbPSF = ', sum_residuals_mre_webb, 'FWHM Residual PSFex', psfex_fwhm_residual, 'FWHM Residual WebbPSF', webbpsf_fwhm_residual, 'Reduced Chi Square PSFex', f'{reduced_chi2_psfex:.1f}', 'Reduced Chi Square WebbPSF', f'{reduced_chi2_webb:.1f}')
            #data_table.add_row(['F'+extract_3_numbers(file)[0]+'W', letter_number_codes[1], round(sum_residuals_abs_psfex, 3), round(sum_residuals_abs_webb, 3), round(sum_residuals_err_psfex, 3), round(sum_residuals_err_webb, 3), round(sum_residuals_chi2_psfex, 3), round(sum_residuals_chi2_webb, 3), round(sum_residuals_mre_psfex, 3), round(sum_residuals_mre_webb, 3), psfex_fwhm_residual, webbpsf_fwhm_residual, round(reduced_chi2_psfex, 3), round(reduced_chi2_webb, 3), round(median_reduced_chi2_psfex, 3), round(median_reduced_chi2_webb, 3)])
    except Exception as e:
        print("An error occurred: ", e)
'''
newfig, newaxs = plt.subplots(1, 1, figsize=(10, 10))
bins = np.logspace(np.log10(0.1), np.log10(100), 100) # Generates 50 bins between 0.1 and 25 on a log scale.
newaxs.hist(chi_square_visit_psfex, bins=bins, label='PSFex', alpha=0.25)
newaxs.hist(chi_square_visit_webbpsf, bins=bins, label='WebbPSF', alpha=0.25)
newaxs.set_xscale('log')
# Set the limits for the x-axis
newaxs.set_xlim(1, 25)
newaxs.set_xlabel('Reduced Chi-Square')
newaxs.set_ylabel('Number of Occurrences')
greater_than_25_psfex = len([value for value in chi_square_visit_psfex if value > 25])
greater_than_100_psfex = len([value for value in chi_square_visit_psfex if value > 100])
greater_than_25_webbpsf = len([value for value in chi_square_visit_webbpsf if value > 25])
greater_than_100_webbpsf = len([value for value in chi_square_visit_webbpsf if value > 100])
#greater_than_25 = len(chi_square_visit_psfex[chi_square_visit_psfex > 25])
newaxs.set_title(f'30 Aggregated Visits Reduced Chi-Square Distribution\nNot shown are {greater_than_100_psfex} PSFex Reduced chi square values greater than 100\nand {greater_than_100_webbpsf} WebbPSF Reduced chi square values greater than 100')
newaxs.legend()
newfig.savefig('validation_thirty_visits/reduced_chi_square_distribution.png')
print("done")
df['Filter'] = df['Filter_Detector'].str.split('_').str[0]
df['Detector'] = df['Filter_Detector'].str.split('_').str[1]

# Sort by Filter and then by Detector
df = df.sort_values(by=['Filter', 'Detector']).reset_index(drop=True)
df = df.drop(columns=['Filter', 'Detector'])

fig, ax = plt.subplots(figsize=(15, 7))

positions = range(len(df['Filter_Detector']))
ax.set_yscale('log')

# Plotting WebbPSF data
ax.errorbar(positions, df['WebbPSF_Mean'], yerr=df['WebbPSF_SEM'], fmt='o', label='WebbPSF', capsize=5, capthick=1, color='blue')

# Plotting PSFex data (with a slight x shift for better visibility)
ax.errorbar([p + 0.3 for p in positions], df['PSFex_Mean'], yerr=df['PSFex_SEM'], fmt='o', label='PSFex', capsize=5, capthick=1, color='red')

# Customization
ax.set_title('Comparison of WebbPSF and PSFex')
ax.set_xlabel('Filter Detector Combinations')
ax.set_ylabel('Mean Reduced Chi Square')
ax.set_xticks([p + 0.15 for p in positions])  
ax.set_xticklabels(df['Filter_Detector'], rotation=45, ha="right")  
ax.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
fig.savefig('mean_comparison.png')

copy_data_table = data_table.copy()
#columns = ['Detector', 'Filter', 'PSFex: Absolute Error', 'WebbPSF: Absolute Error', 'PSFex: Error', 'WebbPSF: Error','PSFex: Chi-Square', 'WebbPSF: Chi-Square','PSFex: Relative Error', 'WebbPSF: Relative Error', 'PSFex: SSIM ', 'WebbPSF: SSIM ', 'PSFex: FWHM Residual', 'WebbPSF: FWHM Residual', 'PSFex: Reduced Chi Square', 'WebbPSF: Reduced Chi Square']
#'PSFex: Reduced Chi Square (avg)', 'WebbPSF: Reduced Chi Square(avg)'
copy_data_table.remove_columns(['PSFex: Chi-Square', 'WebbPSF: Chi-Square','PSFex: Error', 'WebbPSF: Error', 'PSFex: FWHM Residual', 'WebbPSF: FWHM Residual'])

cp_df = copy_data_table.to_pandas()
sorted_cp_df = cp_df.sort_values(by=['Filter', 'Detector'])
copy_data_table = Table.from_pandas(sorted_cp_df)

copy_data_table.write('output_statistics.tab', format='latex', overwrite=True)

df = data_table.to_pandas()
sorted_df = df.sort_values(by=['Filter', 'Detector'])
pd.set_option("display.max_columns", None)  # Display all columns
pd.set_option("display.width", None)  # Disable column width truncation

# Print the DataFrame (which will show all columns)
print(sorted_df)
sorted_df.to_csv('output_statistics.csv', index=False)
# Optionally, convert the DataFrame back to an Astropy table
table = Table.from_pandas(sorted_df)
hdu = fits.PrimaryHDU()
table_hdu = fits.table_to_hdu(table)
hdu = fits.HDUList([hdu, table_hdu])
hdu.writeto('output_statistics.fits', overwrite=True)
newfig, newaxs = plt.subplots(1, 1, figsize=(10, 10))
bins = np.logspace(np.log10(0.1), np.log10(100), 100) # Generates 50 bins between 0.1 and 25 on a log scale.
newaxs.hist(chi_square_visit_psfex, bins=bins, label='PSFex', alpha=0.25)
newaxs.hist(chi_square_visit_webbpsf, bins=bins, label='WebbPSF', alpha=0.25)
newaxs.set_xscale('log')
# Set the limits for the x-axis
newaxs.set_xlim(0.1, 25)
newaxs.set_xlabel('Reduced Chi-Square')
newaxs.set_ylabel('Number of Occurrences')
greater_than_25_psfex = len([value for value in chi_square_visit_psfex if value > 25])
greater_than_100_psfex = len([value for value in chi_square_visit_psfex if value > 100])
greater_than_25_webbpsf = len([value for value in chi_square_visit_webbpsf if value > 25])
greater_than_100_webbpsf = len([value for value in chi_square_visit_webbpsf if value > 100])
#greater_than_25 = len(chi_square_visit_psfex[chi_square_visit_psfex > 25])
newaxs.set_title(f'30 Aggregated Visits Reduced Chi-Square Distribution\nNot shown are {greater_than_100_psfex} PSFex Reduced chi square values greater than 100\nand {greater_than_100_webbpsf} WebbPSF Reduced chi square values greater than 100')
newaxs.legend()
newfig.savefig('validation_thirty_visits/reduced_chi_square_distribution.png')
print('PSFex reduced chi square > 25',chi_square_filter_psfex_25_plus)
print('WebbPSF reduced chi square > 25',chi_square_filter_webbpsf_25_plus)
#box_whisker_plot = newaxs.boxplot([chi_square_visit_psfex, chi_square_visit_webbpsf]) 
#newaxs.set_xticks([1, 2])
#newaxs.set_xticklabels(['PSFex Reduced Chi-Square', 'WebbPSF Reduced Chi-Square'])
#newaxs.set_xlabel('Data Sets')
#newaxs.set_ylabel('Values')
#newaxs.set_title('Visit 153 Reduced Chi-Square Distribution')
#newfig.savefig('reduced_chi_square_distribution.png')
'''
