import webbpsf
import psfex
import piff
from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
import numpy as np
from src.box_cutter import BoxCutter
import fitsio
import yaml
from src.utils import read_yaml
import galsim
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error, normalized_root_mse

'''
Some Utility Functions used for Shopt
'''
def objective_function(p, x, y, degree):
    num_coefficients = (degree + 1) * (degree + 2) // 2
    value = 0
    counter = 0
    for a in range(1, degree + 2):
        for b in range(1, degree + 2):
            if (a - 1) + (b - 1) <= degree:
                value += p[counter] * x**(a - 1) * y**(b - 1)
                counter += 1
    return value

def read_shopt(shoptFile):
    f = fits.open(shoptFile)
    polyMatrix = f[0].data
    try: 
        degree = f[1].data['polynomial degree'][0]
    except:
        degree = f[1].data['POLYNOMIAL_DEGREE'][0]
    return polyMatrix, degree

def p(u,v, polMatrix, degree):
    degree = int(degree)
    psf = np.zeros((polMatrix.shape[0], polMatrix.shape[1]))
    for i in range(polMatrix.shape[0]):
        for j in range(polMatrix.shape[1]):
            psf[i,j] = objective_function(polMatrix[i,j,:], u, v, degree)
    return psf/np.sum(psf)

'''
Define a class for catalog. This seems intuitive as catalogs are databases structures 
that lend themselves nicely to object oriented programming with getters and setters to change the database.
'''
class catalog:
    def __init__(self, catalog):
        self.catalog = catalog
        self.sky_level = 0.0 
        self.sky_std = 0.0
        self.data = Table(fits.open(self.catalog)[2].data)
        self.pixel_scale = 0.03
        self.hsm_sig = []
        self.hsm_g1 = []
        self.hsm_g2 = []
        self.fwhm = []

    def augment(self, psf, outname='new_file.fits'):
        current_catalog = fits.open(self.catalog)
        data = Table(current_catalog[2].data)
        new_column_data = []
        for i in range(len(current_catalog[2].data['XWIN_IMAGE'])):
            ext_name1, ext_name2 = psf.coordinate_columns()
            new_column_data.append(psf.render(data[ext_name1][i], data[ext_name2][i]))
        
        new_column = Column(new_column_data, name=psf.nameColumn())
        data.add_column(new_column)
        current_catalog.close()
        self.data = data

        return
    
    def add_noise_flux(self, psf_list, outname='new_file.fits'):
        catalog = fits.open(self.catalog)
        data = Table(catalog[2].data)
        sky_level = self.sky_level
        sky_std = self.sky_std
        for psf in psf_list:
            for i in range(len(data['FLUX_AUTO'])):
                noise = np.random.normal(loc = sky_level, scale = sky_std, size = data[psf.nameColumn()][i].shape)
                data[psf.nameColumn()][i] *= data['FLUX_AUTO'][i]
                data[psf.nameColumn()][i] = data[psf.nameColumn()][i]/np.nansum(data[psf.nameColumn()][i])
                data[psf.nameColumn()][i] += noise
        catalog.close()
        self.data = data 

        return

    def add_err_cutout(self, config, boxcut, image_file, cat_file, ext='ERR', outname='new_file.fits'):
        config = read_yaml(config)
        cat_hdu = config['input_catalog']['hdu']
        box_size = config['box_size']

        # Read in the fits file so that we can add the column ourselves
        catalog = fits.open(self.catalog)
        table = catalog[2].data
        imcat = Table(table) 

        # We need box_size to be same as star size for chi2 calc
        if (box_size != imcat['VIGNET'][0].shape[0]):
            print(f'supplied box_size={box_size} and vignet size={imcat["VIGNET"][0].shape[0]} differ!!')
            print(f'overriding supplied box_size to {imcat["VIGNET"][0].shape[0]}')
            box_size = imcat['VIGNET'][0].shape[0]
            boxcut.box_size = box_size

        # Call to grab_boxes method
        boxes = boxcut.grab_boxes(image_file=image_file, cat_file=imcat)

        data = np.zeros(len(boxes), dtype=[(f'{ext}_VIGNET', 'f4', (box_size, box_size))])
        for i, box in enumerate(boxes):
            data[f'{ext}_VIGNET'][i] = box

        try:
            # Check if the column already exists, if so, replace it
            if f'{ext}_VIGNET' in imcat.colnames:
                imcat.replace_column(f'{ext}_VIGNET', data[f'{ext}_VIGNET'])
            else:
                imcat.add_column(Column(data=data[f'{ext}_VIGNET'], name=f'{ext}_VIGNET'))
        except Exception as e:
            print(f"Error occurred while adding/replacing the column: {e}")
        finally:
            catalog.close()
        
        self.data = imcat
        return

    def crop(self, psf_list, vignet_size=None, replace_original_psf=False):
        '''
        Crop all the VIGNETS in psf_extname to 75x75
        '''
        vignets_colnames = []
        vignets_colnames.append('VIGNET')
        vignets_colnames.append('ERR_VIGNET')
        
        for psf in psf_list:
            vignets_colnames.append(psf.nameColumn())

        def get_middle_pixels(array_2d, n):
            if len(array_2d) != len(array_2d[0]):
                raise ValueError("The input array must be square (n x n)")

            if n <= 0 or n > len(array_2d):
                raise ValueError("Invalid value of n")
            
            start = (len(array_2d) - n) // 2
            end = start + n
            middle_pixels = [row[start:end] for row in array_2d[start:end]]
            return middle_pixels

        vignets_sizes = []
        catalog = fits.open(self.catalog)
        for colname in vignets_colnames:
            vignets_sizes.append(catalog[2].data[colname].shape[1])

        min_dim = min(vignets_sizes)
        print("Minimum dimension of vignets is", min_dim)
        data = Table(catalog[2].data)
        catalog.close()

        if vignet_size is not None:
            min_dim = vignet_size
        
        print("Minimum dimension of vignets is", min_dim)

        for colname in vignets_colnames:
            new_column_data = []
            for i in range(len(data[colname])):
                if data[colname][i].shape[1] > min_dim:
                    new_column_data.append(get_middle_pixels(data[colname][i], min_dim))
                    if i == range(len(data[colname]))[-1]:
                        if replace_original_psf == True:
                            new_column = Column(new_column_data, name=colname)
                            try:
                                data.add_column(new_column)
                            except:
                                data.remove_column(colname)
                                data.add_column(new_column)
                        else:
                            new_column = Column(new_column_data, name=colname+'_CROPPED')
                            try:
                                data.add_column(new_column)
                            except:
                                data.remove_column(colname+'_CROPPED')
                                data.add_column(new_column)
                            print("Cropped column", colname)
        self.data = data
        return
    
    def save_new(self, outname='new_file.fits'):        
        try:
            with fits.open(self.catalog) as original_hdul:
                new_hdul = fits.HDUList(original_hdul)
                new_table_hdu = fits.BinTableHDU(self.data, name='LDAC_OBJECTS')
                new_hdul[2] = new_table_hdu
                new_hdul.writeto(outname, overwrite=True)
        except (IOError, TypeError) as e:
            print("Error occurred:", e)
        return 

    def concatenate_catalogs(self, catalog_new, outname='new_file.fits'):
        catalog1 = fits.open(self.catalog)
        catalog2 = fits.open(catalog_new.catalog)
        data1 = Table(catalog1[2].data)
        data2 = Table(catalog2[2].data)
        data = vstack([data1, data2])
        
        if set(data1.colnames) != set(data2.colnames):
            raise ValueError("Columns in the two catalogs do not match.")

        catalog1.close()
        catalog2.close()
        self.data = data

        return

'''
Define a super class for PSF with the filename. This is useful because we can have
shared functions for things like adding noise but each subclass will have it's own render
method and coordinate system; (x,y) or (u,v) or (ra,dec) etc. 
'''
class psf:
    def __init__(self, psfFileName):
        self.psfFileName = psfFileName
        self.fwhm = []
        self.pixel_scale = 0.03
        self.hsm_sig = []
        self.hsm_g1 = []
        self.hsm_g2 = []


    def render(self):
        pass


    def nameColumn(self):
        pass
    '''
    More Super Methods to come
    '''

class epsfex(psf):
    def __init__(self, psfFileName):
        super().__init__(psfFileName)
        self.psfex = psfex.PSFEx(self.psfFileName)

    def render(self, x, y, vignet_size = 75):
        psfex_model = self.psfex
        return psfex_model.get_rec(y, x)

    def coordinate_columns(self):
        return 'XWIN_IMAGE', 'YWIN_IMAGE'

    def nameColumn(self):
        return 'VIGNET_PSFEX'

class webb_psf(psf):
    def __init__(self, psfFileName):
        super().__init__(psfFileName)

    def render(self, x, y):
        single_psf = fits.open(self.psfFileName)
        extension_names = [ext.name for ext in single_psf]
        if 'DET_DIST' in extension_names:
            ext = 'DET_DIST'
        elif 'PSF_DATA' in extension_names:
            ext = 'PSF_DATA'
        else:
            raise "Extension not found in WebbPSF/SinglePSF model"
        try:
            # Marko's PSFEx format
            single_im = single_psf[ext].data['PSF_MASK'][0, 0, :,:]
        except:
            single_im = single_psf[ext].data
        print(x,y)
        return single_im

    def nameColumn(self):
        return 'VIGNET_WEBBPSF'
    
    def coordinate_columns(self):
        return 'XWIN_IMAGE', 'YWIN_IMAGE'

class shopt(psf):
    def __init__(self, psfFileName):
        super().__init__(psfFileName)
        self.polMatrix, self.degree = read_shopt(self.psfFileName)

    def render(self, u, v):
        return p(u, v, self.polMatrix, self.degree)

    def coordinate_columns(self):
        return 'ALPHAWIN_J2000', 'DELTAWIN_J2000'

    def nameColumn(self):
        return 'VIGNET_SHOPT'

class piff_psf(psf):
    def __init__(self, psfFileName):
        super().__init__(psfFileName)
        self.piff_psf = piff.read(self.psfFileName)

    def render(self, u, v):
        piff_psf = self.piff_psf
        piff_im = piff_psf.draw(u, v)
        return piff_im.array

    def coordinate_columns(self):
        return 'ALPHAWIN_J2000', 'DELTAWIN_J2000'
    
    def nameColumn(self):
        return 'VIGNET_PIFF'

'''
Test Code by Iteratively Adding Columns for rendering PSF stamps of different fitters.

Note: You should always augment before calling the cropping or noise flux methods, otherwise you will be adding noise flux to stamps that don't yet exist. 
I would also recommend specifying the optional outcame arguments so that you don't overwrite your original catalog file, incase you make a mistake 
(or find a mistake in my code).
'''

'''
In general:

catalog_object = catalog('current_file.fits')

psfex_obect = epsfex('/path/model.psf')
piff_object = piff_psf('/path/model.piff')
shopt_object = shopt('/path/summary.shopt')
webb_object = webb_psf('/path/model.fits')

# Set outname to current name to overwrite your existing file instead of creating a new catalog
catalog_object.augment(psfex_object, outname='current_file.fits')  
catalog_object.augment(piff_object, outname='current_file.fits')
catalog_object.augment(shopt_object, outname='current_file.fits') 
catalog_object.augment(webb_object, outname='current_file.fits') 
'''

'''
Here is some of my testing: 

catalog_object = catalog('new_file.fits')
psfex_object = epsfex('working/psfex-output/jw01727116001_04101_00001_nrca3_cal/jw01727116001_04101_00001_nrca3_cal_starcat.psf')
piff_object = piff_psf('/home/eddieberman/research/mcclearygroup/mock_data/mosaics/COSMOS2020_sims/piff-output/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci.piff')
shopt_object = shopt('/home/eddieberman/research/mcclearygroup/shopt/outdir/2023-07-20T11:42:39.026/summary.shopt')
webb_object = webb_psf('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw01727116001_02101_00004_nrcb4_cal_WebbPSF.fits')

#catalog_object.augment(piff_object) 
#catalog_object.augment(psfex_object)
#catalog_object.augment(shopt_object)
catalog_object.augment(webb_object)

catalog_object.crop([psfex_object, piff_object, shopt_object, webb_object])

catalog_object.add_noise_flux([psfex_object, piff_object, shopt_object, webb_object], outname='new_newest_file.fits')

catalog_object.concatenate_catalogs(catalog_object, outname='newest_file.fits') #concatenates with itself

boxcut = BoxCutter(config_file='/home/eddieberman/research/mcclearygroup/cweb_psf/configs/box_cutter.yaml')
catalog_object.add_err_cutout(config='/home/eddieberman/research/mcclearygroup/cweb_psf/configs/box_cutter.yaml', boxcut=boxcut, image_file='single_exposures/jw01727116001_04101_00001_nrca3_cal.fits', cat_file='working/psfex-output/jw01727116001_04101_00001_nrca3_cal/jw01727116001_04101_00001_nrca3_cal_starcat.fits', outname='Added_err.fits')
'''
#print(catalog('new_file.fits').data['VIGNET'])
