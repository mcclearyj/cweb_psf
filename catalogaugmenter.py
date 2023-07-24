import webbpsf
import psfex
import piff
from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
import numpy as np

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

    def augment(self, psf):
        '''
        Add outdir argument?
        '''
        current_catalog = fits.open(self.catalog)
        data = Table(current_catalog[2].data)
        new_column_data = []
        for i in range(len(current_catalog[2].data['XWIN_IMAGE'])):
            ext_name1, ext_name2 = psf.coordinate_columns()
            new_column_data.append(psf.render(data[ext_name1][i], data[ext_name2][i]))
            #print(data[ext_name1][i], data[ext_name2][i])
        new_column = Column(new_column_data, name=psf.nameColumn())
        data.add_column(new_column)
        current_catalog.close()
        try:
            with fits.open(self.catalog) as original_hdul:
                new_hdul = fits.HDUList(original_hdul)
                new_table_hdu = fits.BinTableHDU(data, name='LDAC_OBJECTS')
                new_hdul[2] = new_table_hdu
                new_hdul.writeto('new_file.fits', overwrite=True)
        except (IOError, TypeError) as e:
            print("Error occurred:", e)
    
    def add_noise(self, psf_extname):
        '''
        Look at Vignets to get noise, then replace all the VIGNETS in psf_extname with themselves + noise
        '''
        pass

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
        try:
            with fits.open(self.catalog) as original_hdul:
                new_hdul = fits.HDUList(original_hdul)
                new_table_hdu = fits.BinTableHDU(data, name='LDAC_OBJECTS')
                new_hdul[2] = new_table_hdu
                new_hdul.writeto('new_file.fits', overwrite=True)
        except (IOError, TypeError) as e:
            print("Error occurred:", e)
        
        

'''
Define a super class for PSF with the filename. This is useful because we can have
shared functions for things like adding noise but each subclass will have it's own render
method and coordinate system; (x,y) or (u,v) or (ra,dec) etc. 
'''
class psf:
    def __init__(self, psfFileName):
        self.psfFileName = psfFileName

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

class webbpsf(psf):
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
Test Code by iteratively adding a VIGNETS column for each PSF Model
'''

catalog_object = catalog('new_file.fits')
psfex_object = epsfex('working/psfex-output/jw01727116001_04101_00001_nrca3_cal/jw01727116001_04101_00001_nrca3_cal_starcat.psf')
piff_object = piff_psf('/home/eddieberman/research/mcclearygroup/mock_data/mosaics/COSMOS2020_sims/piff-output/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci.piff')
shopt_object = shopt('/home/eddieberman/research/mcclearygroup/shopt/outdir/2023-07-20T11:42:39.026/summary.shopt')
webb_object = webbpsf('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw01727116001_02101_00004_nrcb4_cal_WebbPSF.fits')
#catalog_object.augment(piff_object)
#catalog_object.augment(psfex_object)
#catalog_object.augment(shopt_object)
#catalog_object.augment(webb_object)
'''
You should add noise to the image before you crop it, or when you crop it make sure to crop the error image as well
'''
#catalog_object.crop([psfex_object, piff_object, shopt_object, webb_object])
catalog_object.crop([piff_object])
