import webbpsf
import psfex
import piff
from astropy.io import fits
from astropy.table import Table, vstack, hstack, Column
import numpy as np

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
    degree = f[1].data['polynomial degree'][0]
    return polyMatrix, degree

def p(u,v, polMatrix, degree):
    degree = int(degree)
    psf = np.zeros((polMatrix.shape[0], polMatrix.shape[1]))
    for i in range(polMatrix.shape[0]):
        for j in range(polMatrix.shape[1]):
            psf[i,j] = objective_function(polMatrix[i,j,:], u, v, degree)
    return psf/np.sum(psf)

class catalog:
    def __init__(self, catalog):
        self.catalog = catalog

    def augment(self, psf):
        '''
        '''
        current_catalog = fits.open(self.catalog)
        data = Table(current_catalog[2].data)
        new_column_data = []
        for i in range(len(current_catalog[2].data['XWIN_IMAGE'])):
            ext_name1, ext_name2 = psf.coordinate_columns()
            new_column_data.append(psf.render(data[ext_name1][i], data[ext_name2][i]))
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

class psf:
    def __init__(self, psfFileName):
        self.psfFileName = psfFileName

    def render(self):
        pass

    def add_noise(self):
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

    def render(self, x, y):
        psfex_model = self.psfex
        return psfex_model.get_rec(y, x)

    def coordinate_columns(self):
        return 'XWIN_IMAGE', 'YWIN_IMAGE'

    def nameColumn(self):
        return 'VIGNET_PSFEX'

class webbpsf(psf):
    def __init__(self, psfFileName):
        super().__init__(psfFileName)

    def render(self):
        pass

    def nameColumn(self):
        return 'VIGNET_WEBBPSF'

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

catalog_object = catalog('new_file.fits')
psfex_object = epsfex('working/psfex-output/jw01727116001_04101_00001_nrca3_cal/jw01727116001_04101_00001_nrca3_cal_starcat.psf')
piff_object = piff_psf('/home/eddieberman/research/mcclearygroup//mock_data/mosaics/COSMOS2020_sims/piff-output/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci/mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_sci.piff')
catalog_object.augment(piff_object)
