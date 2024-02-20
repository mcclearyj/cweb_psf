import numpy as np
from astropy.io import fits

filename = 'mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_starcat.copy.fits'

# Open the FITS file
with fits.open(filename, mode='update') as hdulist:
    
    # Get the FITS_REC object from the file
    fits_rec = hdulist[2].data
    
    # Create a new array of data with the same length as fits_rec
    new_col = np.random.rand(len(fits_rec), 31, 31)
    
    # Define the new column name and data format
    new_col_name = 'NEW_COLUMN'
    new_col_format = '31E,31E'
    
    # Add the new column to the FITS_REC object
    new_column = fits.Column(name=new_col_name, format=new_col_format, array=new_col)
    
    new_hdu = fits.BinTableHDU.from_column(new_column)
    hdulist[2] = new_hdu
    
    # Write the updated FITS file to disk
    hdulist.flush()

fits.BinTableHDU.from_columns(

for file in *.fits; do
    fitsheader $file | grep "STREHL" | awk -v fname="$file" '{print fname "\t" $0}' >> output.txt
done
