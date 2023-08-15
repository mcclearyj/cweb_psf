import glob
from astropy.io import fits
from astropy.table import Table, vstack
from catalogaugmenter import catalog, psf
from catalogaugmenter import webb_psf, epsfex, shopt, piff_psf 
#from catalogplotter import ResidPlots
import os
import re 
from datetime import datetime, timedelta

#Make augmented catalogs with columns for each psf fitter than use the plotter with these new catalogs

#ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
ims = glob.glob('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/jw0*cal.fits')
print(len(ims))
f115_a1_list=[]
f115_a2_list=[]
f115_a3_list=[]
f115_a4_list=[]
f115_b1_list=[]
f115_b2_list=[]
f115_b3_list=[]
f115_b4_list=[]
f150_a1_list=[]
f150_a2_list=[]
f150_a3_list=[]
f150_a4_list=[]
f150_b1_list=[]
f150_b2_list=[]
f150_b3_list=[]
f150_b4_list=[]
f277_a5_list=[]
f277_b5_list=[]
f444_a5_list=[]
f444_b5_list=[]

for im in ims:
    imhead = fits.getheader(im, ext=0)
    filt = imhead['FILTER']
    det = imhead['DETECTOR'].replace('LONG', '5')
    print(det, filt)
    #string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
    #string2 = string+'_train_starcat.psf'
    #print(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)
    if (det =='NRCA1') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_a1_list.append(catalog_obj)
            #f115_a1_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA1') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_a1_list.append(catalog_obj)
            #f150_a1_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA2') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_a2_list.append(catalog_obj)
            #f115_a2_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA2') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            #print(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_a2_list.append(catalog_obj)
        except:
            pass
        #f150_a2_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCA3') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_a3_list.append(catalog_obj)
            #f115_a3_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA3') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_a3_list.append(catalog_obj)
            #f150_a3_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA4') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_a4_list.append(catalog_obj)
            #f115_a4_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA4') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_a4_list.append(catalog_obj)
            #f150_a4_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA5') & (filt == 'F277W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f277_a5_list.append(catalog_obj)
            #f277_a5_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCA5') & (filt == 'F444W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f444_a5_list.append(catalog_obj)
            #f444_a5_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB1') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b1_list.append(catalog_obj)
            #f115_b1_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB1') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b1_list.append(catalog_obj)
            #f150_b1_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB2') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b2_list.append(catalog_obj)
            #f115_b2_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB2') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b2_list.append(catalog_obj)
            #f150_b2_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB3') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b3_list.append(catalog_obj)
            #f115_b3_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB3') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b3_list.append(catalog_obj)
            #f150_b3_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB4') & (filt == 'F115W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f115_b4_list.append(catalog_obj)
        except:
            pass
        #f115_b4_list.append(fits.getheader(im, ext=0)['DATE'])
    elif (det =='NRCB4') & (filt == 'F150W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f150_b4_list.append(catalog_obj)
            #f150_b4_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB5') & (filt == 'F277W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f277_b5_list.append(catalog_obj)
            #f277_b5_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass
    elif (det =='NRCB5') & (filt == 'F444W'):
        try:
            catalog_obj = catalog(im.replace('single_exposures', 'working').replace('cal.fits' ,'cal_valid_starcat.fits' ))
            catalog_obj.augment(webb_psf(im.replace('_cal', '_cal_WebbPSF')))
            string = im.replace('.fits', '').replace('/home/eddieberman/research/mcclearygroup/cweb_psf/single_exposures/', '')
            string2 = string+'_train_starcat.psf'
            catalog_obj.augment(epsfex(im.replace('single_exposures/', 'working/psfex-output/').replace('_cal.fits', '_cal/')+string2)) #working/psfex-output/
            f444_b5_list.append(catalog_obj)
            #f444_b5_list.append(fits.getheader(im, ext=0)['DATE'])
        except:
            pass

print(len(f115_a1_list))
print(len(f115_a2_list))
print(len(f115_a3_list))
print(len(f115_a4_list))
print(len(f115_b1_list))
print(len(f115_b2_list))
print(len(f115_b3_list))
print(len(f115_b4_list))
print(len(f150_a1_list))
print(len(f150_a2_list))
print(len(f150_a3_list))
print(len(f150_a4_list))
print(len(f150_b1_list))
print(len(f150_b2_list))
print(len(f150_b3_list))
print(len(f150_b4_list))
print(len(f277_a5_list))
print(len(f277_b5_list))
print(len(f444_a5_list))
print(len(f444_b5_list))


#cat = f444_b5_list[0]
#for catalog in f444_b5_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f444_b5_thirty_visits_valid_starcat.fits')


#cat = f444_a5_list[0]
#for catalog in f444_a5_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f444_a5_thirty_visits_valid_psf_starcat.fits')

def process_batch(batch):
    cat = batch[0]
    for catalog in batch[1:]:
        try:
            cat.concatenate_catalogs(catalog)
        except:
            pass
    return cat

# Split f277_b5_list into batches
batch_size = 10
batches = [f444_a5_list[i:i + batch_size] for i in range(0, len(f444_a5_list), batch_size)]

saved_cats = []

for idx, batch in enumerate(batches):
    cat = process_batch(batch)
    filename = f'a5_temporary_batch_{idx}.fits'
    cat.save_new(filename)
    saved_cats.append(filename)

batches = [f444_b5_list[i:i + batch_size] for i in range(0, len(f444_b5_list), batch_size)]

saved_cats = []

for idx, batch in enumerate(batches):
    cat = process_batch(batch)
    filename = f'b5_temporary_batch_{idx}.fits'
    cat.save_new(filename)
    saved_cats.append(filename)


# Once all batches are processed, combine the saved batch files
#combined_cat = catalog(saved_cats[0])
'''
for cat in saved_cats[1:]:
    try:
        combined_cat.concatenate_catalogs(catalog(cat))
    except:
        pass

combined_cat.save_new('f277_b5_thirty_visits_valid_psf_starcat.fits')
'''
#cat = f277_b5_list[0]
#for catalog in f277_b5_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f277_b5_thirty_visits_valid_psf_starcat.fits')



'''
cat = f277_a5_list[0]
for catalog in f277_a5_list[1:]:
    cat.concatenate_catalogs(catalog)
cat.save_new('f277_a5_thirty_visits_valid_psf_starcat.fits')
'''
#cat = f150_b4_list[0]
#for catalog in f150_b4_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_b4_thirty_visits_valid_psf_starcat.fits')

#cat = f115_b4_list[0]
#for catalog in f115_b4_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f115_b4_thirty_visits_valid_psf_starcat.fits')

#cat = f150_b3_list[0]
#for catalog in f150_b3_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_b3_thirty_visits_valid_psf_starcat.fits')

#cat = f115_b3_list[0]
#for catalog in f115_b3_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f115_b3_thirty_visits_valid_psf_starcat.fits')

#cat = f150_b2_list[0]
#for catalog in f150_b2_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_b2_thirty_visits_valid_psf_starcat.fits')

#cat = f115_b2_list[0]
#for catalog in f115_b2_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f115_b2_thirty_visits_valid_psf_starcat.fits')

#cat = f150_b1_list[0]
#for catalog in f150_b1_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_b1_thirty_visits_valid_psf_starcat.fits')

#cat = f115_b1_list[0]
#for catalog in f115_b1_list[1:]:
#    cat.concatenate_catalogs(catalog)
#cat.save_new('f115_b1_thirty_visits_valid_psf_starcat.fits')

#cat = f150_a4_list[0]
#for catalog in f150_a4_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_a4_thirty_visits_valid_psf_starcat.fits')

#cat = f115_a4_list[0]
#for catalog in f115_a4_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f115_a4_thirty_visits_valid_psf_starcat.fits')

#cat = f150_a3_list[0]
#for catalog in f150_a3_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f150_a3_thirty_visits_valid_psf_starcat.fits')

#cat = f115_a3_list[0]
#for catalog in f115_a3_list[1:]:
#    try:
#        cat.concatenate_catalogs(catalog)
#    except:
#        pass
#cat.save_new('f115_a3_thirty_visits_valid_psf_starcat.fits')
'''
cat = f150_a2_list[0]
for catalog in f150_a2_list[1:]:
    try:
        cat.concatenate_catalogs(catalog)
    except:
        pass
cat.save_new('f150_a2_thirty_visits_valid_psf_starcat.fits')

cat = f150_a1_list[0]
for catalog in f150_a1_list[1:]:
    try:
        cat.concatenate_catalogs(catalog)
    except:
        pass
cat.save_new('f150_a1_thirty_visits_valid_psf_starcat.fits')
'''
'''
cat = f115_a2_list[0]
for catalog in f115_a2_list[1:]:
    try:
        cat.concatenate_catalogs(catalog)
    except:
        pass
cat.save_new('f115_a2_thirty_visits_valid_psf_starcat.fits')

cat = f115_a1_list[0]
for catalog in f115_a1_list[1:]:
    try:
        cat.concatenate_catalogs(catalog)
    except:
        pass
cat.save_new('f115_a1_thirty_visits_valid_psf_starcat.fits')
'''
