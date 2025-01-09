export CONFIGDIR='/Users/j.mccleary/Research/jwst_cosmos/astro_config'
export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/mock_data/horizAGN_sims'
export OUTDIR='/Users/j.mccleary/Research/jwst_cosmos/mock_data/horizAGN_sims/working'
python get_jwst_psf.py -basedir $DATADIR -outdir $OUTDIR -configdir $CONFIGDIR


'''
usage: get_jwst_psf.py [-h] [-outdir OUTDIR] [-configdir CONFIGDIR] [-truthstars TRUTHSTARS] [--overwrite] [--vb] images [images ...]
'''

export CODEDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf'
export CONFIGDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config'
export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Apr2023/visit078/f150w'
export OUTDIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Apr2023/visit078/f150w/working'

python $CODEDIR/get_jwst_psf.py -config $CODEDIR/'configs/apr2023_config.yaml'  $DATADIR/jw01727078001_04101_00003_nrc*_cal.fits


export CODEDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf'
export CONFIGDIR='/Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config'
export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/full_i2d'
export OUTDIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/working2'

python $CODEDIR/get_jwst_psf.py -config $CODEDIR/'configs/box_cutter.yaml'  $DATADIR/*30mas*i2d.fits

'''
usage: master_psf_diagnostic.py basedir star_cat [-h] [-min_snr MIN_SNR]
            [-pix_scale PIX_SCALE] [-outdir OUTDIR] [-psfex_name PSFEX_NAME]
            [-im_name IM_NAME] [-piff_name PIFF_NAME] [--epsfex]
            [--gpsfex] [--piff] [--noisefree] [--verbose]
'''

python $CODEDIR/master_psf_diagnostic.py $DATADIR \
        /Users/j.mccleary/Research/jwst_cosmos/real_data/Apr2023/jw01727128001/working/f444w_b5_augmented_psf_starcat.fits   \
        -min_snr 10 -pix_scale 0.03 -outdir $OUTDIR/psf_diagnostic_plots/mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_i2d/ \
        -vignet_size 261\
        -webb_model /Users/j.mccleary/Research/jwst_cosmos/WebbPSF_models/f444w_b5_WebbPSF.fits


        -single_model /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/JWST-PSFEx_out_f150w_JAN_psf_v3.3_bg-sub.psf

python $CODEDIR/master_psf_diagnostic.py $DATADIR \
        $OUTDIR/mosaic_nircam_f277w_COSMOS-Web_30mas_v0_1_pex_stars.fits \
        -min_snr 10 -pix_scale 0.03 -outdir $OUTDIR/psf_diagnostic_plots/mosaic_nircam_f277w_COSMOS-Web_30mas_v0_1_i2d/ \
        -psfex_model $OUTDIR/psfex-output/mosaic_nircam_f277w_COSMOS-Web_30mas_v0_1_i2d/mosaic_nircam_f277w_COSMOS-Web_30mas_v0_1_starcat.psf --epsfex -vignet_size 75\

        -webb_model /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/PSF_F277W_webbpsf-0.03-COSMOS_NIRCam-jit001_v1.1.1_rot2.fits  -single_model /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/JWST-PSFEx_out_f277w_JAN_psf_v3.3_bg-sub.psf

python $CODEDIR/master_psf_diagnostic.py $DATADIR \
        $OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_pex_stars.fits \
        -min_snr 100 -pix_scale 0.03 -outdir $OUTDIR/psf_diagnostic_plots/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_i2d/ \
        -psfex_model $OUTDIR/psfex-output/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_i2d/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_starcat.psf --epsfex -vignet_size 131

        \
        -webb_model /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/PSF_F444W_webbpsf-0.03-COSMOS_NIRCam-jit001_v1.1.1_rot2.fits  -single_model /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/JWST-PSFEx_out_f444w_JAN_psf_v3.3_bg-sub.psf

python $CODEDIR/master_psf_diagnostic.py $DATADIR \
        $OUTDIR/jw01727128001_04101_00004_nrcalong_cal_pex_stars.fits \
        -min_snr 5 -pix_scale 0.06 -outdir $OUTDIR/psf_diagnostic_plots/jw01727128001_04101_00004_nrcalong_cal \
        -psfex_model $OUTDIR/psfex-output/jw01727128001_04101_00004_nrcalong_cal/jw01727128001_04101_00004_nrcalong_cal_starcat.psf \
        --epsfex -vignet_size 75\
        -webb_model $DATADIR/jw01727128001_04101_00004_nrcalong_cal_WebbPSF_rot.fits

## Run sextractor
export DATADIR='/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023'

sex "$DATADIR/full_i2d/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_i2d.fits[1]" \
    -WEIGHT_TYPE NONE -CATALOG_NAME $OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_cat.fits \
    -CHECKIMAGE_NAME $OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_sub.fits,$OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_bg.fits,$OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_bgRMS.fits \
    -CHECKIMAGE_TYPE -BACKGROUND,BACKGROUND,BACKGROUND_RMS -PARAMETERS_NAME $CONFIGDIR/sextractor.param \
    -STARNNW_NAME $CONFIGDIR/default.nnw -FILTER_NAME $CONFIGDIR/gauss_3.0_5x5.conv \
    -SEEING_FWHM 0.058 -c $CONFIGDIR/sextractor.config


sex "$DATADIR/full_i2d/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_i2d.fits[1]" \
    -WEIGHT_TYPE NONE -CATALOG_NAME $OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_cat.fits \
    -CHECKIMAGE_NAME $OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_sub.fits,$OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_bg.fits,$OUTDIR/mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_bgRMS.fits \
    -CHECKIMAGE_TYPE -BACKGROUND,BACKGROUND,BACKGROUND_RMS -PARAMETERS_NAME $CONFIGDIR/sextractor.param \
    -STARNNW_NAME $CONFIGDIR/default.nnw -FILTER_NAME $CONFIGDIR/gauss_3.0_5x5.conv \
    -SEEING_FWHM 0.058 -c $CONFIGDIR/sextractor.config


piffify /Users/j.mccleary/Research/jwst_cosmos/cweb_psf/astro_config/piff.config \
input.image_file_name=$DATADIR/mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_i2d.fits \
input.cat_file_name=$OUTDIR/mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_starcat.fits \
output.file_name=$OUTDIR/mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_sci.piff \
output.dir=$OUTDIR/piff-output/mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_i2d \
input.ra=10.007754753333334 input.dec=2.200973097

    type: SizeMag

    initial_select:
        # These are the default values
        type: SmallBright
        bright_fraction: 0.8
        small_fraction: 0.2
        locus_fraction: 0.5
        max_spread: 0.1

    fit_order: 2

    purity: 0.01

    num_iter: 3

    reserve_frac: 0.15

    fig, ax = plt.subplots(1,1, tight_layout=True, figsize=(10,8))

    ax.plot(im_cat['MAG_AUTO'], im_cat['FWHM_IMAGE'], '.', \
            label='all objects',markersize=2)
    ax.plot(star_cat['MAG_AUTO'], star_cat['FWHM_IMAGE'], '*', \
            label='selected stars',markersize=2)

    ax.set_xlabel('MAG_AUTO (zpt=30.0)', fontsize=14)
    ax.set_ylabel('FWHM_IMAGE (pixels)', fontsize=14)
    ax.set_ylim(-0.05, 40)
    ax.set_xlim(16,35)

    ax.set_title(f'{str(filter_name)} 30 mas SExtractor catalog', fontsize=14)
    ax.legend(markerscale=5, fontsize=14)
    fig.savefig(plot_file_name)


###
### Marko command; the config is his config
###
export $SEcat='/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/COSMOS-Web-PSFEx-SE-Cat_f277w_JAN.fits'

$SEcat -c $CONFIGDIR/marko_psfex.config -PSF_SIZE 261,261 -PSFVAR_DEGREES 2 -PSFVAR_NSNAP 4 -PSF_SAMPLING 1.0 -PSF_RECENTER Y -BASIS_TYPE PIXEL_AUTO -PSF_SUFFIX '_psf_'$ver'.psf' -PSF_ACCURACY 0.05 \
# -SAMPLEVAR_TYPE NONE -SAMPLE_VARIABILITY 0.2 -SAMPLE_AUTOSELECT y -SAMPLE_FWHMRANGE 3.5,5.0 -SAMPLE_MINSN 500 -SAMPLE_MAXELLIP 0.3 -BASIS_NUMBER 20 -PHOTFLUX_KEY "FLUX_APER(1)" -PHOTFLUXERR_KEY "FLUXERR_APER(1)" \
# -CHECKIMAGE_TYPE CHI,SAMPLES,RESIDUALS,SNAPSHOTS -CHECKIMAGE_NAME $outputcat'_ChckImg_chi2.fits',$outputcat'_ChckImg_samp.fits',$outputcat'_ChckImg_res.fits',$outputcat'_ChckImg_snap.fits' \
# -CHECKPLOT_TYPE SELECTION_FWHM,CHI2,COUNTS -CHECKPLOT_NAME $outputcat'_selfwhm',$outputcat'_chi2',$outputcat'_counts' \
# -OUTCAT_TYPE FITS_LDAC -OUTCAT_NAME $outputcat'.cat'

Oh shit, he supplied the whole thing!

export SEcat='/Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/mshuntov_psfex/COSMOS-Web-PSFEx-SE-Cat_f150w_JAN.fits'

psfex $SEcat -c $CONFIGDIR/psfex_marko.conf -PSF_SIZE 201,201 -PSFVAR_DEGREES 2 \
-PSFVAR_NSNAP 4 -PSF_SAMPLING 1.0 -PSF_RECENTER Y -BASIS_TYPE PIXEL_AUTO -PSF_SUFFIX -PSF_ACCURACY 0.05 \
-PSF_DIR /Users/j.mccleary/Research/jwst_cosmos/real_data/Jan2023/working/psfex-output/mosaic_nircam_f150\
-SAMPLEVAR_TYPE NONE -SAMPLE_VARIABILITY 0.2 -SAMPLE_AUTOSELECT Y -SAMPLE_FWHMRANGE 2.0,3.0 -SAMPLE_MINSN 300 \
-SAMPLE_MAXELLIP 0.3 -BASIS_NUMBER 20 -PHOTFLUX_KEY "FLUX_APER(1)" -PHOTFLUXERR_KEY "FLUXERR_APER(1)"
