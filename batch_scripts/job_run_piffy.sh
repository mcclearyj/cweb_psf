#!/bin/sh
#SBATCH -t 6:59:59
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --mem-per-cpu=10g
#SBATCH --partition=short
#SBATCH -J real_piff
#SBATCH -v
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmac.ftw@gmail.com
#SBATCH -o piff_run.out


source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh

conda activate cweb_psf


export DATADIR='/work/mccleary_group/berman.ed/Apr_Data/v0_002'
export CATDIR='/work/mccleary_group/berman.ed/working'
export CONFIGDIR='/work/mccleary_group/jwst_cosmos/cweb_psf/astro_config'

stamps=(A1 A5 A6 A10)

bands=(f115w f150w f277w f444w)

for stamp in "${stamps[@]}"; do
    for band in "${bands[@]}"; do

	IM_BASE="mosaic_nircam_${band}_COSMOS-Web_30mas_${stamp}_v0_002_i2d"
	OUTDIR=${DATADIR}/piff-output/${IM_BASE}
	
	ls "$CATDIR/${IM_BASE}_real_train_75_starcat.fits"

	## Define parent directory
	if [ ! -d "$OUTDIR" ]; then
	    mkdir -p "$OUTDIR"
	fi

	
	piffify $CONFIGDIR/piff.config \
		input.image_file_name=$DATADIR/${IM_BASE}.fits \
		input.ra=10.007754753333334 input.dec=2.200973097 \
		input.cat_file_name=$CATDIR/${IM_BASE}_real_train_75_starcat.fits \
		output.file_name=${OUTDIR}/${IM_DIR}.piff \
		output.dir=$OUTDIR -v 2


   done
done

