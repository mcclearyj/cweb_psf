#!/bin/bash
#SBATCH --job-name=starcats_f150w
#SBATCH --output=./f150w_ims.out
#SBATCH --error=./f150w_ims.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=j.mccleary@northeastern.edu

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data1/mccleary/real_images/20mas_resamp"
#export BANDPASSES="f115w f150w"
export BANDPASSES="f150w"

export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/configs/20mas_mosaics_config.yaml

for bandpass in $BANDPASSES; do
    IMAGES=(${DATADIR}/${bandpass}/visit_*_nircam_${bandpass}_COSMOS-Web_20mas_*_i2d.fits)

    if [ ${#IMAGES[@]} -gt 0 ]; then
        echo ""
        echo "Input images: ${IMAGES[*]}"
        echo "Using config: $CONFIG"
        echo ""

        python ${CODEDIR}/get_galaxy_cutouts.py -c ${CONFIG} -outdir "${DATADIR}/gal_cutout_20mas_working/${bandpass}" "${IMAGES[@]}"

    else
        echo "No valid files found for ${bandpass}"
    fi

    echo "done with ${bandpass}"
done
