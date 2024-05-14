#!/bin/bash
#PBS -S /bin/sh
#PBS -N starcats_diagnostics
#PBS -o ./starcats_diagnostics.out 
#PBS -j oe 
#PBS -l nodes=1,walltime=3:00:00
#PBS -W depend=afterok:297901.cmaster

module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /n23data1/mccleary/miniconda3/bin/activate 
conda activate cweb_psf

export DATADIR="/n23data1/mccleary/real_data/Jan2024/"
export CODEDIR="/n23data1/mccleary/cweb_psf"
export CONFIG=${CODEDIR}/"configs/real_mosaics_psf_diagnostics_config.yaml"

export BANDPASSES="f115w f150w f277w f444w"

for bandpass in $BANDPASSES; do

    export STARCATS=(${DATADIR}/${bandpass}/"mosaic_nircam_*_COSMOS-Web_30mas_*_i2d_valid_starcat.fits")
    export 
    if [ ${#STARCATS[@]} -gt 0 ]; then
        echo ""
        echo "Input images: ${IMAGES[*]}"
        echo "Using config: $CONFIG"
        echo ""
	
	python ${CODEDIR}/run_diagnostics.py -config ${CONFIG} -outdir ${DATADIR}/${bandpass}/"psf-diagnostics" ${STARCATS[@]}
    else
        echo "No valid files found for ${bandpass}"
    fi

    echo "done with ${bandpass}"
done


echo "done"
