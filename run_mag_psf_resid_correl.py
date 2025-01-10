import pdb
import time
from argparse import ArgumentParser

# Local imports
from src.mag_psf_resid_correl import MagPsfResidCorrel
from src import utils

def parse_args():
    """ Parse arguments, duh """
    parser = ArgumentParser()
    parser.add_argument(
        "hsm_cat_file", default=None,
        help="cweb_psf-format star and PSF HSM fit file"
    )
    parser.add_argument(
        "-config","-c", required=True, default=None,
        help="cweb_psf-format config file"
    )

    return parser.parse_args()

def main(args):
    """
    At some point, could implement config-checking here
    and perhaps file checking
    """

    hsm_cat_file = args.hsm_cat_file
    config_file = args.config
    run_config = utils.read_yaml(config_file)

    ## Grab allowed models
    model_map = ['psfex', 'webbpsf', 'single', 'mirage', 'piff', 'shopt']

    tstart = time.perf_counter()

    for model in model_map:
        if run_config['psf_models'].get(model, False):
            mprc = MagPsfResidCorrel(hsm_cat_file, run_config, model, nbins=13)
            mprc.run()

    tend = time.perf_counter()

    print(f"Total run time is {tend-tstart} s")

if __name__ == '__main__':

    args = parse_args()
    main(args)
