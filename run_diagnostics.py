#### run_psf_diagnostics.py:
####	takes an input LDAC-format star catalog
#### 	augmented with PSF vignettes (renderings) and creates a variety of
#### 	diagnostic figures, including chi2, shearsticks, flux residuals, etc.
####
#### syntax:
####	python run_psf_diagnostics.py -config [CONFIG_FILE] starcat1 [starcat2 starcat3...]
####

## Global imports
from argparse import ArgumentParser
import pdb
import re

## Local imports
from src.plotter import size_mag_plots
from src.utils import read_yaml, make_outdir, concatenate_catalogs
from src.box_cutter import BoxCutter
from src.psf_renderer import PSFRenderer
from src.psf_diagnostics_runner import PSFDiagnosticsRunner

def parse_args():

    parser = ArgumentParser()
    parser.add_argument('starcats', nargs='+',
			help = 'Augmented star catalogs to process (wildcards OK)')
    parser.add_argument('-config', default=None,
                        help = 'Configuration file for this script')
    parser.add_argument('-outdir', default=None,
                        help = 'Output directory')
    parser.add_argument('--vb', action='store_true', default=True,
                        help = 'Print detailed messages [does nothing for now]')

    return parser.parse_args()


def main(args):
    # Load config
    config = read_yaml(args.config)

    # Adds an outdir parameter to config if it was missing
    make_outdir(config, arg=args.outdir)

    # For convenience, concatenate individual star catalogs
    starcats = args.starcats
    
    # A little kludgy, but:
    bandpass = re.search(r"f(\d){3}w", starcats[0]).group().upper()
    combined_starcat = concatenate_catalogs(
        starcats, outname=f'combined_validation_starcats_{bandpass}.fits', hdu=2
    )

    # OK, now run on combined catalog 
    psf_diagnostic = PSFDiagnosticsRunner(combined_starcat, config=config)
    psf_diagnostic.run_all()
        
    return 0
if __name__ == '__main__':

   args = parse_args()
   rc = main(args)

   if rc == 0:
      print('run_diagnostics.py successfully completed')
   else:
      print(f'run_diagnostics.py has failed with rc={rc}')
