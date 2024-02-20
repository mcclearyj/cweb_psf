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

    # Run for every star cat
    starcats = args.starcats

    for starcat in starcats:
        psf_diagnostic = PSFDiagnosticsRunner(starcat, config=config)
        psf_diagnostic.make_ouput_subdir()
        try:
            psf_diagnostic.run_all()
        except:
            pdb.set_trace()
    return 0
if __name__ == '__main__':

   args = parse_args()
   rc = main(args)

   if rc == 0:
      print('run_diagnostics.py successfully completed')
   else:
      print(f'run_diagnostics.py has failed with rc={rc}')
