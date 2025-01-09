import os
import json
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import rc,rcParams
rc('font',**{'family':'serif'})
rc('text', usetex=True)


def aggregate_stats(base_dir):
    # Define the data structure for aggregation
    stats = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    # Recursively list files
    pattern = '**/*_diagnostics_results.json'
    files = glob.glob(os.path.join(base_dir, pattern), recursive=True)

    # Read JSON files and extract required information
    for file in files:
        with open(file, 'r') as f:
            data = json.load(f)

        # Parse the file path to get bandpass, tile, and PSF type
        parts = file.split(os.sep)
        bandpass = [part for part in parts if part.upper() in ['F115W', 'F150W', 'F277W', 'F444W']][0].upper()
        tile = [part for part in parts if 'B' in part][0].upper()

        # Iterate over PSF types in the JSON file
        for psf_type, contents in data.items():
            # Navigate to the 'chisq_resids' block and extract 'reduced_chi_square'
            reduced_chi_square = contents.get('chisq_resids', {}).get('reduced_chi_square', None)
            if reduced_chi_square is not None:
                # Aggregate statistics
                stats[bandpass][tile][psf_type].append(reduced_chi_square)


    holder = {
        'F115W': {
            'single': [],
            'PSFEx': []
        },
        'F150W': {
            'single': [],
            'PSFEx': []
        },
        'F277W': {
            'single': [],
            'PSFEx': []
        },
        'F444W': {
            'single': [],
            'PSFEx': []
        }
    }

    # Now this json file has everything. Let's go by keys
    bandpasses = ['F115W', 'F150W', 'F277W', 'F444W']
    for bandpass in bandpasses:
        tiles = stats[bandpass].keys()
        for tile in tiles:
            this_fit = stats[bandpass][tile]
            holder[bandpass]['single'].extend(
                this_fit['single']
            )
            holder[bandpass]['PSFEx'].extend(
                this_fit['psfex']
            )

    # Define a stats holder
    single_stats_holder = {
        'F115W': {
            'single': {'mean': None, 'std': None},
            'PSFEx': {'mean': None, 'std': None}
        },
        'F150W': {
            'single': {'mean': None, 'std': None},
            'PSFEx': {'mean': None, 'std': None}
        },
        'F277W': {
            'single': {'mean': None, 'std': None},
            'PSFEx': {'mean': None, 'std': None}
        },
        'F444W': {
            'single': {'mean': None, 'std': None},
            'PSFEx': {'mean': None, 'std': None}
        }
    }

    for bandpass in bandpasses:
        this_holder = holder[bandpass]
        # Single
        stats_holder[bandpass]['single']['mean'] = \
            np.mean(this_holder['single'])
        stats_holder[bandpass]['single']['std'] = \
            np.std(this_holder['single'])
        # PSFEx
        stats_holder[bandpass]['PSFEx']['mean'] = \
            np.mean(this_holder['PSFEx'])
        stats_holder[bandpass]['PSFEx']['std'] = \
            np.std(this_holder['PSFEx'])

    single_stats = {'mean': [], 'std': []}
    psfex_stats = {'mean': [], 'std': []}


    for bandpass in bandpasses:
        this_holder = holder[bandpass]
        # Single
        single_stats['mean'].append(np.mean(this_holder['single']))
        single_stats['std'].append(np.std(this_holder['single']))
        # PSFEx
        psfex_stats['mean'].append(np.mean(this_holder['PSFEx']))
        psfex_stats['std'].append(np.std(this_holder['PSFEx']))

    return single_stats, psfex_stats

def plot_results():

    plt.rcParams.update({'xtick.major.size': 8})
    plt.rcParams.update({'xtick.major.width': 1.3})
    plt.rcParams.update({'xtick.minor.visible': True})
    plt.rcParams.update({'xtick.minor.width': 1.1})
    plt.rcParams.update({'xtick.minor.size': 6})
    plt.rcParams.update({'xtick.direction': 'out'})
    plt.rcParams.update({'ytick.major.width': 1.3})
    plt.rcParams.update({'ytick.major.size': 8})
    plt.rcParams.update({'ytick.minor.visible': True})
    plt.rcParams.update({'ytick.minor.width': 1.1})
    plt.rcParams.update({'ytick.minor.size':6})
    plt.rcParams.update({'xtick.labelsize': 22})
    plt.rcParams.update({'ytick.labelsize': 22})

    bandpasses = ['F115W', 'F150W', 'F277W', 'F444W']


    fig, ax = plt.subplots(figsize=(15, 7), linewidth=2)
    positions = range(4)
    ax.set_yscale('log')

    # Plotting PSFex data
    ax.errorbar(positions, psfex_stats['mean'], yerr=psfex_stats['std'],
                fmt='o', label='PSFex', capsize=5, ms=8,
                capthick=1.1, color='red', lw=2)

    # Plotting single (Marko-format) data with a slight x shift for better visibility
    ax.errorbar([p + 0.2 for p in positions], single_stats['mean'],
                yerr=single_stats['std'],
                fmt='o', label='Single', capsize=5, ms=8, capthick=1.1,
                color='blue', lw=2)

    # Customization
    ax.set_title('Comparison of PSFex model' + r'reduced $\chi^2$', fontsize=25)
    ax.set_xlabel('Filter', fontsize=25)
    ax.set_ylabel('Average ' + r'$\chi^2$', fontsize=25)
    ax.set_xticks([p + 0.2 for p in positions])
    ax.set_xticklabels(bandpasses, rotation=45, ha="right")

    ax.legend(fontsize=22)
    plt.grid(True, which='both', linestyle='--', linewidth=1)
    plt.tight_layout()
    fig.savefig('reduced_chi2_comparison.png')


def main():

    base_dir = '/n23data1/mccleary/real_data/Jan2024'

    single_stats, psfex_stats = aggregate_stats(base_dir, which_stat)

    plot_results(psfex_stats, single_stats, which_stat)
