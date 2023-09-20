import numpy as np
import matplotlib.pyplot as plt

def polynomial_fit_and_plot(datasets):
    colors = ['blue', 'green', 'purple']
    
    for i, data in enumerate(datasets):
        x, y = data
        # Fit a 2nd degree polynomial to the data
        coefficients = np.polyfit(x, y, 2)
        polynomial = np.poly1d(coefficients)
        
        # Generate x values for plotting the polynomial curve
        x_poly = np.linspace(min(x), max(x), 100)
        y_poly = polynomial(x_poly)
        
        # Plotting the data points
        plt.scatter(x, y, color=colors[i], label=f'Data points {i + 1}')
        
        # Plotting the polynomial fit
        plt.plot(x_poly, y_poly, color=colors[i], linestyle='--', label=f'2nd Degree Fit {i + 1}')

    # Adding fancy touches to the plot
    plt.title('Execution Time vs Number of Pixels')
    plt.xlabel('Number of Pixels')
    plt.ylabel('Execution Time (s)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    
    plt.show()

# Three sets of fake data points
dataset1 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.5, 0.9, 3.1, 12.8, 49.7, 200.2])
dataset2 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.4, 0.8, 2.9, 11.7, 46.6, 180.1])
dataset3 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.6, 1.0, 3.3, 14.0, 51.9, 210.3])

polynomial_fit_and_plot([dataset1, dataset2, dataset3])

import numpy as np
import matplotlib.pyplot as plt

def linear_fit_and_plot(datasets):
    colors = ['blue', 'green', 'purple']
    
    for i, data in enumerate(datasets):
        x, y = data
        # Fit a 1st degree polynomial (linear) to the data
        coefficients = np.polyfit(x, y, 1)
        linear_polynomial = np.poly1d(coefficients)
        
        # Generate x values for plotting the linear curve
        x_poly = np.linspace(min(x), max(x), 100)
        y_poly = linear_polynomial(x_poly)
        
        # Plotting the data points
        plt.scatter(x, y, color=colors[i], label=f'Data points {i + 1}')
        
        # Plotting the linear fit
        plt.plot(x_poly, y_poly, color=colors[i], linestyle='--', label=f'Linear Fit {i + 1}')

    # Adding fancy touches to the plot
    plt.title('Execution Time vs Number of Pixels')
    plt.xlabel('Number of Pixels')
    plt.ylabel('Execution Time (s)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()
    
    plt.show()

# Three sets of fake data points
dataset1 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.5, 0.9, 3.1, 12.8, 49.7, 200.2])
dataset2 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.4, 0.8, 2.9, 11.7, 46.6, 180.1])
dataset3 = ([100000, 200000, 400000, 800000, 1600000, 3200000],
            [0.6, 1.0, 3.3, 14.0, 51.9, 210.3])

linear_fit_and_plot([dataset1, dataset2, dataset3])

import numpy as np
import matplotlib.pyplot as plt

def plot_with_connections(datasets):
    colors = ['blue', 'green', 'purple']

    for i, data in enumerate(datasets):
        x, y = data

        # Plotting the data points
        plt.scatter(x, y, color=colors[i], label=f'Data points {i + 1}')

        # Connecting the data points with a line
        plt.plot(x, y, color=colors[i], linestyle='--', alpha=0.7)

    # Adding fancy touches to the plot
    plt.title('Execution Time vs Number of Pixels with Connections')
    plt.xlabel('Number of Pixels')
    plt.ylabel('Execution Time (s)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()

    plt.show()

# Three sets of fake data points, each with 3 data points
dataset1 = ([100000, 200000, 400000], [0.5, 0.9, 3.1])
dataset2 = ([100000, 200000, 400000], [0.4, 1.2, 3.2])
dataset3 = ([100000, 200000, 400000], [0.6, 1.1, 2.9])

plot_with_connections([dataset1, dataset2, dataset3])

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot_with_interpolation(all_data_sets):
    colors = ['blue', 'green', 'purple']
    filters = ['F1115W', 'F150W', 'F277W', 'F444W']
    methods = ['PIFF', 'ShOpt', 'PSFex']

    fig = plt.figure(figsize=(15, 5))
    spec = gridspec.GridSpec(ncols=4, nrows=1, figure=fig)

    ax1 = fig.add_subplot(spec[0])
    ax1.set_title(filters[0])

    for j, datasets in enumerate(all_data_sets):
        if j != 0:
            ax = fig.add_subplot(spec[j], sharey=ax1)
            ax.set_title(filters[j])
        else:
            ax = ax1

        for i, data in enumerate(datasets):
            x, y = data
            # Interpolation of degree 2
            z = np.polyfit(x, y, 2)
            p = np.poly1d(z)
            
            # Plot interpolated data
            x_interpolated = np.linspace(min(x), max(x), 500)
            y_interpolated = p(x_interpolated)
            ax.plot(x_interpolated, y_interpolated, color=colors[i], linestyle='-', alpha=0.7, label=methods[i])
            
            # Plot original data points
            ax.scatter(x, y, color=colors[i], s=30)
            
            if j == 0:
                ax.legend()
                
            if j != 0:
                plt.setp(ax.get_yticklabels(), visible=False)

        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.set_xlabel('Number of Pixels')

    ax1.set_ylabel('Execution Time (s)')
    fig.suptitle('Execution Time vs Number of Pixels with Interpolation of Degree 2')
    fig.tight_layout()
    plt.show()

# Four sets of data, each containing 3 datasets
data1 = [([100000, 200000, 400000], [0.5, 0.9, 3.1]),
         ([120000, 220000, 420000], [0.4, 1.2, 3.2]),
         ([130000, 230000, 440000], [0.6, 1.1, 2.9])]

data2 = [([105000, 205000, 405000], [0.7, 1.0, 3.3]),
         ([125000, 225000, 425000], [0.5, 1.3, 3.4]),
         ([135000, 235000, 445000], [0.8, 1.2, 3.0])]

data3 = [([110000, 210000, 410000], [0.6, 0.8, 3.2]),
         ([115000, 215000, 415000], [0.4, 1.1, 3.1]),
         ([125000, 225000, 425000], [0.7, 1.0, 2.8])]

data4 = [([102000, 202000, 402000], [0.9, 1.1, 3.4]),
         ([122000, 222000, 422000], [0.3, 1.4, 3.0]),
         ([132000, 232000, 442000], [0.8, 1.3, 3.1])]

plot_with_interpolation([data1, data2, data3, data4])


import numpy as np
import matplotlib.pyplot as plt

# Mock data
data = {
    'F115W': {
        'ShOpt': np.random.rand(10),
        'PIFF': np.random.rand(10),
        'PSFex': np.random.rand(10)
    },
    'F150W': {
        'ShOpt': np.random.rand(10),
        'PIFF': np.random.rand(10),
        'PSFex': np.random.rand(10)
    },
    'F277W': {
        'ShOpt': np.random.rand(10),
        'PIFF': np.random.rand(10),
        'PSFex': np.random.rand(10)
    },
    'F444W': {
        'ShOpt': np.random.rand(10),
        'PIFF': np.random.rand(10),
        'PSFex': np.random.rand(10)
    }
}

x = np.linspace(0, 10, 10)

fig = plt.figure(figsize=(18, 10))
grid = plt.GridSpec(3, 5, wspace=0.4, hspace=0.3)

for row, degree in enumerate([1, 2, 3]):
    for col, filter_name in enumerate(['F115W', 'F150W', 'F277W', 'F444W']):
        ax = fig.add_subplot(grid[row, col])
        for psf_fitter, color in zip(['ShOpt', 'PIFF', 'PSFex'], ['red', 'blue', 'green']):
            y = data[filter_name][psf_fitter]
            ax.scatter(x, y, color=color, label=psf_fitter)
            
            # Fit and plot polynomial
            coefficients = np.polyfit(x, y, degree)
            poly = np.poly1d(coefficients)
            ax.plot(x, poly(x), color=color, linestyle='--')
        
        ax.set_title(filter_name if row == 0 else "")
        ax.set_ylabel(f'Degree {degree}\nExecution Time (seconds)', labelpad=15 if col == 0 else 0)
        if row == 2:
            ax.set_xlabel('# Pixels', labelpad=15)
            
# Setting up legend in its own axis
legend_ax = fig.add_subplot(grid[:, 4])
handles, labels = ax.get_legend_handles_labels()
legend_ax.legend(handles, labels, loc='center', prop={'size': 14})
legend_ax.axis('off')

plt.show()




















import matplotlib.pyplot as plt
import numpy as np

# Set random seed for reproducibility
np.random.seed(0)

# Define polynomial degrees
degrees = [1, 2, 3]

# Define PSF fitters
psf_fitters = ['ShOpt', 'PSFex', 'PIFF']

# Define wavelengths
wavelengths = ['F115W', 'F150W', 'F277W', 'F444W']

# Generate random execution times for each PSF fitter for the given wavelengths and degrees
times = {fitter: {wavelength: np.random.rand(3) for wavelength in wavelengths} for fitter in psf_fitters}

# Create a plot
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(20, 5))

for i, wavelength in enumerate(wavelengths):
    ax = axs[i]
    
    # Plot the execution times for each PSF fitter
    for fitter in psf_fitters:
        ax.plot(degrees, times[fitter][wavelength], '-o', label=fitter)
    
    ax.set_title(f"Wavelength: {wavelength}")
    ax.set_xlabel('Polynomial Degree')
    ax.set_ylabel('Execution Time')
    ax.set_xticks(degrees)
    ax.legend()

plt.tight_layout()
plt.show()

