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

