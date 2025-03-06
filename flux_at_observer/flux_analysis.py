import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Read the CSV file
data = pd.read_csv('fluxmap_data.csv')

# Get unique theta values and calculate mean fraction for each theta
theta_means = data.groupby('theta')['fraction'].mean()
theta_values = theta_means.index.values
fraction_means = theta_means.values

# Define cosine function to fit
def cosine_func(x, a, b, c):
    # Convert degrees to radians for numpy
    return a * np.cos(np.deg2rad(b * x)) + c

# Fit the cosine function to the data
popt, _ = curve_fit(cosine_func, theta_values, fraction_means)

# Generate smooth curve for plotting
theta_smooth = np.linspace(0, 90, 1000)
fit_curve = cosine_func(theta_smooth, *popt)

# Create the plot
plt.figure(figsize=(10, 6))
plt.scatter(theta_values, fraction_means, color='blue', alpha=0.5, label='Data (mean values)')
plt.plot(theta_smooth, fit_curve, 'r-', label=f'Fit: {popt[0]:.3f}*cos({popt[1]:.3f}θ) + {popt[2]:.3f}')

plt.xlabel('θ (degrees)')
plt.ylabel('Fraction')
plt.title('Flux Fraction vs Theta with Cosine Fit')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()