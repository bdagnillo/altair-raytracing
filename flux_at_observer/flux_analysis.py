import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import warnings
import os
from datetime import datetime

# Check command line arguments
if len(sys.argv) < 2:
    print("Usage: python flux_analysis.py <csv_file>")
    sys.exit(1)

path = sys.argv[1]

# Extract metadata from file header
metadata = {}
header_line_count = 0
try:
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_line_count += 1
                # Extract key-value pairs from comments
                if ':' in line:
                    key, value = line[1:].strip().split(':', 1)  # Fixed typo: trip -> strip
                    metadata[key.strip()] = value.strip()
            else:
                break
except FileNotFoundError:
    print(f'File not found: {path}')
    sys.exit(1)

# Read the CSV file, skipping header comments
try:
    # First read the entire file as text to properly filter comments
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Filter out comments at beginning and end
    data_lines = [line for line in lines if not line.startswith('#')]
    
    # Write filtered data to a temporary file
    temp_path = path + '.temp'
    with open(temp_path, 'w') as f:
        f.write(''.join(data_lines))
    
    # Read the filtered data
    data = pd.read_csv(temp_path)
    os.remove(temp_path)  # Clean up temp file
    
    # Force conversion to float
    data['theta'] = pd.to_numeric(data['theta'])
    data['phi'] = pd.to_numeric(data['phi'])
    data['fraction'] = pd.to_numeric(data['fraction'])
    
except Exception as e:
    print(f'Error reading CSV data: {e}')
    sys.exit(1)

# Add 2D heatmap plot function
def plot_2d_heatmap(data, metadata):
    plt.figure(figsize=(10, 8))
    
    # Create pivot table for 2D plotting
    pivot = data.pivot(index='theta', columns='phi', values='fraction')
    
    # Create heatmap
    im = plt.imshow(pivot, 
                    aspect='auto',
                    origin='lower',
                    extent=[0, 360, 0, 90],
                    interpolation='nearest',
                    cmap='viridis')
    
    # Add colorbar
    cbar = plt.colorbar(im)
    cbar.set_label('Fraction of rays detected')
    
    # Add title with metadata
    title = 'Detector Flux Map'
    if 'Number of rays' in metadata:
        title += f" ({metadata['Number of rays']})"
    if 'BRDF Model' in metadata:
        title += f"\n{metadata['BRDF Model']}"
    plt.title(title)
    
    plt.xlabel('φ (degrees)')
    plt.ylabel('θ (degrees)')
    plt.grid(True)
    
    return plt.gcf()

# Get unique theta values and calculate mean fraction for each theta
grouped = data.groupby('theta')
theta_means = grouped['fraction'].mean()
theta_values = np.array(theta_means.index.tolist(), dtype=float)
fraction_means = np.array(theta_means.values, dtype=float)

# Calculate standard errors for error bars
theta_std = grouped['fraction'].std().fillna(0.001)
theta_counts = grouped.size()
theta_errors = theta_std / np.sqrt(theta_counts)

# Define cosine function to fit
def cosine_func(x, a, b, c):
    # Convert degrees to radians for numpy
    return a * np.cos(np.deg2rad(b * x)) + c

# Attempt curve fitting with error handling
try:
    # Initial guess based on data
    p0 = [
        (np.max(fraction_means) - np.min(fraction_means))/2,  # amplitude
        1.0,  # frequency
        np.mean(fraction_means)  # vertical offset
    ]
    
    popt, pcov = curve_fit(cosine_func, theta_values, fraction_means, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    fit_label = f'Fit: {popt[0]:.3f}*cos({popt[1]:.3f}θ) + {popt[2]:.3f}'
    
except Exception as e:
    print(f"Fit error: {e}")
    # Use simple approximation
    popt = [np.mean(fraction_means)/2, 1.0, np.mean(fraction_means)/2]
    perr = [0, 0, 0]
    fit_label = 'Approximate fit: cos(θ) scaled'

# Generate smooth curve for plotting
theta_smooth = np.linspace(min(theta_values), max(theta_values), 1000)
fit_curve = cosine_func(theta_smooth, *popt)

# Create the plot with theta analysis
theta_fig = plt.figure(figsize=(10, 6))  # Store figure handle

# Plot with error bars
plt.errorbar(theta_values, fraction_means, yerr=theta_errors, fmt='o', color='blue', 
             alpha=0.5, label='Data points', capsize=5, elinewidth=1, markeredgewidth=1)

plt.plot(theta_smooth, fit_curve, 'r-', label=fit_label)

# Safely calculate step sizes
if len(theta_values) > 1:
    sorted_theta = np.sort(theta_values)
    theta_step = np.mean(np.diff(sorted_theta))
else:
    theta_step = 0

if len(theta_values) > 0:
    first_theta = theta_values[0]
    phi_values = data[data['theta'] == first_theta]['phi'].unique()
    phi_step = 360 / len(phi_values) if len(phi_values) > 0 else 0
else:
    phi_step = 0

# Calculate R-squared
residuals = fraction_means - cosine_func(theta_values, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((fraction_means - np.mean(fraction_means))**2)
r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

# Add metadata to title
title = 'Flux Fraction vs Theta with Cosine Fit'
if 'Number of rays' in metadata:
    title += f" ({metadata['Number of rays']})"
    
# Add geometry info 
subtitle = []
if 'Detector dimensions' in metadata:
    subtitle.append(f"Detector: {metadata['Detector dimensions']}")
if 'Sphere inner radius' in metadata and 'Sphere outer radius' in metadata:
    subtitle.append(f"Sphere: {metadata['Sphere inner radius']}-{metadata['Sphere outer radius']}")
if 'Gaussian roughness' in metadata:
    subtitle.append(f"Roughness: {metadata['Gaussian roughness']}")

# Add subtitle if we have any info
if subtitle:
    title += f"\n{', '.join(subtitle)}"
    
# Add step sizes and fit quality
title += f"\nSteps: Δθ = {theta_step:.1f}°, Δφ = {phi_step:.1f}°, R² = {r_squared:.4f}"

plt.xlabel('θ (degrees)')
plt.ylabel('Fraction')
plt.title(title)
plt.legend()
plt.grid(True)

# Save plot
output_filename = f"{os.path.splitext(os.path.basename(path))[0]}_plot.png"
# plt.savefig(output_filename, dpi=300)
# print(f"Plot saved as {output_filename}")

# Print fit parameters
print(f"Fit parameters: a={popt[0]:.5f}, b={popt[1]:.5f}, c={popt[2]:.5f}")
print(f"R-squared value: {r_squared:.5f}")

# Show the plot
plt.show()

# Create 2D heatmap
heatmap_fig = plot_2d_heatmap(data, metadata)

# Save both plots
base_filename = os.path.splitext(os.path.basename(path))[0]
theta_fig.savefig(f"{base_filename}_theta_analysis.png", dpi=300, bbox_inches='tight')
heatmap_fig.savefig(f"{base_filename}_heatmap.png", dpi=300, bbox_inches='tight')
print(f"Plots saved as {base_filename}_theta_analysis.png and {base_filename}_heatmap.png")

# Show all plots at once
plt.show()

