import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import warnings
import os
from datetime import datetime

# Function to process a single file and return the data and metadata
def process_file(filepath):
    # Extract metadata from file header
    metadata = {}
    header_line_count = 0
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header_line_count += 1
                    # Extract key-value pairs from comments
                    if ':' in line:
                        key, value = line[1:].strip().split(':', 1)
                        metadata[key.strip()] = value.strip()
                else:
                    break
    except FileNotFoundError:
        print(f'File not found: {filepath}')
        return None, None

    # Read the CSV file, skipping header comments
    try:
        # First read the entire file as text to properly filter comments
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Filter out comments at beginning and end
        data_lines = [line for line in lines if not line.startswith('#')]
        
        # Write filtered data to a temporary file
        temp_path = filepath + '.temp'
        with open(temp_path, 'w') as f:
            f.write(''.join(data_lines))
        
        # Read the filtered data
        data = pd.read_csv(temp_path)
        os.remove(temp_path)  # Clean up temp file
        
        # Force conversion to float
        data['theta'] = pd.to_numeric(data['theta'])
        data['phi'] = pd.to_numeric(data['phi'])
        data['fraction'] = pd.to_numeric(data['fraction'])
        
        return data, metadata
        
    except Exception as e:
        print(f'Error reading CSV data from {filepath}: {e}')
        return None, None

# Define cosine function to fit
def cosine_func(x, a, b, c):
    # Convert degrees to radians for numpy
    return a * np.cos(np.deg2rad(b * x)) + c

# Check command line arguments
if len(sys.argv) < 2:
    print("Usage: python flux_analysis.py <csv_file_or_folder>")
    sys.exit(1)

path = sys.argv[1]

# Determine if path is a file or directory
files_to_process = []
if os.path.isdir(path):
    # Find all CSV files in the directory
    for filename in os.listdir(path):
        if filename.endswith('.csv'):
            files_to_process.append(os.path.join(path, filename))
    if not files_to_process:
        print(f"No CSV files found in directory: {path}")
        sys.exit(1)
else:
    # Single file mode
    files_to_process.append(path)

# Create plots
theta_fig = plt.figure(figsize=(12, 8))  # Main theta analysis figure
heatmap_grid = plt.figure(figsize=(15, 10))  # Heatmap figure with subplots

# Process each file
all_data = []
colors = plt.cm.tab10.colors  # Use tab10 colormap for distinct colors
marker_styles = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

for i, filepath in enumerate(files_to_process):
    # Get filename for legend
    filename = os.path.basename(filepath)
    color = colors[i % len(colors)]
    marker = marker_styles[i % len(marker_styles)]
    
    # Process the file
    data, metadata = process_file(filepath)
    if data is None:
        continue
        
    all_data.append((data, metadata, filename, color, marker))
    
    # Add subplot for heatmap (if we have multiple files)
    if len(files_to_process) > 1:
        ax = heatmap_grid.add_subplot(len(files_to_process) // 2 + len(files_to_process) % 2, 
                                      2 if len(files_to_process) > 1 else 1, 
                                      i + 1)
        # Create pivot table for 2D plotting
        pivot = data.pivot(index='theta', columns='phi', values='fraction')
        
        # Create heatmap
        im = ax.imshow(pivot, aspect='auto', origin='lower',
                       extent=[0, 360, 0, 90], interpolation='nearest', cmap='viridis')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Fraction of rays detected')
        
        ax.set_title(f"{filename}\n{metadata.get('BRDF Model', '')}")
        ax.set_xlabel('φ (degrees)')
        ax.set_ylabel('θ (degrees)')
        ax.grid(True)

# Now plot all theta analyses on the main figure
plt.figure(theta_fig.number)  # Switch to theta figure

for data, metadata, filename, color, marker in all_data:
    # Get unique theta values and calculate mean fraction for each theta
    grouped = data.groupby('theta')
    theta_means = grouped['fraction'].mean()
    theta_values = np.array(theta_means.index.tolist(), dtype=float)
    fraction_means = np.array(theta_means.values, dtype=float)
    
    # Calculate standard errors for error bars
    theta_std = grouped['fraction'].std().fillna(0.001)
    theta_counts = grouped.size()
    theta_errors = theta_std / np.sqrt(theta_counts)
    
    # Attempt curve fitting with error handling
    try:
        # Initial guess based on data
        p0 = [(np.max(fraction_means) - np.min(fraction_means))/2, 1.0, np.mean(fraction_means)]
        
        popt, pcov = curve_fit(cosine_func, theta_values, fraction_means, p0=p0)
        perr = np.sqrt(np.diag(pcov))
        fit_label = f'{filename}: {popt[0]:.3f}*cos({popt[1]:.3f}θ) + {popt[2]:.3f}'
        
    except Exception as e:
        print(f"Fit error for {filename}: {e}")
        # Use simple approximation
        popt = [np.mean(fraction_means)/2, 1.0, np.mean(fraction_means)/2]
        perr = [0, 0, 0]
        fit_label = f'{filename}: Approximate fit'
    
    # Generate smooth curve for plotting
    theta_smooth = np.linspace(min(theta_values), max(theta_values), 1000)
    fit_curve = cosine_func(theta_smooth, *popt)
    
    # Plot with error bars
    plt.errorbar(theta_values, fraction_means, yerr=theta_errors, fmt=marker, color=color, 
                 alpha=0.5, label=f'Data: {filename}', capsize=5, elinewidth=1, markeredgewidth=1)
    
    plt.plot(theta_smooth, fit_curve, '-', color=color, label=fit_label)
    
    # Calculate R-squared
    residuals = fraction_means - cosine_func(theta_values, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((fraction_means - np.mean(fraction_means))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Print fit parameters
    print(f"File: {filename}")
    print(f"  Fit parameters: a={popt[0]:.5f}, b={popt[1]:.5f}, c={popt[2]:.5f}")
    print(f"  R-squared value: {r_squared:.5f}")
    
    # Add exit port angle annotation if available (only for the first few files to avoid clutter)
    if i < 3:  # Limit to first 3 files to avoid cluttering
        if 'Exit port angle' in metadata:
            # Extract numeric part from string like "169 degrees"
            exit_angle_str = metadata['Exit port angle']
            exit_angle = float(exit_angle_str.split()[0])  # Take first part before any spaces
            calculated_angle = (180 - exit_angle) * 2
            # Position text in upper right area of plot
            plt.annotate(f"{filename}\nExit port angle: {exit_angle}°\nCalculated angle: {calculated_angle:.1f}°", 
                        xy=(0.95, 0.9 - i*0.1),  # Stack vertically
                        xycoords='axes fraction',
                        horizontalalignment='right',
                        verticalalignment='top',
                        color=color,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.8))
        elif 'Port angle' in metadata:  # Alternative key that might be used
            # Extract numeric part from string like "169 degrees"
            exit_angle_str = metadata['Port angle']
            exit_angle = float(exit_angle_str.split()[0])  # Take first part before any spaces
            calculated_angle = (180 - exit_angle) * 2
            plt.annotate(f"{filename}\nExit port angle: {exit_angle}°\nCalculated angle: {calculated_angle:.1f}°", 
                        xy=(0.95, 0.9 - i*0.1),  # Stack vertically
                        xycoords='axes fraction',
                        horizontalalignment='right',
                        verticalalignment='top',
                        color=color,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.8))

# Set up the main theta plot
plt.figure(theta_fig.number)
plt.xlabel('θ (degrees)')
plt.ylabel('Fraction')
plt.title('Flux Fraction vs Theta with Cosine Fit - Multiple Files Comparison')
plt.legend(loc='best', fontsize='small')
plt.grid(True)

# Adjust layout for both figures
theta_fig.tight_layout()
heatmap_grid.tight_layout()

# Save plots
if os.path.isdir(path):
    base_filename = os.path.basename(os.path.normpath(path))
else:
    base_filename = os.path.splitext(os.path.basename(path))[0]

# theta_fig.savefig(f"{base_filename}_theta_comparison.png", dpi=300, bbox_inches='tight')
# heatmap_grid.savefig(f"{base_filename}_heatmap_comparison.png", dpi=300, bbox_inches='tight')
# print(f"Plots saved as {base_filename}_theta_comparison.png and {base_filename}_heatmap_comparison.png")

# Show all plots
plt.show()

