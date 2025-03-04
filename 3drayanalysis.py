import numpy as np
import matplotlib.pyplot as plt

# Read the ray data from the file
data = np.loadtxt('/home/bdagn/altair-raytracing/3dRayLog.txt')  # Removed delimiter=',' 

# Extract x and z components
x_components = data[:, 0]  # First column contains x components
z_components = data[:, 2]  # Third column contains z components

# Filter rays near x=0 (±1)
mask = np.abs(x_components) <= 1.0
filtered_z = z_components[mask]

# Calculate angles for filtered rays
z_angle_dist = np.arccos(filtered_z) * 180 / np.pi - 180

# Create histogram
plt.figure(figsize=(10, 6))
plt.hist(z_angle_dist, bins=100, edgecolor='black')
plt.xlabel('Z Angle (degrees)')
plt.ylabel('Frequency')
plt.title('Distribution of Ray Z Angles (at x = 0 ± 1)')
plt.grid(True, alpha=0.3)

# Save the plot
# plt.savefig('/home/bdagn/altair-raytracing/z_distribution.png')
plt.show()
plt.close()
