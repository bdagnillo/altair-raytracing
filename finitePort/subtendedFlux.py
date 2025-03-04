import numpy as np
import matplotlib.pyplot as plt

# Constants
R = 1.0  # sphere radius
a_values = np.linspace(0.1, 0.9, 5)  # port radii as fractions of R
theta = np.linspace(0, np.pi/2, 100)  # angles from 0 to 90 degrees
Phi_input = 1
rho_values = [0.95, 0.99, 1.00]  # different reflectance values

# Create figure with three subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Light Flux Distribution as a Function of Angle for Different Reflectances')

# Plot for each reflectance value
for idx, rho in enumerate(rho_values):
    for a in a_values:
        alpha = np.arcsin(a / R)  # port half-angle
        f = (a / R) ** 2  # area ratio
        Phi_theta = (Phi_input / (1 - rho * (1 - f))) * (0.5 * np.sin(alpha) ** 2) * np.cos(theta)
        
        axes[idx].plot(np.degrees(theta), Phi_theta, label=f'a/R = {a:.1f}')
    
    axes[idx].set_xlabel('Observation Angle θ (degrees)')
    axes[idx].set_ylabel('Relative Flux $\Phi(θ)$')
    axes[idx].set_title(f'ρ = {rho:.2f}')
    axes[idx].legend()
    axes[idx].grid(True)

plt.tight_layout()
plt.savefig('subtendedFlux.pdf')
plt.show()

plt.clf()

fig1 = plt.figure(figsize=(6, 5))
ax = fig1.add_subplot(111)  # Creates a single subplot that spans the entire figure
ax.set_title('Light Flux Distribution as a Function of Angle for Different Radius Ratios')
# Add a text annotation for the subtitle
ax.text(0.5, 0.95, '$\Phi(θ)$ $\propto$ $1/2 * sin^2(α)cos(θ)$', 
        horizontalalignment='center', transform=ax.transAxes)
ax.set_xlabel('Observation Angle θ (degrees)')
ax.set_ylabel('Relative Flux $\Phi(θ)$')

for a in a_values:
        alpha = np.arcsin(a / R)  # port half-angle
        f = (a / R) ** 2  # area ratio
        
        Phi_theta = f* (0.5 * np.sin(alpha) ** 2) * np.cos(theta)
        Phi_theta_no_correction = f * np.cos(theta)
        ax.plot(np.degrees(theta), Phi_theta, label=f'a/R = {a:.1f}')
        ax.plot(np.degrees(theta), Phi_theta_no_correction, label=f'a/R = {a:.1f} (no correction)', linestyle='--')

plt.legend() 
plt.show()