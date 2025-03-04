import numpy as np
import matplotlib.pyplot as plt

# Define parameters
R = 1.0  # Radius of the integrating sphere (arbitrary units)
r_p = 0.1  # Radius of the exit port

R = float(input("Radius of integrating sphere:"))
r_p = float(input("Radius of exit port:"))

# Define theta range (observer angle)
theta_deg = np.linspace(0, 90, 100)  # From 0 to 90 degrees
theta_rad = np.radians(theta_deg)  # Convert to radians

# Improved version with better handling of numerical stability

def safe_projection_factor(theta, R, r_p, num_points=100):
    """
    Compute the flux projection factor with improved numerical stability.
    """
    r_vals = np.linspace(0, r_p, num_points)  # Radial points on the port
    phi_vals = np.linspace(0, 2 * np.pi, num_points)  # Angular points
    R_grid, Phi_grid = np.meshgrid(r_vals, phi_vals)  # Create grid

    # Compute the denominator of cos(theta') safely
    denominator = np.sqrt(np.maximum(R**2 + R_grid**2 - 2 * R * R_grid * np.sin(Phi_grid) * np.tan(theta), 1e-10))

    # Compute cos(theta') for each point on the port
    cos_theta_prime = (R - R_grid * np.sin(Phi_grid) * np.tan(theta)) / denominator
    
    print(np.max(cos_theta_prime))

    # Ensure values are in valid range [-1, 1]
    cos_theta_prime = np.clip(cos_theta_prime, -1, 1)

    # Integrate over the port surface
    dA = R_grid * (r_p / num_points) * (2 * np.pi / num_points)  # Differential area elements
    flux = np.sum(cos_theta_prime * dA)  # Sum contributions

    return flux

# Recompute the flux distribution with the improved function
flux_distribution_improved = np.array([safe_projection_factor(t, R, r_p) for t in theta_rad])

# Normalize flux for better visualization
flux_distribution_improved /= np.max(flux_distribution_improved)

# Plot the improved results
plt.figure(figsize=(8, 6))
plt.plot(theta_deg, flux_distribution_improved, label="Finite Port Projection Factor (Improved)", color="blue", linewidth=2)
plt.plot(theta_deg, np.cos(theta_rad), label="Ideal Lambertian (cosθ)", linestyle="dashed", color="red")
plt.xlabel("Observer Angle θ (degrees)")
plt.ylabel("Normalized Observed Flux")
plt.title("Flux Distribution from an Integrating Sphere Exit Port (Improved)")
plt.legend()
plt.grid(True)
plt.show()

