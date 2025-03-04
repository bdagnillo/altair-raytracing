import numpy as np
import scipy.integrate as spi

# Define constants
I0 = 1  # Example value for I0, adjust as needed
R = 1   # Example value for R, adjust as needed
a = 1   # Upper limit for r
epsilon = 1e-8  # Small value to prevent division by zero

# Define the integrand function with stability considerations
def integrand(r, phi, theta):
    sin_phi_tan_theta = np.sin(phi) * np.tan(theta)
    
    # Compute the argument of the square root safely
    radical = R**2 + r**2 - 2 * R * r * sin_phi_tan_theta
    radical = max(radical, epsilon)  # Ensure non-negative values
    
    numerator = (R - r * sin_phi_tan_theta)
    denominator = np.sqrt(radical)
    
    return (numerator / denominator) * r

# Compute the double integral
def compute_integral(theta):
    if theta >= np.pi / 2:  # Prevent singularity at 90 degrees
        raise ValueError("Theta must be less than 90 degrees (π/2 radians) to avoid instability.")

    result, error = spi.dblquad(integrand, 0, 2 * np.pi, lambda phi: 0, lambda phi: a, args=(theta,))
    return I0 * result

# Example: Compute the integral for theta values in a safe range
theta_values = np.linspace(0, np.pi / 2 - 0.01, 10)  # Avoid exactly 90 degrees
results = [compute_integral(theta) for theta in theta_values]

for theta, res in zip(theta_values, results):
    print(f"Integral result for theta={np.degrees(theta):.2f} degrees: {res:.6f}")

