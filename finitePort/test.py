import numpy as np
import matplotlib.pyplot as plt

# Define parameters
rho = 0.95  # Sphere reflectance
f = 0.3     # Fraction of surface occupied by the port
theta = np.linspace(0, np.pi/2, 100)  # Angle from 0 to 90 degrees
Phi_input = 1  # Arbitrary input flux normalization

# Compute effective flux inside the sphere
Phi_eff = Phi_input / (1 - rho * (1 - f))

# Compute observed flux distribution
Phi_theta = Phi_eff * f * np.cos(theta)

# Plot the distribution
plt.figure(figsize=(8,5))
plt.plot(np.degrees(theta), Phi_theta, label=r'$\Phi(\theta) \propto \cos\theta$')
plt.xlabel('Observation Angle θ (degrees)')
plt.ylabel('Relative Flux $\Phi(\theta)$')
plt.title('Light Flux Distribution as a Function of Angle')
plt.legend()
plt.grid()
plt.show()

