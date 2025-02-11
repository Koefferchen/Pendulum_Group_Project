import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy import optimize

def normalize_theta(theta):
    return np.mod(theta + np.pi, 2 * np.pi) - np.pi

# Analytical solution to the simple pendulum (small angle approximation)
def simple_pendulum_analytical(theta_0, l, t):
    return theta_0 * np.cos(np.sqrt(const.g / l) * t)

# Euler forward method for solving simple pendulum
def Euler(theta, omega, h, l):
    theta = normalize_theta(theta)
    a = -const.g / l * theta  # Angular acceleration
    theta_next = theta + omega * h
    omega_next = omega + a * h
    return theta_next, omega_next

# Parameters
l = 1               # Length of the pendulum (m)
t_max = 100          # Total simulation time (s)
theta_0 = 0.5 * np.pi  # Initial angle (rad)
omega_0 = 0 * np.pi    # Initial angular velocity (rad/s)
h = 0.0001

t = np.arange(0, t_max, h)  # Adjust time array according to current step size

theta_analytical = simple_pendulum_analytical(theta_0, l, t)
    
# Initialize numerical solution arrays
theta_num = np.zeros(len(t))
omega_num = np.zeros(len(t))
theta_num[0] = theta_0
omega_num[0] = omega_0

# Numerical solution using Euler's method
for i in range(1, len(t)):
    theta_num[i], omega_num[i] = Euler(theta_num[i-1], omega_num[i-1], h, l)

plt.plot(t, theta_num, label='Euler solution')
plt.plot(t, theta_analytical, label='Analytical solution')
plt.legend()
plt.show()


