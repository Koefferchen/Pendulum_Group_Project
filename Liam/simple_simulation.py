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

def RK4():
    return theta_next, omega_next

# Parameters
g = const.g
l = 1                     # Length of the pendulum (m)
t_max = 10                  # Total simulation time (s)
theta_0 = 0.05 * np.pi       # Initial angle (rad)
omega_0_stand = np.sqrt(g / l) * np.sqrt(2 * (1+np.cos(theta_0))) # Calculates velocity needed to make pendulum stand up
omega_0 = 0     # Initial angular velocity (rad/s)
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

# Plotting
plt.plot(t, theta_num, label='Numerical Solution')
plt.plot(t, theta_analytical, label='Analytical solution')
plt.xlabel('Time (s)')
plt.ylabel('Theta (rad)')
plt.title('Solutions to Simple Pendulum')
plt.legend()
plt.show()

plt.title('Phasespace')
plt.plot(theta_num, omega_num)
plt.xlabel('Theta (rad)')
plt.ylabel('Omega (rad/s)')
plt.show()


