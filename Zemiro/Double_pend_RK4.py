import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
from scipy import constants


def RK4_double_pendulum(t1, theta1_0, theta2_0, theta1_dot_0, theta2_dot_0, L1, L2, m1, m2, g, h):
    """
    ODE Solver for the double pendulum using Runge-Kutta 4th order Method.

    Args:
        t1: Time to solve the equation up to
        theta1_0: Initial angle of the first pendulum
        theta2_0: Initial angle of the second pendulum
        theta1_dot_0: Initial angular velocity of the first pendulum
        theta2_dot_0: Initial angular velocity of the second pendulum
        L1, L2: Lengths of the pendulum arms
        m1, m2: Masses of the pendulum bobs
        g: Gravitational acceleration
        h: Step size of the method

    Returns:
        t: Array of time points
        theta1: Array of angle values for the first pendulum
        theta2: Array of angle values for the second pendulum
        theta1_dot: Array of angular velocity values for the first pendulum
        theta2_dot: Array of angular velocity values for the second pendulum
    """
    def derivatives(y):
        """
        Compute the derivatives for the double pendulum.
        Args:
            y: State vector [theta1, theta2, theta1_dot, theta2_dot]

        Returns:
            dydt: Derivatives [theta1_dot, theta2_dot, theta1_ddot, theta2_ddot]
        """
        theta1, theta2, theta1_dot, theta2_dot = y

        mass_diff=m2/(m1+m2)
        delta = theta1 - theta2
        denom1 = (1-((np.cos(delta)*np.cos(delta)) * mass_diff))
        denom2 = (1-((np.cos(delta)*np.cos(delta)) * mass_diff))

        theta1_ddot = (
            g/L1 *( (mass_diff * np.sin(theta2) * np.cos(delta) ) - np.sin(theta1) )
            - (np.sin(delta) * mass_diff * (theta2_dot*theta2_dot * (L2/L1) + theta1_dot*theta1_dot * np.cos(delta)))
        ) / denom1

        theta2_ddot = (
                g/L2 *( (mass_diff * np.sin(theta1) * np.cos(delta)) - np.sin(theta2) )
            + np.sin(delta) * (theta1_dot*theta1_dot * (L1/L2) + theta2_dot*theta2_dot * mass_diff * np.cos(delta))
       
            ) / denom2

        return np.array([theta1_dot, theta2_dot, theta1_ddot, theta2_ddot])

    # Initialize arrays
    t = np.arange(0,t1+h,h)
    n = np.shape(t)[0]
    theta1 = np.zeros(n)
    theta2 = np.zeros(n)
    theta1_dot = np.zeros(n)
    theta2_dot = np.zeros(n)

    # Set initial conditions
    theta1[0] = theta1_0
    theta2[0] = theta2_0
    theta1_dot[0] = theta1_dot_0
    theta2_dot[0] = theta2_dot_0

    # Main RK4 loop
    for i in range(1, n):
        y = np.array([theta1[i - 1], theta2[i - 1], theta1_dot[i - 1], theta2_dot[i - 1]])

        k1 = h * derivatives(y)
        k2 = h * derivatives(y + 0.5 * k1)
        k3 = h * derivatives(y + 0.5 * k2)
        k4 = h * derivatives(y + k3)

        y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        theta1[i], theta2[i], theta1_dot[i], theta2_dot[i] = y_next

    theta1 = np.mod(theta1+np.pi, 2*np.pi)-np.pi
    theta2 = np.mod(theta2+np.pi, 2*np.pi)-np.pi
    
    return t, theta1, theta2, theta1_dot, theta2_dot

# Parameters
L1, L2 = 1.0, 1.0  # Lengths of the pendulums
m1, m2 = 1.0, 1.0  # Masses of the pendulum bobs
g =sp.constants.g         # Gravitational acceleration
theta1_0, theta2_0 = 1*np.pi, 0.01*np.pi  # Initial angles
theta1_dot_0, theta2_dot_0 = 0.0, 0.0         # Initial angular velocities
t1 = 10.0          # Simulation time
h = 0.001           # Time step

# Simulate
t, theta1, theta2, theta1_dot, theta2_dot = RK4_double_pendulum(
    t1, theta1_0, theta2_0, theta1_dot_0, theta2_dot_0, L1, L2, m1, m2, g, h
)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(t, theta1, label=r'$\theta_1$')
plt.plot(t, theta2, label=r'$\theta_2$')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.legend()
plt.title('Double Pendulum Simulation with RK4')
plt.grid()
plt.show()
