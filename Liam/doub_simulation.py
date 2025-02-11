import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import time

# Numerical ODE solvers

def normalize_theta(theta):
    return np.mod(theta + np.pi, 2 * np.pi) - np.pi

def Euler(theta, omega, h, a):
    """Euler method"""
    theta = normalize_theta(theta)
    theta_next = theta + omega * h
    omega_next = omega + a * h
    return theta_next, omega_next

def RK4(theta, omega, h, l1, l2, m1, m2, g, a_func):
    """RK4 method for double pendulum"""
    theta = normalize_theta(theta)

    # Compute k1
    k1_theta = omega
    k1_omega = a_func(theta, omega)

    # Compute k2 using k1
    k2_theta = omega + 0.5 * h * k1_omega
    k2_omega = a_func(theta + 0.5 * h * k1_theta, omega + 0.5 * h * k1_omega)

    # Compute k3 using k2
    k3_theta = omega + 0.5 * h * k2_omega
    k3_omega = a_func(theta + 0.5 * h * k2_theta, omega + 0.5 * h * k2_omega)

    # Compute k4 using k3
    k4_theta = omega + h * k3_omega
    k4_omega = a_func(theta + h * k3_theta, omega + h * k3_omega)

    # Update values
    theta_next = normalize_theta(theta + (h / 6) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta))
    omega_next = omega + (h / 6) * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega)

    return theta_next, omega_next

def compute_energy(theta1, theta2, omega1, omega2, m1, m2, l1, l2, g):
    """Computes total energy."""
    # Kinetic energy
    T = 0.5 * m1 * l1**2 * omega1**2 + 0.5 * m2 * (
        l1**2 * omega1**2
        + l2**2 * omega2**2
        + 2 * l1 * l2 * omega1 * omega2 * np.cos(theta1 - theta2)
    )
    # Potential energy
    V = -(m1 + m2) * g * l1 * np.cos(theta1) - m2 * g * l2 * np.cos(theta2)

    return T + V

def second_order_ode(t, theta1_0, theta2_0, omega1_0, omega2_0, l1, l2, m1, m2, h, method="Euler"):
    """Numerically solves double pendulum using Euler or RK4."""
    g = const.g
    n = t.shape[0]
    batch_size = len(theta1_0)

    theta1 = np.zeros((batch_size, n))
    theta2 = np.zeros((batch_size, n))
    omega1 = np.zeros((batch_size, n))
    omega2 = np.zeros((batch_size, n))
    energy = np.zeros((batch_size, n))

    # Set initial conditions
    theta1[:, 0] = normalize_theta(theta1_0)
    theta2[:, 0] = normalize_theta(theta2_0)
    omega1[:, 0] = omega1_0
    omega2[:, 0] = omega2_0
    energy[:, 0] = compute_energy(theta1[:, 0], theta2[:, 0], omega1[:, 0], omega2[:, 0], m1, m2, l1, l2, g)

    # Poincare section arrays
    poincare_theta2 = [[] for _ in range(batch_size)]
    poincare_omega2 = [[] for _ in range(batch_size)]

    for i in range(1, n):
        delta_theta = theta1[:, i - 1] - theta2[:, i - 1]
        sin_delta = np.sin(delta_theta)
        cos_delta = np.cos(delta_theta)
        denominator = m1 + m2 * sin_delta**2

        def a1_func(theta, omega):
            return (-g * (2 * m1 + m2) * np.sin(theta)
                    - m2 * g * np.sin(theta - 2 * theta2[:, i - 1])
                    - 2 * sin_delta * m2 * (omega2[:, i - 1]**2 * l2 + omega**2 * l1 * cos_delta)) \
                   / (l1 * denominator)

        def a2_func(theta, omega):
            return (2 * sin_delta * (omega1[:, i - 1]**2 * l1 * (m1 + m2)
                    + g * (m1 + m2) * np.cos(theta1[:, i - 1])
                    + omega**2 * l2 * m2 * cos_delta)) \
                   / (l2 * denominator)

        if method == "Euler":
            theta1[:, i], omega1[:, i] = Euler(theta1[:, i - 1], omega1[:, i - 1], h, a1_func(theta1[:, i - 1], omega1[:, i - 1]))
            theta2[:, i], omega2[:, i] = Euler(theta2[:, i - 1], omega2[:, i - 1], h, a2_func(theta2[:, i - 1], omega2[:, i - 1]))
        elif method == "RK4":
            theta1[:, i], omega1[:, i] = RK4(theta1[:, i - 1], omega1[:, i - 1], h, l1, l2, m1, m2, g, a1_func)
            theta2[:, i], omega2[:, i] = RK4(theta2[:, i - 1], omega2[:, i - 1], h, l1, l2, m1, m2, g, a2_func)

        energy[:, i] = compute_energy(
            normalize_theta(theta1[:, i]),
            normalize_theta(theta2[:, i]),
            omega1[:, i],
            omega2[:, i],
            m1,
            m2,
            l1,
            l2,
            g
        )

        # Poincaré section
        mask = (theta1[:, i - 1] < 0) & (theta1[:, i] > 0) & (omega1[:, i] > 0)
        for j in range(batch_size):
            if mask[j]:
                poincare_theta2[j].append(theta2[j, i])
                poincare_omega2[j].append(omega2[j, i])

    return theta1, theta2, omega1, omega2, energy, poincare_theta2, poincare_omega2

# Parameters
t_max = 30
theta1_0 = 0
omega1_0 = 0
l1 = 1
l2 = 1
m1 = 1
m2 = 100
h = 0.0001

# Generates a range of theta2 and omega2 that give the same constant energy, to generate poincare sections
t = np.arange(0.0, t_max, h)
theta2_0_vals = 0
omega2_0_vals = 0.01 * np.pi

# Timer start
t_start = time.time()

# Solve ODE for all initial conditions using RK4
theta1, theta2, omega1, omega2, energy, poincare_theta2, poincare_omega2 = second_order_ode(
    t,
    np.full(1, theta1_0),
    theta2_0_vals,
    np.full(1, omega1_0),
    omega2_0_vals,
    l1,
    l2,
    m1,
    m2,
    h,
    method="RK4"
)

# Timer end
t_end = time.time()
print(f"Execution time: {t_end - t_start:.2f} seconds")

# Plot Theta1 and Theta2 vs Time
plt.figure(figsize=(10, 6))
plt.plot(t, theta1[0], label='Theta1 (rad)')
plt.plot(t, theta2[0], label='Theta2 (rad)')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.title('Theta1 and Theta2 vs Time')
plt.legend()
plt.grid()
plt.show()

# Plot Position
plt.figure(figsize=(10, 6))
plt.plot(theta2[0], theta1[0])
plt.xlabel('Theta1 (rad)')
plt.ylabel('Theta2 (rad)')
plt.title('Double Pendulum Position')
plt.grid()
plt.show()

# Plot Phasespace
plt.figure(figsize=(10, 6))
plt.plot(theta1[0], omega1[0], label='1')
plt.plot(theta2[0], omega2[0], label='2')
plt.xlabel('Theta1 (rad)')
plt.ylabel('Theta2 (rad)')
plt.title('Double Pendulum Position')
plt.legend()
plt.grid()
plt.show()

# Plot Energy vs Time
plt.figure(figsize=(10, 6))
plt.plot(t, energy[0], label='Total Energy')
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.title('Energy vs Time')
plt.legend()
plt.grid()
plt.show()