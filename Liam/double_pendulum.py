import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import time

def normalize_theta(theta):
    return np.mod(theta + np.pi, 2 * np.pi) - np.pi

def Euler(theta, omega, h, a):
    """Euler method"""
    theta = normalize_theta(theta)
    theta_next = theta + omega * h
    omega_next = omega + a * h
    return theta_next, omega_next

def RK4(theta, omega, h, a):
    """RK4 method"""
    theta = normalize_theta(theta)

    k1_theta = omega
    k1_omega = a

    k2_theta = omega + 0.5 * h * k1_omega
    k2_omega = a  

    k3_theta = omega + 0.5 * h * k2_omega
    k3_omega = a  

    k4_theta = omega + h * k3_omega
    k4_omega = a  

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

<<<<<<< HEAD
    return T + V
=======
    return T + 2 * V
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e

def second_order_ode(t, theta1_0, theta2_0, omega1_0, omega2_0, l1, l2, m1, m2, h, method="Euler"):
    """Numerically solves double pendulum."""
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

    # Poincaré section
<<<<<<< HEAD
    poincare_theta2 = []
    poincare_omega2 = []
=======
    poincare_theta2 = [[] for _ in range(batch_size)]
    poincare_omega2 = [[] for _ in range(batch_size)]
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e

    for i in range(1, n):
        delta_theta = theta1[:, i - 1] - theta2[:, i - 1]
        sin_delta = np.sin(delta_theta)
        cos_delta = np.cos(delta_theta)
        denominator = m1 + m2 * sin_delta**2

        a1 = (-g * (2 * m1 + m2) * np.sin(theta1[:, i - 1])
              - m2 * g * np.sin(theta1[:, i - 1] - 2 * theta2[:, i - 1])
              - 2 * sin_delta * m2 * (omega2[:, i - 1]**2 * l2 + omega1[:, i - 1]**2 * l1 * cos_delta)) \
             / (l1 * denominator)

        a2 = (2 * sin_delta * (omega1[:, i - 1]**2 * l1 * (m1 + m2)
              + g * (m1 + m2) * np.cos(theta1[:, i - 1])
              + omega2[:, i - 1]**2 * l2 * m2 * cos_delta)) \
             / (l2 * denominator)

        if method == "Euler":
            theta1[:, i], omega1[:, i] = Euler(theta1[:, i - 1], omega1[:, i - 1], h, a1)
            theta2[:, i], omega2[:, i] = Euler(theta2[:, i - 1], omega2[:, i - 1], h, a2)
        elif method == "RK4":
            theta1[:, i], omega1[:, i] = RK4(theta1[:, i - 1], omega1[:, i - 1], h, a1)
            theta2[:, i], omega2[:, i] = RK4(theta2[:, i - 1], omega2[:, i - 1], h, a2)

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

        # Poincaré section: theta1 crosses zero from negative to positive and omega1 positive
        mask = (theta1[:, i - 1] < 0) & (theta1[:, i] > 0) & (omega1[:, i] > 0)
<<<<<<< HEAD
        poincare_theta2.extend(theta2[mask, i])
        poincare_omega2.extend(omega2[mask, i])
=======
        for j in range(batch_size):
            if mask[j]:
                poincare_theta2[j].append(theta2[j, i])
                poincare_omega2[j].append(omega2[j, i])

>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e

    return theta1, theta2, omega1, omega2, energy, poincare_theta2, poincare_omega2

# Parameters
<<<<<<< HEAD
t_max = 5
=======
t_max = 1000
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e
theta1_0 = 0
omega1_0 = 0
l1 = 1
l2 = 1
m1 = 1
m2 = 1
<<<<<<< HEAD
h = 0.0001
constant_energy = 10
=======
h = 0.001
constant_energy = 2000
>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e

# Prompt user for Poincare section generation
generate_poincare = input("Do you want to generate a Poincare map? (yes/no): ").strip().lower() == "yes"
num_points = int(input("Enter the number of initial points for the Poincare map: ")) if generate_poincare else 1

# Gets initial theta2 and omega2 values for Poincare sections 
t = np.arange(0.0, t_max, h)

theta2_0_vals = np.linspace(0, 2 * np.pi, num_points)

def compute_initial_omega2(constant_energy, theta2_0):
    return np.sqrt(
        np.maximum(
            0,
            (2 * constant_energy
             + (m1 + m2) * const.g * l1 * np.cos(theta1_0)
             + m2 * const.g * l2 * np.cos(theta2_0)
             - 0.5 * (m1 + m2) * l1**2 * omega1_0**2)
            / (0.5 * m2 * l2**2)
        )
    )

omega2_0_vals = compute_initial_omega2(constant_energy, theta2_0_vals)

# Timer start
t_start = time.time()

# Solve ODE for all initial conditions using RK4
theta1, theta2, omega1, omega2, energy, poincare_theta2, poincare_omega2 = second_order_ode(
    t,
    np.full(num_points, theta1_0),
    theta2_0_vals,
    np.full(num_points, omega1_0),
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

# Plot Energy vs Time
plt.figure(figsize=(10, 6))
plt.plot(t, energy[0], label='Total Energy')
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.title('Energy vs Time')
plt.legend()
plt.grid()
plt.show()

# Plot Poincare section if requested
if generate_poincare:
    plt.figure(figsize=(8, 6))
<<<<<<< HEAD
    plt.scatter(poincare_theta2, poincare_omega2, s=1, color='black', alpha=0.2)
    plt.xlabel('Theta2 (rad)')
    plt.ylabel('Omega2 (rad/s)')
    plt.title('Poincaré Section')
    plt.grid()
    plt.show()
=======
    colors = plt.cm.viridis(np.linspace(0, 1, num_points))  # Generate colors
    for j in range(num_points):
        plt.scatter(poincare_theta2[j], poincare_omega2[j], s=5, alpha=0.6, color=colors[j])
    plt.xlabel('Theta2 (rad)')
    plt.ylabel('Omega2 (rad/s)')
    plt.xlim(-np.pi, np.pi)
    plt.title('Poincaré Section')
    plt.legend()
    plt.grid()
    plt.show()

>>>>>>> 3a169a0a2e11b347a0465350e6c38331db4a045e
