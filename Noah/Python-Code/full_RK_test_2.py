import numpy as np
import matplotlib.pyplot as plt

# Constants
g_grav      = 9.81
length      = 1.0
mass        = 1.0
t_max       = 1.0
theta0      = 0.1 * np.pi  
theta_dot0  = 0.0

# Right-hand side of the differential equation for the simple pendulum
def y_dot_simp_pend(t, y):
    dy = np.zeros(2)
    dy[0] = y[1]
    dy[1] = -g_grav/length * y[0]  # Small-angle approximation (linearized)
    return dy

# Second-order Runge-Kutta method (Midpoint Method)
def RK_2(ODE, y0, h):
    steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, steps)
    y = np.zeros((steps, len(y0)))
    y[0] = y0

    for i in range(steps - 1):
        k1 = h * ODE(t[i], y[i])
        k2 = h * ODE(t[i] + h/2, y[i] + k1/2)
        y[i+1] = y[i] + (k1 + k2) / 2

    return y

# Fourth-order Runge-Kutta method 
def RK_4(ODE, y0, h):
    steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, steps)
    y = np.zeros((steps, len(y0)))
    y[0] = y0

    for i in range(steps - 1):
        k1      = h * ODE(t[i],         y[i])
        k2      = h * ODE(t[i] + h/2,   y[i] + k1/2)
        k3      = h * ODE(t[i] + h/2,   y[i] + k2/2)
        k4      = h * ODE(t[i] + h/1,   y[i] + k3)
        y[i+1] = y[i] + k1/6 + k2/3 + k3/3 + k4/6

    return y


# Analytical solution of the simple pendulum (small angle approximation)
def analyt_pend_sol(y0, h):
    steps = int(t_max / h) + 1  # Ensure consistency
    t = np.linspace(0, t_max, steps)
    omega = np.sqrt(g_grav / length)
    theta_t = y0[0] * np.cos(omega * t) + y0[1] / omega * np.sin(omega * t)
    return theta_t

# Compute deviation between numerical and analytical solutions for various step sizes
def calc_deviation(h_array, procedure):
    deviation = np.zeros(len(h_array))

    for i, h in enumerate(h_array):
        steps = int(t_max / h) + 1
        num_theta_t = procedure(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 0]
        ana_theta_t = analyt_pend_sol([theta0, theta_dot0], h)

        # Root Mean Square (RMS) Error
        deviation[i] = np.sqrt(np.sum((num_theta_t - ana_theta_t)**2) / len(ana_theta_t))

    return deviation

# Define step sizes
h_array = np.logspace(-5, -1, 20)  # Logarithmic spacing for better visualization

# Compute deviation
deviation_RK2 = calc_deviation(h_array, RK_2)
deviation_RK4 = calc_deviation(h_array, RK_4)


# Plot deviation vs step size on a log-log scale
plt.figure(figsize=(8,6))
plt.loglog(h_array, deviation_RK2, 'o-', label='RK2 Error')
plt.loglog(h_array, deviation_RK4, 'o-', label='RK4 Error')
plt.loglog(h_array, h_array**2, '--', label='$h^2$ Reference')  # Expected slope
plt.loglog(h_array, h_array**4, '--', label='$h^4$ Reference')  # Expected slope
plt.xlabel("Step Size (h)")
plt.ylabel("Global Error")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.title("RK2 Global Error Scaling with h")
plt.show()
