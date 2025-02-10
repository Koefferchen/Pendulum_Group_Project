
from ultimate_plotting_v7 import *

# Simple Pendulum Constants
g_grav      = 9.81
length      = 1.0
mass        = 1.0
t_max       = 2.0
theta0      = 0.1 * np.pi  
theta_dot0  = 0.0

# Exponetial Constants
alpha       = 2.0
x_0         = 1.0

# Gauss function Constants
beta        = -0.1
z_0         = 1.0


# Right-hand side of the differential equation for the simple pendulum
def ODE_simp_pend(t, y):
    dy = np.zeros(2)
    dy[0] = y[1]
    dy[1] = -g_grav/length * y[0]  # Small-angle approximation (linearized)
    return dy


# Analytical solution of the simple pendulum (small angle approximation)
def analyt_pend_sol(y0, h):
    steps = int(t_max / h) + 1  # Ensure consistency
    t = np.linspace(0, t_max, steps)
    omega = np.sqrt(g_grav / length)
    theta_t = y0[0] * np.cos(omega * t) + y0[1] / omega * np.sin(omega * t)
    return theta_t


# Right-hand side of the differential equation for the exponential function
def ODE_exp(t, y):
    dy = np.zeros(2, dtype=np.float128)
    dy[0] = alpha * y[0]
    dy[1] = 0.0
    return dy


# Analytical solution of the exp ODE
def analyt_exp_sol(y0, h):
    steps = int(t_max / h) + 1  # Ensure consistency
    t = np.linspace(0, t_max, steps)
    sol_t = x_0 * np.exp( alpha * t )
    return sol_t


# Right-hand side of the differential equation for the gauss function
def ODE_gauss(t, y):
    dy = np.zeros(2, dtype=np.float128)
    dy[0] = beta * 2 * t * y[0]
    dy[1] = 0.0
    return dy


# Analytical solution of the Gauss ODE
def analyt_gauss_sol(y0, h):
    steps = int(t_max / h) + 1  # Ensure consistency
    t = np.linspace(0, t_max, steps)
    sol_t = z_0 * np.exp( beta * t**2 )
    return sol_t


# Second-order Runge-Kutta method (Midpoint Method)
def RK_2(ODE, y0, h):
    steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, steps, dtype=np.float128)
    y = np.zeros((steps, len(y0)), dtype=np.float128)
    y[0] = y0

    for i in range(steps - 1):
        k1 = h * ODE(t[i], y[i])
        k2 = h * ODE(t[i] + h/2, y[i] + k1/2)
        y[i+1] = y[i] + (k1 + k2) / 2

    return y


# Fourth-order Runge-Kutta method 
def RK_4(ODE, y0, h):
    steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, steps, dtype=np.float128)
    y = np.zeros((steps, len(y0)), dtype=np.float128)
    y[0] = y0

    for i in range(steps - 1):
        k1      = h * ODE(t[i],         y[i])
        k2      = h * ODE(t[i] + h/2,   y[i] + k1/2)
        k3      = h * ODE(t[i] + h/2,   y[i] + k2/2)
        k4      = h * ODE(t[i] + h/1,   y[i] + k3)
        y[i+1] = y[i] + k1/6 + k2/3 + k3/3 + k4/6

    return y


# Sixth-order Runge-Kutta method 
def RK_6(ODE, y0, h):
    steps = int(t_max / h) + 1  
    t = np.linspace(0, t_max, steps, dtype=np.float128)
    y = np.zeros((steps, len(y0)), dtype=np.float128)
    y[0] = y0

    for i in range(steps - 1):
        k1      = h * ODE(t[i],         y[i])
        k2      = h * ODE(t[i] + h*1/3, y[i] + k1*1/3 )
        k3      = h * ODE(t[i] + h*2/3, y[i] + k1*0 + k2*2/3 )
        k4      = h * ODE(t[i] + h*1/3, y[i] + k1*(1/12) + k2*(1/3) + k3*(-1/12) )
        k5      = h * ODE(t[i] + h*1/2, y[i] + k1*(-1/16) + k2*(9/8) + k3*(-3/16) + k4*(-3/8) )
        k6      = h * ODE(t[i] + h*1/2, y[i] + k1*0 + k2*(9/8) + k3*(-3/8) + k4*(-3/4) + k5*1/2 )
        k7      = h * ODE(t[i] + h*1/1, y[i] + k1*9/44 + k2*(-9/11) + k3*63/44 + k4*18/11 + k5*0 + k6*(-16/11) ) 
        y[i+1] = y[i] + k1*11/120 + k2*0 + k3*27/40 + k4*27/40 + k5*(-4/15) + k6*(-4/15) + k7*11/120

    return y


# Compute deviation between numerical and analytical solutions for various step sizes
def calc_deviation_pend(h_array, procedure):
    deviation = np.zeros(len(h_array), dtype=np.float128)

    for i, h in enumerate(h_array):
        steps = int(t_max / h) + 1
        num_theta_t = procedure(ODE_simp_pend, [theta0, theta_dot0], h)[:, 0]
        ana_theta_t = analyt_pend_sol([theta0, theta_dot0], h)

        # Root Mean Square (RMS) Error
        deviation[i] = np.sqrt(np.sum((num_theta_t - ana_theta_t)**2) / len(ana_theta_t))

    return deviation

# Compute deviation between numerical and analytical solutions for various step sizes
def calc_deviation_exp(h_array, procedure):
    deviation = np.zeros(len(h_array), dtype=np.float128)

    for i, h in enumerate(h_array):
        steps = int(t_max / h) + 1
        num_sol = procedure(ODE_exp, [x_0, 0], h)[:, 0]
        ana_sol = analyt_exp_sol( [x_0, 0], h )

        print(num_sol)
        # Root Mean Square (RMS) Error
        deviation[i] = np.sqrt(np.sum((num_sol - ana_sol)**2) / len(ana_sol))

    return deviation

# Compute deviation between numerical and analytical solutions for various step sizes
def calc_deviation_gauss(h_array, procedure):
    deviation = np.zeros(len(h_array), dtype=np.float128)

    for i, h in enumerate(h_array):
        steps = int(t_max / h) + 1
        num_sol = procedure(ODE_gauss, [z_0, 0], h)[:, 0]
        ana_sol = analyt_gauss_sol( [z_0, 0], h )

        print(num_sol)
        # Root Mean Square (RMS) Error
        deviation[i] = np.sqrt(np.sum((num_sol - ana_sol)**2) / len(ana_sol))

    return deviation


# Define step sizes
h_array = np.logspace(-4, -2, 20, dtype=np.float128)  # Logarithmic spacing for better visualization

# Compute deviation
deviation_RK2 = calc_deviation_gauss(h_array, RK_2)
deviation_RK4 = calc_deviation_gauss(h_array, RK_4)
deviation_RK6 = calc_deviation_gauss(h_array, RK_6)

x_fit_RK2, y_fit_RK2, a_RK2, b_RK2 =  linear_fit( np.log(h_array), np.log(deviation_RK2), False, None, True )
#x_fit_RK4, y_fit_RK4, a_RK4, b_RK4 =  linear_fit( np.log(h_array), np.log(deviation_RK4) )
#x_fit_RK6, y_fit_RK6, a_RK6, b_RK6 =  linear_fit( np.log(h_array), np.log(deviation_RK6) )

# Plot deviation vs step size on a log-log scale    
plt.figure(figsize=(8,6))
plt.loglog(h_array, deviation_RK2, 'o-', label='RK2 Error')
plt.loglog(h_array, deviation_RK4, 'o-', label='RK4 Error')
plt.loglog(h_array, deviation_RK6, 'o-', label='RK6 Error')
plt.loglog(h_array, h_array**2, '--', label='$h^2$ Reference')  # Expected slope
plt.loglog(h_array, h_array**4, '--', label='$h^4$ Reference')  # Expected slope
plt.loglog(h_array, h_array**6, '--', label='$h^6$ Reference')  # Expected slope
plt.xlabel("Step Size (h)")
plt.ylabel("Global Error")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.title("RK2 Global Error Scaling with h")
plt.savefig("../plots/plot_DEV.jpg")