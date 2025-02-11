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

def RK4(x, k, h):
    k1 = h * (k * x)
    k2 = h * (k * (x + k1 / 2))
    k3 = h * (k * (x + k2 / 2))
    k4 = h * (k * (x + k3))

    return x + (1/6) * (k1 + 2*k2 + 2*k3 + k4)

# Parameters
l = 1               # Length of the pendulum (m)
t_max = 10          # Total simulation time (s)
theta_0 = 0.5 * np.pi  # Initial angle (rad)
omega_0 = 0 * np.pi    # Initial angular velocity (rad/s)

# Array of step sizes to test
h_arr = ( np.linspace(0.0001, 0.01, 100) )
err = np.zeros(len(h_arr))  # Error array to store average errors

# Loop over different step sizes
for j in range(len(h_arr)):
    h = h_arr[j]
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

    # Compute error between analytical and numerical solutions
    diff = np.abs(theta_analytical - theta_num)
    err[j] = np.mean(diff)


def powerLaw(x,a,b): # as the error increases linearly on a log-log plot we can expect a power function to fit nicely
    return a*x**b

params, params_covariance = optimize.curve_fit(powerLaw, h_arr, err, p0=[1,2])

# ----------------------------------------------------


from ultimate_plotting_v7 import *

def ultimate_plot():
    
    sample_format_dict_1 = {
        "label"      : "numeric",          
        "fmt"        : '+', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 3, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"fit: $\\Delta \\propto h^{{ {params[1]:.2f} }} $",          
        "fmt"        : '--', 
        "color"      : "grey",                               
        "markersize" : 2, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "title"       : r"Error on Euler's method",
        "x_ax_label"  : r"h [ ]",
        "y_ax_label"  : r"absolute deviation [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    general_format_dict["log_scaling_xy"] = [True, True, 10]


    data_set_1  = h_arr,    None,   err,   None          
    data_set_2  = h_arr,    None,   powerLaw(h_arr, params[0],params[1] ), None
    all_data                = data_set_1 + data_set_2                              
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]

    save_plot = True, "./plot_simple_err.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)

ultimate_plot()
ultimate_plot()

