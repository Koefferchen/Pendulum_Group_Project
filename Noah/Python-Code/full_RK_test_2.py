
from ultimate_plotting_v7 import *

    # Global Constants
g_grav      = 9.81
length      = 1.0
mass        = 1.0
t_max       = 1.0
theta0      = 0.3 * np.pi  
theta_dot0  = 0.0

    # Right-hand side of the differential equation for the simple pendulum
def y_dot_simp_pend(t, y):
    dy = np.zeros(2, dtype=np.float128)
    dy[0] = y[1]
    dy[1] = -g_grav/length * y[0]  # Small-angle approximation (linearized)
    return dy

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
        k1      = h * ODE(t[i] + h*0/1,   y[i])
        k2      = h * ODE(t[i] + h*1/3,   y[i] + k1*(1/3))
        k3      = h * ODE(t[i] + h*2/3,   y[i] + k1*(0/1) + k2*(2/3) )
        k4      = h * ODE(t[i] + h*1/3,   y[i] + k1*(1/12) + k2*(1/3) + k3*(-1/12)  )
        k5      = h * ODE(t[i] + h*1/2,   y[i] + k1*(-1/16) + k2*(9/8) + k3*(-3/16) + k4*(-3/8) )
        k6      = h * ODE(t[i] + h*1/2,   y[i] + k1*(0/1) + k2*(9/8) + k3*(-3/8) + k4*(-3/4) + k5*(1/2) )
        k7      = h * ODE(t[i] + h*1/1,   y[i] + k1*(9/44) + k2*(-9/11) + k3*(63/44) + k4*(18/11) + k5*(0/1) + k6*(-16/11) )
        y[i+1] = y[i] + k1*(11/120) + k2*(0/1) + k3*(27/40) + k4*(27/40) + k5*(-4/15) + k6*(-4/15) + k7*(11/120)
    return y

    # Analytical solution of the simple pendulum (small angle approximation)
def analyt_pend_sol(y0, h):
    steps       = int(t_max / h) + 1  # Ensure consistency
    t           = np.linspace(0, t_max, steps, dtype=np.float128)
    omega       = np.sqrt(g_grav / length)
    theta_t     = y0[0] * np.cos(omega * t) + y0[1] / omega * np.sin(omega * t)
    theta_dot_t = -omega* y0[0] * np.sin(omega * t) + y0[1] * np.cos(omega * t)
    return theta_t, theta_dot_t

def calc_energy(theta, theta_dot):
    return 0.5 * mass * length**2 * theta_dot**2 + g_grav * mass * length *(1 - np.cos(theta))

    # Compute deviation between numerical and analytical solutions for various step sizes
def calc_deviation(h_array, procedure):
    deviation_theta     = np.zeros(len(h_array), dtype=np.float128)
    deviation_theta_dot = np.zeros(len(h_array), dtype=np.float128)
    deviation_energy    = np.zeros(len(h_array), dtype=np.float128)

    for i, h in enumerate(h_array):
        steps = int(t_max / h) + 1
        
        num_theta_t     = procedure(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 0]
        num_theta_dot_t = procedure(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 1]
        ana_theta_t     = analyt_pend_sol([theta0, theta_dot0], h)[0]
        ana_theta_dot_t = analyt_pend_sol([theta0, theta_dot0], h)[1]

            # Root Mean Square (RMS) Error
        deviation_theta[i]      = np.sqrt(np.sum((num_theta_t - ana_theta_t)**2) / len(ana_theta_t))
        deviation_theta_dot[i]  = np.sqrt(np.sum((num_theta_dot_t - ana_theta_dot_t)**2) / len(ana_theta_dot_t))
        deviation_energy[i]     = np.sqrt(np.sum((calc_energy(theta0, theta_dot0) - calc_energy(num_theta_t, num_theta_dot_t) )**2) / len(num_theta_t) )

    return deviation_theta, deviation_theta_dot, deviation_energy

    # Define step sizes
h_array = np.logspace(-5, 0, 20, dtype=np.float128)  # Logarithmic spacing for better visualization

    # Compute deviation
#deviation_RK2_x, deviation_RK2_v, deviation_RK2_E = calc_deviation(h_array, RK_2)
#deviation_RK4_x, deviation_RK4_v, deviation_RK4_E = calc_deviation(h_array, RK_4)
#deviation_RK6_x, deviation_RK6_v, deviation_RK6_E = calc_deviation(h_array, RK_6)


def comparison_plot( dev_RK2, dev_RK4, dev_RK6, error_attrib, save_as ):
    
    sample_format_dict_1 = {
        "label"      : r"RK2 Error",          
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"RK4 Error",                       
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[6],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  
    sample_format_dict_3 = {
        "label"      : r"RK6 Error",                       
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : "Numerical Global RMS Error in " + error_attrib,
        "x_ax_label"  : r"Step Size $h$ [a.u.]",
        "y_ax_label"  : r"Global RMS Error $\Delta_{RMS}$ [a.u.]"
    }
    
    general_format_dict = standard_format_dict.copy()
    general_format_dict["log_scaling_xy"] = [True, True, 10]
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = h_array,    None,  dev_RK2,  None       
    data_set_2  = h_array,    None,  dev_RK4,  None 
    data_set_3  = h_array,    None,  dev_RK6,  None 
    
    all_data                = data_set_1 + data_set_2 + data_set_3                                
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]

    save_plot = True, save_as                                    
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)


def energy_divergence_plot( h, t_max_new ):
    
    global t_max 
    t_max = t_max_new
    steps       = int(t_max / h) + 1  # Ensure consistency
    t           = np.linspace(0, t_max, steps, dtype=np.float128)
    
    num_theta_t     = RK_2(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 0]
    num_theta_dot_t = RK_2(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 1]
    E_2             = calc_energy(theta0, theta_dot0) - calc_energy(num_theta_t, num_theta_dot_t)

    num_theta_t     = RK_4(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 0]
    num_theta_dot_t = RK_4(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 1]
    E_4             = calc_energy(theta0, theta_dot0) - calc_energy(num_theta_t, num_theta_dot_t)

    num_theta_t     = RK_6(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 0]
    num_theta_dot_t = RK_6(y_dot_simp_pend, [theta0, theta_dot0], h)[:, 1]
    E_6             = calc_energy(theta0, theta_dot0) - calc_energy(num_theta_t, num_theta_dot_t)

    
    sample_format_dict_1 = {
        "label"      : r"RK2 Error",          
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[0],                               
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"RK4 Error",                       
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[6],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  
    sample_format_dict_3 = {
        "label"      : r"RK6 Error",                       
        "fmt"        : '-o', 
        "color"      : sns.color_palette("dark")[9],        
        "markersize" : 4, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1
    }  

    writtings = {
        "title"       : r"Numerical Energy Divergence $\Delta E$ ",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Local Error $\Delta E$ [$J$]"
    }
    
    general_format_dict = standard_format_dict.copy()
    general_format_dict["log_scaling_xy"] = [False, False, 10]
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = t,    None,  E_2,  None       
    data_set_2  = t,    None,  E_4,  None 
    data_set_3  = t,    None,  E_6,  None 
    
    all_data                = data_set_1 + data_set_2 + data_set_3                                
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]

    save_plot = True, "../plots/plot_deviation_diverg.jpg"                                    
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)


#comparison_plot( deviation_RK2_x, deviation_RK4_x, deviation_RK6_x, r"Angle $\theta$", "../plots/plot_deviation_x.jpg" )
#comparison_plot( deviation_RK2_x, deviation_RK4_x, deviation_RK6_x, r"Angle $\theta$", "../plots/plot_deviation_x.jpg" )
#comparison_plot( deviation_RK2_v, deviation_RK4_v, deviation_RK6_v, r"Frequency $\dot{\theta}$", "../plots/plot_deviation_v.jpg" )
#comparison_plot( deviation_RK2_E, deviation_RK4_E, deviation_RK6_E, r"Energy $E$", "../plots/plot_deviation_E.jpg" )

energy_divergence_plot( 0.01, 0.5 )
energy_divergence_plot( 0.01, 0.5 )

print("Numerical Deviation plotted")