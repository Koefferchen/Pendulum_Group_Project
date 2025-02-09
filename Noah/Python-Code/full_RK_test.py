from ultimate_plotting_v7 import *

g_grav      = 9.81
length      = 1.0
mass        = 1.0
t_max       = 1.0
theta0      = 0.1 * np.pi
theta_dot0  = 0.0

    # right hand side of the differential equation for the simple pendulum
def y_dot_simp_pend(t, y):

    dy = np.zeros(2)
    dy[0] = y[1]
    dy[1] = -g_grav/length * y[0]
    
    return dy

    # full runge kutta method
def RK_2( ODE, y0, h ):
    steps = round(t_max/h)
        # t is the time array
    t       = np.linspace(0, t_max, steps)
        # y[i] is the state of the system at time t[i] 
    y       = np.zeros((steps, len(y0)))
        # y[0] is the initial state of the system
    y[0]    = y0
    for i in range(0, steps-1):
        k1 = h * ODE(t[i], y[i])
        k2 = h * ODE(t[i] + h/2, y[i] + k1/2)
        y[i+1] = y[i] + 0.5 * (k1 + k2)
        # returns array of system state y(t)  
    return y

    # calculates analytical solution to simple pendulum (small angle approx)
def analyt_pend_sol( y0, h ):
    steps = round( t_max/h )
    t           = np.linspace(0, t_max, steps)
    theta_t     = +y0[0] * np.cos(np.sqrt(g_grav/length) * t) + y0[1]/np.sqrt(g_grav/length) * np.sin(np.sqrt(g_grav/length) * t)
    theta_dot_t = -y0[0] * np.sqrt(g_grav/length) * np.sin(np.sqrt(g_grav/length) * t) + y0[1] * np.cos(np.sqrt(g_grav/length) * t)
    return theta_t, theta_dot_t

    # calculates the deviation between analyt and numeric method for h in h_array
def calc_deviation( h_array, procedure ):
    
    deviation   = np.zeros( len(h_array) )


    for i in range( len(h_array) ):
        h = h_array[i]
        steps       = round( t_max/h)

        num_theta_t     = procedure(y_dot_simp_pend, [theta0, theta_dot0], h)[:,0]
        num_theta_dot_t = procedure(y_dot_simp_pend, [theta0, theta_dot0], h)[:,1]

        ana_theta_t, ana_theta_dot_t = analyt_pend_sol([theta0, theta_dot0], h)

        #num_theta_t = np.mod( num_theta_t, 2 * np.pi)
        #ana_theta_t = np.mod( ana_theta_t, 2 * np.pi)

        print(num_theta_dot_t-ana_theta_dot_t)

        deviation[i] = np.sqrt( np.sum( (num_theta_t - ana_theta_t)**2 ) / len(ana_theta_t) )
        #deviation[i] = np.sum( abs( num_theta_t - ana_theta_t ) )/len(num_theta_t)
        #deviation[i] = abs( num_theta_t[steps-1] - ana_theta_t[steps-1] ) 

    return deviation


h_array         = np.linspace( 0.1, 0.001, 100 )
deviation_RK2   = calc_deviation( h_array, RK_2 )



# ------------------------- Plotting -------------------------


#x_fit_RK2, y_fit_RK2, a_RK2, b_RK2 =  linear_fit( np.log(h_array), np.log(deviation_RK2) )
#x_fit_RK4, y_fit_RK4, a_RK4, b_RK4 =  linear_fit( np.log(h), np.log(deviation_RK4) )
#x_fit_RK6, y_fit_RK6, a_RK6, b_RK6 =  linear_fit( np.log(h), np.log(deviation_RK6) )

    # plotting numeric/analytic solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : f"Runge-Kutta-2 method ",#~{a_RK2:.2f}",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }


    all_sample_format_dicts = [ sample_format_dict_1]#, sample_format_dict_2b, sample_format_dict_3b ]

    writtings = {
        "title"       : r"Deviation of the numerical from the analytical solution",
        "x_ax_label"  : r"log. stepsize $\log{(h)}$ [ ]",
        "y_ax_label"  : r"log. deviation $\log{(\theta_{num} - \theta_{analyt})}$)"
    }
    
    general_format_dict = standard_format_dict
    general_format_dict["log_scaling_xy"] = [False, False, 10]
    zoom_params= no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
#    data_set_1  = np.log(h_array), None, np.log(deviation_RK2), None 
#    data_set_2  = np.log(h), None, np.log(deviation_RK4), None 
#    data_set_3  = np.log(h), None, np.log(deviation_RK6), None 
#    data_set_1b = x_fit_RK2, None, y_fit_RK2, None 
#    data_set_2b = x_fit_RK4, None, y_fit_RK4, None 
#    data_set_3b = x_fit_RK6, None, y_fit_RK6, None 

    data_set_1  = h_array, None, deviation_RK2, None 
#    data_set_2  = h, None, deviation_RK4, None 
#    data_set_3  = h, None, deviation_RK6, None 
#    data_set_1b = x_fit_RK2, None, y_fit_RK2, None 
#    data_set_2b = x_fit_RK4, None, y_fit_RK4, None 
#    data_set_3b = x_fit_RK6, None, y_fit_RK6, None 

    
    all_data    = data_set_1  #+ data_set_2b + data_set_3b                         
    save_plot = True, "../plots/plot_num_deviation.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()

print("Numerical Solvers Deviation plotted")

