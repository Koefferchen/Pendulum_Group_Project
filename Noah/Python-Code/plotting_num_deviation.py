
from ultimate_plotting_v7 import *


    # importing the data generated
data_bsp        = np.loadtxt("../data/data_test_num_solver.txt", skiprows=1 )
h               = data_bsp[ : , 0 ]
params          = data_bsp[ : , 1 ]
deviation_Eul   = data_bsp[ : , 2 ]
deviation_RK4   = data_bsp[ : , 3 ]
deviation_RK6   = data_bsp[ : , 4 ]

x_fit_Eul, y_fit_Eul, a_Eul, b_Eul = linear_fit( np.log(h), np.log(deviation_Eul) )
x_fit_RK4, y_fit_RK4, a_RK4, b_RK4 = linear_fit( np.log(h), np.log(deviation_RK4) )
x_fit_RK6, y_fit_RK6, a_RK6, b_RK6 = linear_fit( np.log(h), np.log(deviation_RK6) )

    # plotting numeric/analytic solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : f"Euler's method ~{a_Eul:.2f}",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : f"Runge-Kutta-4 method ~{a_RK4:.2f}",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : f"Runge-Kutta_6 method ~{a_RK6:.2f}",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_1b = {
        "label"      : None,         
        "fmt"        : '--', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2b = {
        "label"      : None,         
        "fmt"        : '--', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3b = {
        "label"      : None,         
        "fmt"        : '--', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_1b, sample_format_dict_2b, sample_format_dict_3b ]

    writtings = {
        "title"       : r"Deviation of the numerical from the analytical solution",
        "x_ax_label"  : r"log. stepsize $\log{(h)}$ [ ]",
        "y_ax_label"  : r"log. deviation $\log{(\theta_{num} - \theta_{analyt})}$)"
    }
    
    general_format_dict = standard_format_dict
    general_format_dict["log_scaling_xy"] = [False, False, np.e]
    zoom_params= no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = np.log(h), None, np.log(deviation_Eul), None 
    data_set_2  = np.log(h), None, np.log(deviation_RK4), None 
    data_set_3  = np.log(h), None, np.log(deviation_RK6), None 
    data_set_1b = x_fit_Eul, None, y_fit_Eul, None 
    data_set_2b = x_fit_RK4, None, y_fit_RK4, None 
    data_set_3b = x_fit_RK6, None, y_fit_RK6, None 
    
    all_data    = data_set_1 + data_set_2 + data_set_3 + data_set_1b + data_set_2b + data_set_3b                         
    save_plot = True, "../plots/plot_num_deviation.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()

print("Numerical Solvers Deviation plotted")