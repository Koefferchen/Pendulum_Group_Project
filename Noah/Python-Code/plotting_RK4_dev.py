
from ultimate_plotting_v7 import *


    # importing the data generated
data_bsp        = np.loadtxt("../data/data_test_num_solver.txt", skiprows=1 )
h               = data_bsp[ : , 0 ]
params          = data_bsp[ : , 1 ]
deviation_Eul   = data_bsp[ : , 2 ]
deviation_RK4   = data_bsp[ : , 3 ]
deviation_RK6   = data_bsp[ : , 4 ]


    # plotting numeric/analytic solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"Euler's method",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"Runge-Kutta-4 method",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : r"Runge-Kutta_6 method",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]

    writtings = {
        "title"           : r"Deviation of the numerical from the analytical solution",
        "x_ax_label"  : r"stepsize $h$ [ ]",
        "y_ax_label"  : r"absolute deviation ($\theta_{num} - \theta_{analyt}$)"
    }
    
    general_format_dict = standard_format_dict
    general_format_dict["log_scaling_xy"] = [True, True, 10]
    zoom_params= no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = h, None, deviation_Eul, None 
    data_set_2  = h, None, deviation_RK4, None 
    data_set_3  = h, None, deviation_RK6, None 
    
    all_data    = data_set_1+data_set_2+data_set_3                             
    save_plot = True, "../plots/plot_num_deviation.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()
