
from ultimate_plotting_v5 import *


    # importing the data generated
data_bsp        = np.loadtxt("../data/data_test_RK4.txt", skiprows=1 )
h               = data_bsp[ : , 0 ]
deviation       = data_bsp[ : , 1 ]
params          = data_bsp[ : , 2 ]

    # plotting numeric/analytic solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"RK4 method",         
        "fmt"        : 'o', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1 ]

    writtings = {
        "titel"           : r"Deviation of the numerical from the analytical solution",
        "x_beschriftung"  : r"stepsize $h$ []",
        "y_beschriftung"  : r"absolute deviation $\theta_{num} - \theta_{analyt}$"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    general_format_dict["log_scaling_xy"] = [True, True, 10]

    
    data_set_1  = h, None, deviation, None 
    
    all_data    = data_set_1                             
    save_plot = True, "../plots/plot_RK4_deviation.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend()
ultimate_plot_pend()
