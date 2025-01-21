
from ultimate_plotting_v5 import *


    # importing the data generated
data_bsp        = np.loadtxt("../data/data_trip_pend.txt", skiprows=1 )
time            = data_bsp[ : , 0 ]
theta1          = data_bsp[ : , 1 ]
theta1_dot      = data_bsp[ : , 2 ]
theta2          = data_bsp[ : , 3 ]
theta2_dot      = data_bsp[ : , 4 ]
theta3          = data_bsp[ : , 5 ]
theta3_dot      = data_bsp[ : , 6 ]


    # theta should always be within [-pi, +pi]
theta1   = (theta1 + np.pi) % (2*np.pi) - np.pi
theta2   = (theta2 + np.pi) % (2*np.pi) - np.pi
theta3   = (theta3 + np.pi) % (2*np.pi) - np.pi


    # plotting double pendulum solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"upper$(\theta_{1})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"middle mass $(\theta_{2})$",    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : r"lower mass $(\theta_{3})$",    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "titel"           : r"Solving the Triple Pendulum",
        "x_beschriftung"  : r"Time $t$ [$s$]",
        "y_beschriftung"  : r"Angle $\theta_{j}$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = time, None, theta1, None 
    data_set_2  = time, None, theta2, None 
    data_set_3  = time, None, theta3, None 
    
    all_data    = data_set_1 + data_set_2 + data_set_3                          
    save_plot = True, "../plots/plot_trip_pend.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend()


def ultimate_plot_pend_phasespace():
    
    sample_format_dict_1 = {
        "label"      : r"upper mass $(\theta_{1},\dot{\theta}_1)$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"middle mass $(\theta_{2},\dot{\theta}_2)$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : r"lower mass $(\theta_{3},\dot{\theta}_3)$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "titel"           : r"Numerical Phasespace of the Triple Pendulum",
        "x_beschriftung"  : r"Angle $\theta_{j}$ [$rad$]",
        "y_beschriftung"  : r"Frequency $\dot{\theta_{j}}$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = theta1, None, theta1_dot, None 
    data_set_2  = theta2, None, theta2_dot, None 
    data_set_3  = theta3, None, theta3_dot, None 
    all_data    = data_set_1 + data_set_2 + data_set_3                            
    save_plot = True, "../plots/plot_trip_pend_phasesspace.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend_phasespace()

print("Triple Pendulum plotted")