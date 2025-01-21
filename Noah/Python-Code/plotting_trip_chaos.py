
from ultimate_plotting_v5 import *


    # importing the data generated
data_a          = np.loadtxt("../data/data_trip_chaos_a.txt", skiprows=1 )
time            = data_a[ : , 0 ]
theta1a         = data_a[ : , 1 ]
theta2a         = data_a[ : , 3 ]
theta3a         = data_a[ : , 5 ]
theta1a_dot     = data_a[ : , 2 ]
theta2a_dot     = data_a[ : , 4 ]
theta3a_dot     = data_a[ : , 6 ]

data_b          = np.loadtxt("../data/data_trip_chaos_b.txt", skiprows=1 )
theta1b         = data_b[ : , 1 ]
theta2b         = data_b[ : , 3 ]
theta3b         = data_b[ : , 5 ]   
theta1b_dot     = data_b[ : , 2 ]
theta2b_dot     = data_b[ : , 4 ]
theta3b_dot     = data_b[ : , 6 ]

data_c          = np.loadtxt("../data/data_trip_chaos_c.txt", skiprows=1 )
theta1c         = data_c[ : , 1 ]
theta2c         = data_c[ : , 3 ]
theta3c         = data_c[ : , 5 ]   
theta1c_dot     = data_c[ : , 2 ]
theta2c_dot     = data_c[ : , 4 ]
theta3c_dot     = data_c[ : , 6 ]



    # theta should always be within [-pi, +pi]
theta1a   = (theta1a + np.pi) % (2*np.pi) - np.pi
theta2a   = (theta2a + np.pi) % (2*np.pi) - np.pi
theta3a   = (theta3a + np.pi) % (2*np.pi) - np.pi
theta1b   = (theta1b + np.pi) % (2*np.pi) - np.pi
theta2b   = (theta2b + np.pi) % (2*np.pi) - np.pi
theta3b   = (theta3b + np.pi) % (2*np.pi) - np.pi
theta1c   = (theta1c + np.pi) % (2*np.pi) - np.pi
theta2c   = (theta2c + np.pi) % (2*np.pi) - np.pi
theta3c   = (theta3c + np.pi) % (2*np.pi) - np.pi


    # plotting double pendulum solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"lower mass $(\theta_{3a})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1,
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"lower mass $(\theta_{3b})$",    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : r"lower mass $(\theta_{3c})$",    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                    
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "titel"           : r"Triple Pendulum: Sensitivity to Initial Conditions",
        "x_beschriftung"  : r"Time $t$ [$s$]",
        "y_beschriftung"  : r"Angle $\theta_{3j}$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = time, None, theta3a, None 
    data_set_2  = time, None, theta3b, None 
    data_set_3  = time, None, theta3c, None  
    
    all_data    = data_set_1 + data_set_2 + data_set_3                           
    save_plot = True, "../plots/plot_trip_chaos.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend()


def ultimate_plot_pend_phasespace():
    
    sample_format_dict_1 = {
        "label"      : r"lower mass $(\theta_{3a},\dot{\theta}_{3a})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"lower mass $(\theta_{3b},\dot{\theta}_{3b})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    sample_format_dict_3 = {
        "label"      : r"lower mass $(\theta_{3c},\dot{\theta}_{3c})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "titel"           : r"Triple Pendulum: Sensitivity to Initial Conditions (Phase Space)",
        "x_beschriftung"  : r"Angle $\theta_{3j}$ [$rad$]",
        "y_beschriftung"  : r"Frequency $\dot{\theta_{3j}}$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = theta3a, None, theta3a_dot, None 
    data_set_2  = theta3b, None, theta3b_dot, None 
    data_set_3  = theta3c, None, theta3c_dot, None 
    all_data    = data_set_1 + data_set_2 + data_set_3                             
    save_plot = True, "../plots/plot_trip_chaos_phasesspace.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend_phasespace()

print("Triple Pendulum Chaos plotted")