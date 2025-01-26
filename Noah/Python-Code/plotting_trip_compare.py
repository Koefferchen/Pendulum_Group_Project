
from ultimate_plotting_v5 import *


    # importing the data generated
data_31         = np.loadtxt("../data/data_trip_compare_3-1.txt", skiprows=1 )
time            = data_31[ : , 0 ]
simp_theta1     = data_31[ : , 1 ]
trip_theta1     = data_31[ : , 2 ]


    # theta should always be within [-pi, +pi]
simp_theta1   = (simp_theta1 + np.pi) % (2*np.pi) - np.pi
trip_theta1   = (trip_theta1 + np.pi) % (2*np.pi) - np.pi


    # plotting simple and triple pendulum solution
def ultimate_plot_31():
    
    sample_format_dict_1 = {
        "label"      : r"Simple Pendulum ($\theta$)",         
        "fmt"        : '-', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"Triple Pendulum ($\theta_1$)",    
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
   
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]
    
    writtings = {
        "titel"           : r"Triple Pendulum for ($m_2+m_3 \ll m_1$, $l_2+l_3 \ll l_1$)",
        "x_beschriftung"  : r"Time $t$ [$s$]",
        "y_beschriftung"  : r"Angle $\theta$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = time, None, simp_theta1, None 
    data_set_2  = time, None, trip_theta1, None 
    
    all_data    = data_set_1 + data_set_2                           
    save_plot = True, "../plots/plot_trip_compare_3-1.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_31()
ultimate_plot_31()




    # importing the data generated
data_32         = np.loadtxt("../data/data_trip_compare_3-2.txt", skiprows=1 )
time            = data_32[ : , 0 ]
doub_theta1     = data_32[ : , 1 ]
doub_theta2     = data_32[ : , 2 ]
trip_theta1     = data_32[ : , 3 ]
trip_theta2     = data_32[ : , 4 ]

    # theta should always be within [-pi, +pi]
doub_theta1   = (doub_theta1 + np.pi) % (2*np.pi) - np.pi
doub_theta2   = (doub_theta2 + np.pi) % (2*np.pi) - np.pi
trip_theta1   = (trip_theta1 + np.pi) % (2*np.pi) - np.pi
trip_theta2   = (trip_theta2 + np.pi) % (2*np.pi) - np.pi

    # plotting double and triple pendulum solution
def ultimate_plot_32():
    
    sample_format_dict_1 = {
        "label"      : r"Simple Pendulum ($\theta_1$)",         
        "fmt"        : '-', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"Simple Pendulum ($\theta_2$)",    
        "fmt"        : '-', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : r"Triple Pendulum ($\theta_1$)",    
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_4 = {
        "label"      : r"Triple Pendulum ($\theta_2$)",    
        "fmt"        : '--', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4 ]
    
    writtings = {
        "titel"           : r"Triple Pendulum for ($m_3 \ll m_2$, $l_3 \ll l_2$)",
        "x_beschriftung"  : r"Time $t$ [$s$]",
        "y_beschriftung"  : r"Angle $\theta_j$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1  = time, None, doub_theta1, None 
    data_set_2  = time, None, doub_theta2, None 
    data_set_3  = time, None, trip_theta1, None 
    data_set_4  = time, None, trip_theta2, None 
    
    all_data    = data_set_1 + data_set_2 + data_set_3 + data_set_4                           
    save_plot = True, "../plots/plot_trip_compare_3-2.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_32()


print("Triple Pendulum Comparison plotted")