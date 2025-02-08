
from ultimate_plotting_v7 import *


    # importing the data generated
data_bsp        = np.loadtxt("../data/data_simp_pend.txt", skiprows=1 )
time            = data_bsp[ : , 0 ]
theta_num       = data_bsp[ : , 1 ]
theta_dot_num   = data_bsp[ : , 2 ]
theta_ana       = data_bsp[ : , 3 ]

    # theta should always be within [-pi, +pi]
theta_num   = (theta_num + np.pi) % (2*np.pi) - np.pi


    # plotting numeric/analytic solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"numerical soultion ",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"analytical soultion",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[1],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]
    
    writtings = {
        "title"           : r"Solving $\ddot{\theta} = - \sqrt{ \frac{g}{l} } \cdot \sin{(\theta)}$",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = time, None, theta_num, None 
    data_set_2  = time, None, theta_ana, None 
    
    all_data    = data_set_1 + data_set_2                            
    save_plot = True, "../plots/plot_simp_pend.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()


def ultimate_plot_pend_phasespace():
    
    sample_format_dict_1 = {
        "label"      : r"numerical soultion",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    
    all_sample_format_dicts = [ sample_format_dict_1]
    
    writtings = {
        "title"           : r"Numerical Phasespace for $\ddot{\theta} = - \omega^{2} \cdot \theta$",
        "x_ax_label"  : r"Angle $\theta$ [$rad$]",
        "y_ax_label"  : r"Frequency $\dot{\theta}$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = theta_num, None, theta_dot_num, None 
    all_data    = data_set_1                            
    save_plot = True, "../plots/plot_simp_pend_phasesspace.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend_phasespace()

print("Simple Pendulum plotted")