
from ultimate_plotting_v7 import *


    # names for plots:
identifier      = "06"
saveas_pend     = "../plots/plot_doub_"+identifier+"_pendl.jpg"
saveas_phase    = "../plots/plot_doub_"+identifier+"_phase.jpg"
saveas_pos      = "../plots/plot_doub_"+identifier+"_posit.jpg"

    # importing the data generated
data_bsp        = np.loadtxt("../data/data_doub_pend.txt", skiprows=1 )
time            = data_bsp[ : , 0 ]
theta1          = data_bsp[ : , 1 ]
theta1_dot      = data_bsp[ : , 2 ]
theta2          = data_bsp[ : , 3 ]
theta2_dot      = data_bsp[ : , 4 ]
params          = data_bsp[ : , 5 ]

    # unpack the parameters of the simulation
t_end       = params[0]
h           = params[1]
g_grav      = params[2]
mass_1      = params[3]
mass_2      = params[4]
length_1    = params[5]
length_2    = params[6]
theta1_0    = params[7]
theta1_dot_0= params[8]
theta2_0    = params[9]
theta2_dot_0= params[10]

m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n\n" 
i_label1 = f"$\\theta_1(0) = {theta1_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.2f}\\pi$ \n"
i_label2 = f"$\\theta_2(0) = {theta2_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_2(0) = {theta2_dot_0/np.pi:.2f}\\pi$ \n"

    # theta should always be within [-pi, +pi]
theta1   = (theta1 + np.pi) % (2*np.pi) - np.pi
theta2   = (theta2 + np.pi) % (2*np.pi) - np.pi


    # plotting double pendulum solution
def ultimate_plot_pend():
    
    sample_format_dict_1 = {
        "label"      : r"upper mass $(\theta_{1})$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"lower mass $(\theta_{2})$",    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]
    
    writtings = {
        "title"       : r"Solving the Double Pendulum",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta_{j}$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label1+i_label2)
    }
    
    data_set_1  = time, None, theta1, None 
    data_set_2  = time, None, theta2, None 
    
    all_data    = data_set_1 + data_set_2                            
    save_plot = True, saveas_pend                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()


def ultimate_plot_pend_phasespace():
    
    sample_format_dict_1 = {
        "label"      : r"upper mass $(\theta_{1},\dot{\theta}_2)$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"lower mass $(\theta_{2},\dot{\theta}_2)$",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2 ]
    
    writtings = {
        "title"       : r"Numerical Phasespace of the Double Pendulum",
        "x_ax_label"  : r"Angle $\theta_{j}$ [$rad$]",
        "y_ax_label"  : r"Frequency $\dot{\theta_{j}}$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label1+i_label2)
    }
    
    
    data_set_1  = theta1, None, theta1_dot, None 
    data_set_2  = theta2, None, theta2_dot, None 
    all_data    = data_set_1 + data_set_2                            
    save_plot = True, saveas_phase                                     
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend_phasespace()

def ultimate_plot_pend_positionspace():
    
    sample_format_dict_1 = {
        "label"      : None,         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    all_sample_format_dicts = [ sample_format_dict_1 ]
    
    writtings = {
        "title"       : r"Positional trajectory of the Double Pendulum",
        "x_ax_label"  : r"Angle $\theta_{1}$ [$rad$]",
        "y_ax_label"  : r"Angle $\theta_{2}$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label1+i_label2)
    }
    
    
    data_set_1  = theta1, None, theta2, None 
    all_data    = data_set_1                            
    save_plot = True, saveas_pos                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend_positionspace()

print("Double Pendulum plotted")