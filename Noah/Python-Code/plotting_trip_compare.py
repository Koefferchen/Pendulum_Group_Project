
from ultimate_plotting_v7 import *


    # importing the data generated
data_31         = np.loadtxt("../data/data_trip_compare_3-1.txt", skiprows=1 )
time            = data_31[ : , 0 ]
simp_theta1     = data_31[ : , 1 ]
trip_theta1     = data_31[ : , 2 ]
params_31       = data_31[ : , 3 ]


    # unpack the parameters of the simulation
t_end       = params_31[0]
h           = params_31[1]
g_grav      = params_31[2]
mass_1      = params_31[3]
mass_2      = params_31[4]
mass_3      = params_31[5]
length_1    = params_31[6]
length_2    = params_31[7]
length_3    = params_31[8]
theta1_0    = params_31[9]
theta1_dot_0= params_31[10]

m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n$m_3 = {mass_3:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n$l_3 = {length_3:.2f}$m  \n\n" 
i_label = f"$\\theta_1(0) = {theta1_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.2f}\\pi$"


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
        "title"       : r"Triple Pendulum for ($m_2+m_3 \ll m_1$, $l_2+l_3 \ll l_1$)",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label)
    }

    data_set_1  = time, None, simp_theta1, None 
    data_set_2  = time, None, trip_theta1, None 
    
    all_data    = data_set_1 + data_set_2                           
    save_plot = True, "../plots/plot_trip_compare_3-1.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)

ultimate_plot_31()
ultimate_plot_31()


# ---------------------------------------------------------------------------------


    # importing the data generated
data_32         = np.loadtxt("../data/data_trip_compare_3-2.txt", skiprows=1 )
time            = data_32[ : , 0 ]
doub_theta1     = data_32[ : , 1 ]
doub_theta2     = data_32[ : , 2 ]
trip_theta1     = data_32[ : , 3 ]
trip_theta2     = data_32[ : , 4 ]
params_32       = data_32[ : , 5 ]


    # unpack the parameters of the simulation
t_end       = params_32[0]
h           = params_32[1]
g_grav      = params_32[2]
mass_1      = params_32[3]
mass_2      = params_32[4]
mass_3      = params_32[5]
length_1    = params_32[6]
length_2    = params_32[7]
length_3    = params_32[8]
theta1_0    = params_32[9]
theta1_dot_0= params_32[10]
theta2_0    = params_32[11]
theta2_dot_0= params_32[12]

m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n$m_3 = {mass_3:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n$l_3 = {length_3:.2f}$m  \n\n" 
i_label = f"$\\theta_1(0) = {theta1_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.2f}\\pi$ \n$\\theta_2(0) = {theta2_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_2(0) = {theta2_dot_0/np.pi:.2f}\\pi$"

    # theta should always be within [-pi, +pi]
doub_theta1   = (doub_theta1 + np.pi) % (2*np.pi) - np.pi
doub_theta2   = (doub_theta2 + np.pi) % (2*np.pi) - np.pi
trip_theta1   = (trip_theta1 + np.pi) % (2*np.pi) - np.pi
trip_theta2   = (trip_theta2 + np.pi) % (2*np.pi) - np.pi

    # plotting double and triple pendulum solution
def ultimate_plot_32():
    
    sample_format_dict_1 = {
        "label"      : r"Double Pendulum ($\theta_1$)",         
        "fmt"        : '-', 
        "color"      : "grey",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : r"Double Pendulum ($\theta_2$)",    
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
    writtings = {
        "title"       : r"Triple Pendulum for ($m_3 \ll m_2$, $l_3 \ll l_2$)",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta_j$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label)
    }


    
    data_set_1  = time, None, doub_theta1, None 
    data_set_2  = time, None, doub_theta2, None 
    data_set_3  = time, None, trip_theta1, None 
    data_set_4  = time, None, trip_theta2, None 
    
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4 ]
    all_data    = data_set_1 + data_set_2 + data_set_3 + data_set_4                           
    save_plot = True, "../plots/plot_trip_compare_3-2.jpg"                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_32()


print("Triple Pendulum Comparison plotted")