
from ultimate_plotting_v7 import *


    # names for plots:
identifier      = "03"
saveas_pend     = "../plots/plot_trip_"+identifier+"_pendl.jpg"
saveas_phase    = "../plots/plot_trip_"+identifier+"_phase.jpg"
saveas_pos_xy   = "../plots/plot_trip_"+identifier+"_posit_xy.jpg"

    # importing the data generated
data        = np.loadtxt("../data/data_trip_pend.txt", skiprows=1 )
time            = data[ : , 0 ]
theta1          = data[ : , 1 ]
theta1_dot      = data[ : , 2 ]
theta2          = data[ : , 3 ]
theta2_dot      = data[ : , 4 ]
theta3          = data[ : , 5 ]
theta3_dot      = data[ : , 6 ]
params          = data[ : , 7 ]

    # unpack the parameters of the simulation
t_end       = params[0]
h           = params[1]
g_grav      = params[2]
mass_1      = params[3]
mass_2      = params[4]
mass_3      = params[5]
length_1    = params[6]
length_2    = params[7]
length_3    = params[8]
theta1_0    = params[9]
theta1_dot_0= params[10]
theta2_0    = params[11]
theta2_dot_0= params[12]
theta3_0    = params[13]
theta3_dot_0= params[14]

m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n$m_3 = {mass_3:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n$l_3 = {length_3:.2f}$m  \n\n" 
i_label1 = f"$\\theta_1(0) = {theta1_0/np.pi:.3f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.3f}\\pi$ \n"
i_label2 = f"$\\theta_2(0) = {theta2_0/np.pi:.3f}\\pi$ \n$\\dot{{\\theta}}_2(0) = {theta2_dot_0/np.pi:.3f}\\pi$ \n"
i_label3 = f"$\\theta_3(0) = {theta3_0/np.pi:.3f}\\pi$ \n$\\dot{{\\theta}}_3(0) = {theta3_dot_0/np.pi:.3f}\\pi$ \n"

    # theta should always be within [-pi, +pi]
theta1   = (theta1 + np.pi) % (2*np.pi) - np.pi
theta2   = (theta2 + np.pi) % (2*np.pi) - np.pi
theta3   = (theta3 + np.pi) % (2*np.pi) - np.pi


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
        "title"       : r"Solving the Triple Pendulum",
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
        "content"   :   (m_label+l_label+i_label1+i_label2+i_label3)
    }

    data_set_1  = time, None, theta1, None 
    data_set_2  = time, None, theta2, None 
    data_set_3  = time, None, theta3, None 
    
    all_data    = data_set_1 + data_set_2 + data_set_3                          
    save_plot = True, saveas_pend                                     
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
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
        "title"       : r"Numerical Phasespace of the Triple Pendulum",
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
        "content"   :   (m_label+l_label+i_label1+i_label2+i_label3)
    }
    
    data_set_1  = theta1, None, theta1_dot, None 
    data_set_2  = theta2, None, theta2_dot, None 
    data_set_3  = theta3, None, theta3_dot, None 
    all_data    = data_set_1 + data_set_2 + data_set_3                            
    save_plot = True, saveas_phase                                   
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend_phasespace()



    # helper function to produce diff. plots with little effort
def ret_array( theta_number ):
    if( theta_number == 1 ):
        return theta1
    elif( theta_number == 2 ):
        return theta2
    elif( theta_number == 3 ):
        return theta3
    else:
        print("Input of 'ret_array() in 'plotting_trip_pend.py' must be 1,2,3 ")
    
def ultimate_plot_pend_positionspace( thetas ):


    
    sample_format_dict_1 = {
        "label"      : r"($\theta_"+str(thetas[0])+r"(t), \theta_"+str(thetas[1])+"(t)$)",         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1 ]
    
    writtings = {
        "title"       : r"Positional Trajectory of the Triple Pendulum ("+str(thetas[0])+","+str(thetas[1])+")",
        "x_ax_label"  : r"Angle $\theta_{"+str(thetas[0])+r"}$ [$rad$]",
        "y_ax_label"  : r"Angle $\theta_{"+str(thetas[1])+r"}$ [$rad$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = {
        "do_label"  :   True,
        "position"  :   [1.03, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label1+i_label2+i_label3)
    }

    data_set_1  = ret_array(thetas[0]), None, ret_array(thetas[1]), None 
    all_data    = data_set_1                             
    save_plot = True, saveas_pos_xy.replace("xy",str(thetas[0])+str(thetas[1]) )                                  
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)

ultimate_plot_pend_positionspace([1,2])
ultimate_plot_pend_positionspace([2,3])
ultimate_plot_pend_positionspace([3,1])

print("Triple Pendulum plotted")