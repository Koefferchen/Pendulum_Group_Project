
from ultimate_plotting_v7 import *
import sys

    # making sure an argument is passed
if len(sys.argv) != 2:
    print("Usage: python3 script.py <1 or 0>")
    sys.exit(1)


usage = int(sys.argv[1])
if( usage == 0 ):
        # plotting sensitivity on initial conditions 
    title1  = r"Triple Pendulum: Sensitivity to Initial Conditions"
    title2  = r"Triple Pendulum: Sensitivity to Initial Conditions (Phase Space)"
    save_as1= "../plots/plot_trip_chaos.jpg" 
    save_as2= "../plots/plot_trip_chaos_phasesspace.jpg"
    label1  = r"lower mass $(\theta_{3a})$"
    label2  = r"lower mass $(\theta_{3b})$"
    label3  = r"lower mass $(\theta_{3c})$"
    info    = "Triple Pendulum Chaos plotted"
    file_a  = "../data/data_trip_chaos_a.txt"
    file_b  = "../data/data_trip_chaos_b.txt"
    file_c  = "../data/data_trip_chaos_c.txt"
elif( usage == 1 ):
        # plotting sensitivity on numerical solver
    title1  = r"Triple Pendulum: Sensitivity to Numerical Solvers"
    title2  = r"Triple Pendulum: Sensitivity to Numerical Solvers (Phase Space)"
    save_as1= "../plots/plot_trip_nums.jpg" 
    save_as2= "../plots/plot_trip_nums_phasesspace.jpg"
    label1  = r"Euler Procedure $\mathcal{O}(h^2)$"
    label2  = r"Runge-Kutta-4 Procedure $\mathcal{O}(h^4)$"
    label3  = r"Runge-Kutta-6 Procedure $\mathcal{O}(h^6)$"
    info    = "Triple Pendulum Numerical Solvers plotted"
    file_a  = "../data/data_trip_Euler.txt"
    file_b  = "../data/data_trip_RuKu4.txt"
    file_c  = "../data/data_trip_RuKu6.txt"
else:
    print("Could not resolve argument in 'plotting_trip_chaos.py'.")
    sys.exit(1)


    # importing the data generated
data_a          = np.loadtxt(file_a, skiprows=1 )
time            = data_a[ : , 0 ]

theta1a         = data_a[ : , 1 ]
theta2a         = data_a[ : , 3 ]
theta3a         = data_a[ : , 5 ]
theta1a_dot     = data_a[ : , 2 ]
theta2a_dot     = data_a[ : , 4 ]
theta3a_dot     = data_a[ : , 6 ]
params          = data_a[ : , 7 ]

data_b          = np.loadtxt(file_b, skiprows=1 )
theta1b         = data_b[ : , 1 ]
theta2b         = data_b[ : , 3 ]
theta3b         = data_b[ : , 5 ]   
theta1b_dot     = data_b[ : , 2 ]
theta2b_dot     = data_b[ : , 4 ]
theta3b_dot     = data_b[ : , 6 ]

data_c          = np.loadtxt(file_c, skiprows=1 )
theta1c         = data_c[ : , 1 ]
theta2c         = data_c[ : , 3 ]
theta3c         = data_c[ : , 5 ]   
theta1c_dot     = data_c[ : , 2 ]
theta2c_dot     = data_c[ : , 4 ]
theta3c_dot     = data_c[ : , 6 ]

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

    # create extra_label
m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n$m_3 = {mass_3:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n$l_3 = {length_3:.2f}$m  \n\n" 
i_label1 = f"$\\theta_1(0) = {theta1_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.2f}\\pi$ \n"
i_label2 = f"$\\theta_2(0) = {theta2_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_2(0) = {theta2_dot_0/np.pi:.2f}\\pi$ \n"
i_label3 = f"$\\theta_3(0) = {theta3_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_3(0) = {theta3_dot_0/np.pi:.2f}\\pi$ \n"


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
        "label"      : label1,         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1,
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : label2,    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_3 = {
        "label"      : label3,    
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                    
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "title"       : title1,
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta_{3j}$ [$rad$]"
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
    
    data_set_1  = time, None, theta3a, None 
    data_set_2  = time, None, theta3b, None 
    data_set_3  = time, None, theta3c, None  
    
    all_data    = data_set_1 + data_set_2 + data_set_3                           
    save_plot = True, save_as1                                     
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()


def ultimate_plot_pend_phasespace():
    
    sample_format_dict_1 = {
        "label"      : label1,       
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    sample_format_dict_2 = {
        "label"      : label2,         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[6],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
    
    sample_format_dict_3 = {
        "label"      : label3,         
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],                            
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    all_sample_format_dicts = [ sample_format_dict_1, sample_format_dict_2, sample_format_dict_3 ]
    
    writtings = {
        "title"       : title2,
        "x_ax_label"  : r"Angle $\theta_{3j}$ [$rad$]",
        "y_ax_label"  : r"Frequency $\dot{\theta_{3j}}$ [$s^{-1}$]"
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
    
    data_set_1  = theta3a, None, theta3a_dot, None 
    data_set_2  = theta3b, None, theta3b_dot, None 
    data_set_3  = theta3c, None, theta3c_dot, None 
    all_data    = data_set_1 + data_set_2 + data_set_3                             
    save_plot = True, save_as2                                      
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend_phasespace()

print(info)