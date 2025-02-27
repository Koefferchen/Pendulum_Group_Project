
from ultimate_plotting_v7 import *


    # names for plots:
identifier  = "01"
saveas      = "../plots/plot_doub_poinc_"+identifier+".jpg" 

    # importing the data generated
data        = np.loadtxt("../data/data_doub_poincare.txt", skiprows=0 )

params      = data[ : , 3 ]
repitions   = int(params[13])                # how many sets of data pairs there are
E_value     = params[12]
theta2_0    = data[ 10 , 3::4]
cmap = "magma"
cycl_colormap    = sns.color_palette( cmap, n_colors=repitions ) # + sns.color_palette("hls", repitions//2)[::-1]

    # unpack the parameters of the simulation
t_end       = params[1]
h           = params[2]
g_grav      = params[3]
mass_1      = params[4]
mass_2      = params[5]
length_1    = params[6]
length_2    = params[7]
theta1_0    = params[8]
theta1_dot_0= params[9]

m_label = f"$m_1 = {mass_1:.2f}$kg \n$m_2 = {mass_2:.2f}$kg \n"
l_label = f"$l_1 = {length_1:.2f}$m \n$l_2 = {length_2:.2f}$m \n" 
i_label1 = f"$\\theta_1(0) = {theta1_0/np.pi:.2f}\\pi$ \n$\\dot{{\\theta}}_1(0) = {theta1_dot_0/np.pi:.2f}\\pi$ \n"


    # plotting poincare section at ( theta1_0 = 0 ; theta1_dot_0 > 0 ) for the double pendulum
def ultimate_plot_pend():
    
    sample_format_dict = {
        "label"      : None,         
        "fmt"        : '.', 
        "color"      : "",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "title"       : f"{len(theta2_0):.0f} Poincare Sections for the Double Pendulum  at E = {E_value:.2f}J",
        "x_ax_label"  : r"Angle $\theta_2$ [$rad$]",
        "y_ax_label"  : r"Frequency $\dot {\theta}_2$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params = no_zooming
    extra_label         = {
        "do_label"  :   False,
        "position"  :   [1.16, 0.97],
        "font_size" :   12,
        "content"   :   (m_label+l_label+i_label1)
    }
    
    colorbar_params = {
        "do_cbar"       : True,
        "position"      : [0.92, 0.15],
        "size"          : [0.03, 0.70],
        "scale_range"   : [min(theta2_0), max(theta2_0)],
        "title"         : r"$\theta_2 (t=0)$",
        "colormap"      : cmap
    }

    all_sample_format_dicts = []
    all_data = ()

    for i in range(repitions):
        
        sample_format_dict["color"] = cycl_colormap[i] 
        sample_size     = int(data[ 0 , 1 +4*i ])
        theta2          = data[ 1:sample_size , 1 +4*i ]
        theta2          = (theta2 + np.pi) % (2*np.pi) - np.pi   
        theta2_dot      = data[ 1:sample_size , 2 +4*i ]
        data_set        = theta2, None, theta2_dot, None    

        all_sample_format_dicts.append(sample_format_dict.copy())
        all_data = all_data + data_set
                           
    save_plot = True, saveas                                     
        
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)
ultimate_plot_pend()
ultimate_plot_pend()

print("Double Pendulum Poincare plotted")