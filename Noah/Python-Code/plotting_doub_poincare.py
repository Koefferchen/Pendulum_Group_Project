




from ultimate_plotting_v6 import *


    # importing the data generated
data        = np.loadtxt("../data/data_doub_poincare.txt", skiprows=0 )

params      = data[ : , 3 ]
repitions   = int(params[13])                # how many sets of data pairs there are
E_value     = params[12]
theta2_0    = data[ 9 , 3::4]
cmap = "magma"
cycl_colormap    = sns.color_palette( cmap, n_colors=repitions ) # + sns.color_palette("hls", repitions//2)[::-1]

    # plotting poincare section at ( theta1_0 = 0 ; theta1_dot_0 > 0 ) for the double pendulum
def ultimate_plot_pend():
    
    sample_format_dict = {
        "label"      : None,         
        "fmt"        : 'o', 
        "color"      : "",                            
        "markersize" : 1, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }

    writtings = {
        "titel"           : r"The Poincare Section of the double pendulum for ($E = {:.2f}$)".format(E_value),
        "x_beschriftung"  : r"Angle $\theta_2$ [$rad$]",
        "y_beschriftung"  : r"Frequency $\dot {\theta}_2$ [$s^{-1}$]"
    }
    
    general_format_dict = standard_format_dict
    zoom_params = no_zooming
    
    colorbar_params = {
        "do_cbar"       : True,
        "position"      : [0.92, 0.15],
        "size"          : [0.03, 0.70],
        "scale_range"   : [0.0, theta2_0[::-1][0]],
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
                           
    save_plot = True, "../plots/plot_doub_poincare.jpg"                                      
        
    ultimate_plot_advanced( all_data, writtings, zoom_params, colorbar_params, save_plot, all_sample_format_dicts, general_format_dict )
ultimate_plot_pend()
ultimate_plot_pend()

print("Double Pendulum Poincare plotted")