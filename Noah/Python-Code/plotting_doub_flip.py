from ultimate_plotting_v7 import *

data_1 = np.loadtxt("../data/data_doub_flip_1.txt")
data_2 = np.loadtxt("../data/data_doub_flip_2.txt")

params          = data_1[ 1: , 0 ]
t_flip1_matrix  = np.transpose( data_1[ : , 1: ] )
t_flip2_matrix  = np.transpose( data_2[ : , 1: ] )

theta1_min      = params[7]
theta1_max      = params[8]
theta2_min      = params[10]
theta2_max      = params[11]
density         = params[13]

theta1_array    = np.linspace(theta1_min, theta1_max, int(density))
theta2_array    = np.linspace(theta2_min, theta2_max, int(density))
theta1_matrix, theta2_matrix = np.meshgrid(theta1_array, theta2_array)
theta1_array    = theta1_matrix.flatten()
theta2_array    = theta2_matrix.flatten()
t_flip1_array   = t_flip1_matrix.flatten()
t_flip2_array   = t_flip2_matrix.flatten()


mask_flipped1 = (t_flip1_array != 0.0)
mask_flipped2 = (t_flip2_array != 0.0)


# -------------------------------------------------------------

def smooth_plot():
    plt.figure(figsize=(8, 6))
    sns.kdeplot(x=theta1_array, y=theta2_array, weights=t_flip2_array, fill=True, cmap="coolwarm", levels=10, thresh=0)

    plt.xlabel("theta1_0")
    plt.ylabel("theta2_0")
    plt.title("Heatmap of time for flip over")
    plt.show()
smooth_plot()

# -------------------------------------------------------------


def ultimate_plot():

    sample_format_dict = {
        "label"      : None,          
        "fmt"        : 'o', 
        "color"      : "to be changed",                               
        "markersize" : 1, 
        "linewidth"  : 1,
        "capsize"    : 0,
        "alpha"      : 1                                   
    }
 
    writtings = {
        "title"       : r"Titel",
        "x_ax_label"  : r"Time $t$ [$s$]",
        "y_ax_label"  : r"Angle $\theta$ [$rad$]"
    }
    all_sample_format_dicts = []    
    all_data = ()
    general_format_dict = standard_format_dict
    zoom_params         = no_zooming
    colorbar_params     = no_colorbar
    extra_label         = no_extra_label
    
    data_set_1  = np.zeros(4), None, np.zeros(4), None       

    cmap = cm.get_cmap("magma")
    norm = mcolors.Normalize(vmin=0, vmax=params[0])  
    
    for i in range(len(theta1_array)):

        sample_format_dict["color"] = cmap(norm(t_flip2_array[i])) 
        data_set                    = np.array([theta1_array[i]]), None, np.array([theta2_array[i]]), None    

        all_sample_format_dicts.append(sample_format_dict.copy())
        all_data = all_data + data_set
                    

    save_plot = True, "../plots/plot_doub_flip1_01.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)

ultimate_plot()
ultimate_plot()