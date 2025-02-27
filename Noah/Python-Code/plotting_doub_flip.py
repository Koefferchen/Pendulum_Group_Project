
from ultimate_plotting_v7 import *
from matplotlib.pyplot import get_cmap

# ------------------------------------------------------------------------- import data

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


# ------------------------------------------------------------------------- ultimate plot approach

def ultimate_plot():

    sample_format_dict = {
        "label"      : None,          
        "fmt"        : 'o', 
        "color"      : "to be changed",                               
        "markersize" : 2, 
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

    #cmap = cm.get_cmap("viridis")
    cmap = get_cmap("viridis")
    norm = mcolors.Normalize(vmin=0, vmax=params[0])  
    
    for i in range(len(theta1_array)):

        sample_format_dict["color"] = cmap(norm(t_flip2_array[i])) 
        data_set                    = np.array([theta1_array[i]]), None, np.array([theta2_array[i]]), None    

        all_sample_format_dicts.append(sample_format_dict.copy())
        all_data = all_data + data_set
                    

    save_plot = True, "../plots/plot_doub_flip1_02.jpg"                                     
    ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict)

#ultimate_plot()
#ultimate_plot()


# ------------------------------------------------------------------------- interpolation approach


from scipy.interpolate import griddata

def landscape_plot(x_data, y_data, z_data):

        # set colormap for z-data
    cmap = get_cmap("viridis")
    norm = mcolors.Normalize(vmin=0, vmax=params[0]) 

        # Interpolate Data
    interpol_density        = 1000
    x_grid                  = np.linspace( min(x_data), max(x_data), interpol_density)
    y_grid                  = np.linspace( min(y_data), max(y_data), interpol_density)
    X_meshgrid, Y_meshgrid  = np.meshgrid(x_grid, y_grid)
    interpol_data = griddata((x_data, y_data), cmap(norm(t_flip2_array)), (X_meshgrid, Y_meshgrid), method='nearest')

    fig, ax     = plt.subplots(figsize = din_norm(4), dpi = 300)   

        # text settings 
    plt.rc ('text',   usetex    = True)
    plt.rc ('font',   family    = 'computer modern')                                     
    plt.rc ('font',   size      = 14) 

        # set title and ax-labels
    ax.set_title(r"The Flip-Over Behaviour of the Double Pendulum")                                                      
    ax.set_xlabel(r"Initial Angle $\theta_1 (0)$ [$rad$]")
    ax.set_ylabel(r"Initial Angle $\theta_2 (0)$ [$rad$]")

        # plotting
    image = ax.imshow(interpol_data, extent=(min(x_data), max(x_data), min(y_data), max(y_data)), origin='lower', cmap='viridis', aspect='auto')
    cbar = fig.colorbar(image, ax=ax, label=r"Relative Frequency of Flip-Over for $m_2$")
    plt.savefig("../plots/plot_doub_flip_new_01.jpg")

landscape_plot(theta1_array, theta2_array, t_flip2_array)
landscape_plot(theta1_array, theta2_array, t_flip2_array)