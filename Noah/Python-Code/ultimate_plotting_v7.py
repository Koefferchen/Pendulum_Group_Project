
# ------------------ ultimate_plotting_v7.py ---- Version 07 ---- Last Update: 07.02.25 -----------------------



    # if some <library> not installed: "!pip install <library>"
import numpy as np                      
import matplotlib.pyplot as plt         
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import seaborn as sns                   
from scipy.optimize import curve_fit    



    # generalised main program for plotting -- needs to by called by ultimate_plot
def ultimate_plot_advanced (all_data, writtings, zoom_params, colorbar_params, extra_label, save_plot, all_sample_format_dicts, general_format_dict):
    
        # abbreviation
    gfd = general_format_dict      
    
        # check if data is given in correct form
    ultimate_data_check(all_data, writtings, zoom_params, save_plot, all_sample_format_dicts, general_format_dict)
    
        # initialize plot of given size & resolution
    fig, ax     = plt.subplots(figsize = gfd["fig_side_lengh"], dpi = gfd["dpi_resoltion"])   
    init_legend = False

        # add grid lines to the background
    plt.grid(visible=True, which='major', color='black', linestyle='-', alpha=0.5)           
    plt.grid(visible=True, which='minor', color='black', linestyle='-', alpha=0.1)
    plt.minorticks_on()

        # limit plot window in x & y direction
    if (gfd["custom_x_range"][0] ):                                               
        plt.xlim(gfd["custom_x_range"][1], gfd["custom_x_range"][2])
    if (gfd["custom_y_range"][0] ):
        plt.ylim(gfd["custom_y_range"][1], gfd["custom_y_range"][2])    
    
        # realise font & label settings
    plt.rc ('text',   usetex    = True)
    plt.rc ('font',   family    = gfd["font_family"])                                     
    plt.rc ('font',   size      = gfd["standard_font_size"])        
    plt.rc ('axes',   titlesize = gfd["title_size"])       
    plt.rc ('axes',   labelsize = gfd["ax_label_font_size"])           
    plt.rc ('xtick',  labelsize = gfd["x_tick_font_size"])
    plt.rc ('ytick',  labelsize = gfd["y_tick_font_size"]) 
    plt.rc ('legend', fontsize  = gfd["legend_font_size"])     

        # scale axes logarithmic 
    if( gfd["log_scaling_xy"][0] ):                                              
        ax.set_xscale("log", base=gfd["log_scaling_xy"][2])
    if( gfd["log_scaling_xy"][1] ):
        ax.set_yscale("log", base=gfd["log_scaling_xy"][2])

        # set title and ax-labels
    ax.set_title(writtings["title"])                                                      
    ax.set_xlabel(writtings["x_ax_label"])
    ax.set_ylabel(writtings["y_ax_label"])
    
        # add the given extra label for additional information to the plot
    if(extra_label["do_label"]):
        ax.text(extra_label["position"][0], extra_label["position"][1], extra_label["content"], transform=ax.transAxes, fontsize=extra_label["font_size"], ha='left', va="top")

        # add the given colorbar to the side of the plot
    if (colorbar_params["do_cbar"]):
        add_colorbar( fig, colorbar_params )
    
        # create a window within the plot which shows a magnified section of the plot
    if (zoom_params["do_zooming"] == True):                                          
        sub_axes = plt.axes(zoom_params["window_position"] + zoom_params["window_size"])
        plt.xlim(zoom_params["x_range"][0], zoom_params["x_range"][1])                                                         
        plt.ylim(zoom_params["y_range"][0], zoom_params["y_range"][1])
        plt.grid(visible=True, which='major', color='black', linestyle='-', alpha=0.5)         
        plt.grid(visible=True, which='minor', color='black', linestyle='-', alpha=0.1)
        plt.minorticks_on()

        # for each tuple of (data_set, style_instructions)
    for i in range( len(all_sample_format_dicts) ):
        
            # unpack style_instructions
        sample_format_dict = all_sample_format_dicts[i]             

            # unpack data_set
        x_data     = all_data[0 + 4*i]
        x_data_err = all_data[1 + 4*i]
        y_data     = all_data[2 + 4*i]
        y_data_err = all_data[3 + 4*i]
        
            # if a float "err" is given instead of an error_array: expand the error to be (err, err, err, ...)
        if (isinstance(x_data_err, float) == True):                 
              x_data_err = np.full( len(x_data), x_data_err )
        if (isinstance(y_data_err, float) == True):
              y_data_err = np.full( len(y_data), y_data_err )
        
            # check if any labels were given
        if not( sample_format_dict["label"] == None ):
            init_legend = True

            # plot the data_set
        ax.errorbar(x_data, y_data , xerr = x_data_err, yerr = y_data_err, **sample_format_dict)
         
            # plot the zoom_window
        if (zoom_params["do_zooming"] == True):                 
            sub_axes.errorbar(x_data, y_data , xerr = x_data_err, yerr = y_data_err, **sample_format_dict)

        # plots the labels    
    if(init_legend):
        plt.legend()                                                
    
        # save the plot at given place in given format
    if (save_plot[0] == True):                                      
         plt.savefig( save_plot[1], bbox_inches='tight')
    
    
# ----------------------------------- Support Programs -----------------------------------
    

    # check if data is entered correctly
def ultimate_data_check(all_data, writtings, zoom_params, save_plot, all_sample_format_dicts, general_format_dict):
    
    if (len(all_data) != 4 * len(all_sample_format_dicts)):
        print("Die Anzahl der Argumente ist kein Vielfaches von 4!")
        print("Prüfe, ob pro Datenset alles vorhanden ist:  x_values, x_error, y_values, y_error .")
        print("Anzahl Daten-Arrays: ", len(all_data))
        print("Prüfe, ob die für jedes Datenset ein formater übergeben wird.")
        print("Anzahl Formater: ", len(all_sample_format_dicts))
    
    for i in range( len(all_sample_format_dicts) ):
              
        x     = all_data[4*i]
        x_err = all_data[4*i + 1]
        y     = all_data[4*i + 2]
        y_err = all_data[4*i + 3]  
        
        if not (isinstance(x, np.ndarray) == True and isinstance(y, np.ndarray) == True ):
            print("Alle X- und Y-Werte müssen in Numpy-arrays gegeben werden.")
            print("Probiere es mit ' X = np.array(X) '. ")


    # format plot according to the DIN-norm    
def din_norm(scale):
    return [ 29.7/scale, 21.0/scale ]


    # adds colorbar at the side of the plot
def add_colorbar(fig, colorbar_params):
    cbar_ax = fig.add_axes(colorbar_params["position"] + colorbar_params["size"]) 
    sm = ScalarMappable(Normalize(vmin=colorbar_params["scale_range"][0], vmax=colorbar_params["scale_range"][1]), cmap=colorbar_params["colormap"])  # Create ScalarMappable
    sm.set_array([]) 
    fig.colorbar(sm, cax=cbar_ax).set_label(colorbar_params["title"]) 



# --------------------- Standard dictionaries to be used in "ultimate_plot" ---------------------


    # option for "general_format_dict"
standard_format_dict = {
    "fig_side_lengh"       : din_norm(4),           # Format des Bildes: [x_länge, y_länge]
    "dpi_resoltion"        : 300,                   # Auflösungsfaktor des Bildes
    "font_family"          : 'computer modern',     # Schriftart
    "standard_font_size"   : 11,                    # steuert die Standardtextgröße
    "title_size"           : 15,                    # Schriftgröße des titles
    "legend_font_size"     : 10,                    # Schriftgröße der Legende
    "ax_label_font_size"   : 12,                    # Schriftgröße der x- und y-Beschriftungen
    "x_tick_font_size"     : 8,                     # Schriftgröße der x-Tick-Labels
    "y_tick_font_size"     : 8,                     # Schriftgröße der y-Tick-Labels
    "custom_x_range"       : [ False, 0, 100 ],     # Wahl des Bildausschnitts vom Koordinatensystem (Hauptplot)
    "custom_y_range"       : [ False, 0, 100 ],     # Format: [ ja/nein, unteres Limit, oberes Limit ]
    "log_scaling_xy"       : [ False, False, 10 ]   # Loarithmische Skalierung der [X-Achse, Y-Achse, Basis]
}  

    # option for "zoom_params"
no_zooming = {
    "do_zooming"      : False,                      # True ---> Zoomausschnitt aktiviert
    "x_range"         : [ 0, 0 ],                   # Wahl des Bildausschnitts vom Koordinatensystem (Zoom-plot)
    "y_range"         : [ 0, 0 ],
    "window_position" : [ 0, 0 ],                   # relative Position des sub_windows: minimal [0, 0] maximal [1, 1]
    "window_size"     : [ 0, 0 ]                    # relative Größe des sub_windwows:   minimal [0, 0] maximal [1, 1]
}

    # option for "colorbar_params"
no_colorbar = {
    "do_cbar"       : False,
    "position"      : [0.92, 0.15],
    "size"          : [0.03, 0.70],
    "scale_range"   : [-np.pi, +np.pi],
    "title"         : "Colorbar",
    "colormap"      : "magma"
}

    # option for "extra_label"
no_extra_label = {
    "do_label"  :   False,
    "position"  :   [1.03, 0.97],
    "font_size" :   12,
    "content"   :   "I used these parameters: \n m = 1kg"
}

    # general form of each "sample_format_dict"
standard_sample_dict = {
    "label"      : r"Measurement 1",          
    "fmt"        : 'o', 
    "color"      : "red",                               
    "markersize" : 4, 
    "linewidth"  : 1,
    "capsize"    : 0,
    "alpha"      : 1  
}

    # option for "writtings"
standard_writtings = {
    "title"           : r"title",
    "x_ax_label"  : r"X-Werte [Einheit]",
    "y_ax_label"  : r"Y-Werte [Einheit]"
}


#--------------------- Wahl der Datenpunkt-Form ------------------------------------
# 
#    -	Durchgezogene Linie
#    --	Gestrichelte Linie
#    -.	Abwechselnd gestrichelte und gepunktete Linie
#    :	Gepunktete Linie
#    o	Einzelne Punkte, Darstellung als farbige Kreise
#    s	Einzelne Punkte, Darstellung als farbige Rechtecke
#    D	Einzelne Punkte, Darstellung als Diamant-Form
#    ^	Einzelne Punkte, Darstellung als farbige Dreiecke
#    x	Einzelne Punkte, Darstellung als farbige x-Zeichen
#    *  Einzelne Punkte, Darstellung als farbige *-Zeichen
#    +	Einzelne Punkte, Darstellung als farbige +-Zeichen"
#
# ---> option for "fmt"



# -------------- Grundlegender Aufbau eines beliebigen Funktionen-Fits ---------------------  

def least_sqare_fit():
    
    def function(x, a, b):
        return a * ( 1 - np.exp(-x/b) )
    
    initial_guess = [10, 80]
    x_data = np.linspace(0,100,5)
    y_data = x_data**2
    y_err  = None
    
    parameters, kovarianz_matrix = curve_fit(function, x_data, y_data, absolute_sigma = True, p0 = initial_guess, sigma = y_err)
    unsicherheiten = np.sqrt( np.diag(kovarianz_matrix) )
    
    print("---------------------------------------------------")
    print("  a = ", parameters[0], " +- ", unsicherheiten[0] )
    print("  b = ", parameters[1], " +- ", unsicherheiten[1] )
    print("---------------------------------------------------")
    
    x_fit = np.linspace(0, 100, 300)
    y_fit = function(x_fit, parameters[0], parameters[1])
    
    return x_fit, y_fit


def linear_fit(x_data, y_data, fit_range=False , y_err = None, info=False):
    
    def function(x, a, b):
        return a * x + b
    
    if not(fit_range):
        fit_range = [ x_data[0], x_data[::-1][0] ]

    parameters, kovarianz_matrix = curve_fit(function, x_data, y_data, absolute_sigma = True, sigma = y_err)
    unsicherheiten = np.sqrt( np.diag(kovarianz_matrix) )
    if(info):
        print("---------------------------------------------------")
        print("  Best linear fit, where f(x) = a * x + b  :")
        print("  a = ", parameters[0], " +- ", unsicherheiten[0] )
        print("  b = ", parameters[1], " +- ", unsicherheiten[1] )
        print("---------------------------------------------------")
    x_fit = np.linspace(fit_range[0], fit_range[1], 2)
    y_fit = function(x_fit, parameters[0], parameters[1])
    
    return x_fit, y_fit, parameters[0], parameters[1]



