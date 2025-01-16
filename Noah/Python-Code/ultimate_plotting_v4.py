
# ------------------ ultimate_plotting.py ---- Version 04 ---- Last Update: 17.07.24 --------------------


import numpy as np                      # if not installed, run:   !pip install numpy
import matplotlib.pyplot as plt         # if not installed, run:   !pip install matplotlib
import seaborn as sns                   # if not installed, run:   !pip install seaborn
from scipy.optimize import curve_fit    # if not installed, run:   !pip install scipy


# ----------------------- This is the main program to create plots from data points ---------------------------
#-------------------------------- It is meant to be casted by "ultimate_plot" ---------------------------------

def ultimate_plot_advanced (all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict):
    
    gfd = general_format_dict      # eine Abkürzung 
    
    ultimate_data_check(all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict)
    
    fig, ax = plt.subplots(figsize = gfd["fig_side_lengh"], dpi = gfd["dpi_resoltion"])   # erstelle Hauptplot

    if (gfd["custom_x_range"][0] == True):                                                # falls gewünscht: beschränke den Sicht-auschnitt
        plt.xlim(gfd["custom_x_range"][1], gfd["custom_x_range"][2])
    if (gfd["custom_y_range"][0] == True):
        plt.ylim(gfd["custom_y_range"][1], gfd["custom_y_range"][2])    
    
    plt.grid(visible=True, which='major', color='black', linestyle='-', alpha=0.5)        # erstelle ein Hintergrund-Raster      
    plt.grid(visible=True, which='minor', color='black', linestyle='-', alpha=0.1)
    plt.minorticks_on()
    
    plt.rc ('text',   usetex    = True)
    plt.rc ('font',   family    = gfd["font_family"])                                     # stelle Schriftgrößen & -stile ein
    plt.rc ('font',   size      = gfd["standard_font_size"])        
    plt.rc ('axes',   titlesize = gfd["title_size"])       
    plt.rc ('axes',   labelsize = gfd["ax_label_font_size"])           
    plt.rc ('xtick',  labelsize = gfd["x_tick_font_size"])
    plt.rc ('ytick',  labelsize = gfd["y_tick_font_size"]) 
    plt.rc ('legend', fontsize  = gfd["legend_font_size"])     

    if( gfd["log_scaling_xy"][0] == True ):                                               # Stellt logarithmische Skalierung der Achsen ein (Basis 10)
        ax.set_xscale("log", base=gfd["log_scaling_xy"][2])
    if( gfd["log_scaling_xy"][1] == True ):
        ax.set_yscale("log", base=gfd["log_scaling_xy"][2])

    ax.set_title(writtings["titel"])                                                      # Beschrifte die Achsen mit Text
    ax.set_xlabel(writtings["x_beschriftung"])
    ax.set_ylabel(writtings["y_beschriftung"])
    
    
    if (zoom_parameters["do_zooming"] == True):                                           # erstelle ein Subwindow im Plot mit Zoom
        sub_axes = plt.axes(zoom_parameters["window_position"] + zoom_parameters["window_size"])
        plt.xlim(zoom_parameters["x_range"][0], zoom_parameters["x_range"][1])                                                         
        plt.ylim(zoom_parameters["y_range"][0], zoom_parameters["y_range"][1])
        plt.grid(visible=True, which='major', color='black', linestyle='-', alpha=0.5)         
        plt.grid(visible=True, which='minor', color='black', linestyle='-', alpha=0.1)
        plt.minorticks_on()

    for i in range( len(all_sample_format_dicts) ):
        
        sample_format_dict = all_sample_format_dicts[i]             # trage nacheinander alle Messdaten mit ihren Farben & Größen ein
               
        x_data     = all_data[0 + 4*i]
        x_data_err = all_data[1 + 4*i]
        y_data     = all_data[2 + 4*i]
        y_data_err = all_data[3 + 4*i]
        
        if (isinstance(x_data_err, float) == True):                 # implementiere shortcut zur Eingabe von überall gleichen Fehlern
              x_data_err = np.full( len(x_data), x_data_err )
        if (isinstance(y_data_err, float) == True):
              y_data_err = np.full( len(y_data), y_data_err )
        
        if ( x_data_err == False):                                  # damit errorbar versteht, dass keine Fehlerbalken gezeichnet werden (shortcut)       
            x_data_err = None
        if ( y_data_err == False):
            y_data_err = None

        ax.errorbar(x_data, y_data , xerr = x_data_err, yerr = y_data_err, **sample_format_dict)
         
        if (zoom_parameters["do_zooming"] == True):                # trage in das Subwindow mit Zoom alle Messdaten ein
            sub_axes.errorbar(x_data, y_data , xerr = x_data_err, yerr = y_data_err, **sample_format_dict)
            
    
    plt.legend()                                                   # erstelle dort, wo noch Platz ist, ein Label für jedes Datenset
    
    
    if (save_plot[0] == True):                                     # speicher den Plot
         plt.savefig( save_plot[1], bbox_inches='tight')
    
    
# ---------- Hilfsprogramm für Eingabeshortcuts und Überprüfung richtiger Eingabe -----------    
    
def ultimate_data_check(all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict):
    
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
        
        
    
              
# ---------- Hilfsprgramm zur Formatierung von Plots ----------------------------------------
    
def din_norm(scale):
    return [ 29.7/scale, 21.0/scale ]


# ------------ hilfreiche Standard-dictonaries : --------------------------------------------

standard_format_dict = {
    "fig_side_lengh"       : din_norm(4),           # Format des Bildes: [x_länge, y_länge]
    "dpi_resoltion"        : 300,                   # Auflösungsfaktor des Bildes
    "font_family"          : 'computer modern',     # Schriftart
    "standard_font_size"   : 11,                    # steuert die Standardtextgröße
    "title_size"           : 15,                    # Schriftgröße des Titels
    "legend_font_size"     : 10,                    # Schriftgröße der Legende
    "ax_label_font_size"   : 12,                    # Schriftgröße der x- und y-Beschriftungen
    "x_tick_font_size"     : 8,                     # Schriftgröße der x-Tick-Labels
    "y_tick_font_size"     : 8,                     # Schriftgröße der y-Tick-Labels
    "custom_x_range"       : [ False, 0, 100 ],     # Wahl des Bildausschnitts vom Koordinatensystem (Hauptplot)
    "custom_y_range"       : [ False, 0, 100 ],     # Format: [ ja/nein, unteres Limit, oberes Limit ]
    "log_scaling_xy"       : [ False, False, 10 ]   # Loarithmische Skalierung der [X-Achse, Y-Achse, Basis]
}  
# ---> als 'general_format_dict'

no_zooming = {
    "do_zooming"      : False,                   # True ---> Zoomausschnitt aktiviert
    "x_range"         : [ 0, 0 ],                # Wahl des Bildausschnitts vom Koordinatensystem (Zoom-plot)
    "y_range"         : [ 0, 0 ],
    "window_position" : [ 0, 0 ],                # relative Position des sub_windows: minimal [0, 0] maximal [1, 1]
    "window_size"     : [ 0, 0 ]                 # relative Größe des sub_windwows:   minimal [0, 0] maximal [1, 1]
}
# ---> als 'zoom_parameters'

standard_sample_dict = {
    "label"      : r"Messwerte",          # r" $ Hier ist Mathmode $ - hier nicht" ---> mit $$ einzäunen
    "fmt"        : 'o', 
    "color"      : "red",                               # sämtliche farben werdem auf englisch unterstützt
    "markersize" : 4, 
    "linewidth"  : 1,
    "capsize"    : 0,
    "alpha"      : 1  
}
# ---> Als Vorlage für ein Element der Liste 'all_sample_format_dicts'

standard_writtings = {
    "titel"           : r"Titel",
    "x_beschriftung"  : r"X-Werte [Einheit]",
    "y_beschriftung"  : r"Y-Werte [Einheit]"
}
# ---> Als Vorlage für Achsenbeschriftungen in 'writtings'



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
# ---> als Eingabeparameter für 'fmt'



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


def linear_fit(x_data, y_data, fit_range , y_err = None):
    
    def function(x, a, b):
        return a * x + b
    
    parameters, kovarianz_matrix = curve_fit(function, x_data, y_data, absolute_sigma = True, sigma = y_err)
    unsicherheiten = np.sqrt( np.diag(kovarianz_matrix) )
    
    print("---------------------------------------------------")
    print("  Best linear fit, where f(x) = a * x + b  :")
    print("  a = ", parameters[0], " +- ", unsicherheiten[0] )
    print("  b = ", parameters[1], " +- ", unsicherheiten[1] )
    print("---------------------------------------------------")
    x_fit = np.linspace(fit_range[0], fit_range[1], 300)
    y_fit = function(x_fit, parameters[0], parameters[1])
    
    return x_fit, y_fit, parameters[0], parameters[1]