# "pyhton3 plotting_data.py"                                        (kompilieren & ausführen)
# Tobias Neuhoff, Noah Reinhardt                                    (Autoren)
# Zum Kompilieren & Ausführen sind hier diverse packages notwendig, insbesonder latex, texlive, seaborn, matplotlib, numpy

from ultimate_plotting_v2 import *      # selbst geschriebenes Programm zum Plotten von Daten



# ------------------------- Aufgabe 1.1 --------------------------

data_SRI_1_365 = np.loadtxt("./data/SIR_1_data365.txt", skiprows=1)      

time        = data_SRI_1_365[: , 0]     # Einlesen der Daten
susceptible = data_SRI_1_365[: , 1]              
infected    = data_SRI_1_365[: , 2]
recovered   = data_SRI_1_365[: , 3]


def ultimate_plot_SIR_1_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
   
    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3
    
    writtings = {
        "titel"           : r"SIR 1 Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False   
    all_data = data_set_1 + data_set_2 + data_set_3                              
                                                                    
    save_plot = True, "./plots/SIR_1_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIR_1_365()  # Probleme bei jeweils 1. Plot mit Latex (??)
ultimate_plot_SIR_1_365()


# ------------------------- Aufgabe 1.2 --------------------------

data_SRI_2_365 = np.loadtxt("./data/SIR_2_data365.txt", skiprows=1)      

time        = data_SRI_2_365[: , 0]     # Einlesen der Daten
susceptible = data_SRI_2_365[: , 1]              
infected    = data_SRI_2_365[: , 2]
recovered   = data_SRI_2_365[: , 3]


def ultimate_plot_SIR_2_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
   
    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3
    
    writtings = {
        "titel"           : r"SIR 2 Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False   
    all_data = data_set_1 + data_set_2 + data_set_3                              
                                                                    
    save_plot = True, "./plots/SIR_2_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIR_2_365()



# ------------------------- Aufgabe 1.3 --------------------------

data_SRIV_365 = np.loadtxt("./data/SIRV_data365.txt", skiprows=1)      

time        = data_SRIV_365[: , 0]      # Einlesen der Daten
susceptible = data_SRIV_365[: , 1]              
infected    = data_SRIV_365[: , 2]
recovered   = data_SRIV_365[: , 3]
vaccinated  = data_SRIV_365[: , 4]


def ultimate_plot_SIRV_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
    
    sample_format_dict_4 = {
        "label"      : r"vaccinated",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4
    
    writtings = {
        "titel"           : r"SIRV Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False 
    data_set_4 = time, False, 100* vaccinated, False   
    all_data = data_set_1 + data_set_2 + data_set_3 + data_set_4                        
                                                                    
    save_plot = True, "./plots/SIRV_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIRV_365()



# ------------------------- Aufgabe 1.4 --------------------------

data_SRIVD_365 = np.loadtxt("./data/SIRVD_data365.txt", skiprows=1)      

time        = data_SRIVD_365[: , 0]     # Einlesen der Daten
susceptible = data_SRIVD_365[: , 1]              
infected    = data_SRIVD_365[: , 2]
recovered   = data_SRIVD_365[: , 3]
vaccinated  = data_SRIVD_365[: , 4]
dead        = data_SRIVD_365[: , 5]


def ultimate_plot_SIRVD_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
    
    sample_format_dict_4 = {
        "label"      : r"vaccinated",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_5 = {
        "label"      : r"dead",                     
        "fmt"        : '-', 
        "color"      : "black",      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4, sample_format_dict_5
    
    writtings = {
        "titel"           : r"SIRVD Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False 
    data_set_4 = time, False, 100* vaccinated, False 
    data_set_5 = time, False, 100* dead, False   
    all_data = data_set_1 + data_set_2 + data_set_3 + data_set_4 + data_set_5                        
                                                                    
    save_plot = True, "./plots/SIRVD_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIRVD_365()



# ------------------------- Aufgabe 2 ----------------------------

data_SRIVD_2_365 = np.loadtxt("./data/SIRVD_2_data365.txt", skiprows=1)      

time        = data_SRIVD_2_365[: , 0]       # Einlesen der Daten
susceptible = data_SRIVD_2_365[: , 1]              
infected    = data_SRIVD_2_365[: , 2]
recovered   = data_SRIVD_2_365[: , 3]
vaccinated  = data_SRIVD_2_365[: , 4]
dead        = data_SRIVD_2_365[: , 5]


def ultimate_plot_SIRVD_2_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
    
    sample_format_dict_4 = {
        "label"      : r"vaccinated",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_5 = {
        "label"      : r"dead",                     
        "fmt"        : '-', 
        "color"      : "black",      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4, sample_format_dict_5
    
    writtings = {
        "titel"           : r"SIRVD 2 Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False 
    data_set_4 = time, False, 100* vaccinated, False 
    data_set_5 = time, False, 100* dead, False   
    all_data = data_set_1 + data_set_2 + data_set_3 + data_set_4 + data_set_5                        
                                                                    
    save_plot = True, "./plots/SIRVD_2_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIRVD_2_365()



# ------------------------- Aufgabe 3 ----------------------------

data_SRIVDZ_365 = np.loadtxt("./data/SIRVDZ_data365.txt", skiprows=1)      

time        = data_SRIVDZ_365[: , 0]        # Einlesen der Daten
susceptible = data_SRIVDZ_365[: , 1]              
infected    = data_SRIVDZ_365[: , 2]
recovered   = data_SRIVDZ_365[: , 3]
vaccinated  = data_SRIVDZ_365[: , 4]
dead        = data_SRIVDZ_365[: , 5]
zombified   = data_SRIVDZ_365[: , 6]


def ultimate_plot_SIRVDZ_365():
    
    sample_format_dict_1 = {
        "label"      : r"susceptible",          
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[0],                              
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1                                    
    }
    
    sample_format_dict_2 = {
        "label"      : r"infected",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[3],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_3 = {
        "label"      : r"recovered",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[4],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }
    
    sample_format_dict_4 = {
        "label"      : r"vaccinated",                     
        "fmt"        : '-', 
        "color"      : sns.color_palette("dark")[9],      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_5 = {
        "label"      : r"dead",                     
        "fmt"        : '-', 
        "color"      : "black",      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    sample_format_dict_6 = {
        "label"      : r"zombified",                     
        "fmt"        : '-', 
        "color"      : "green",      
        "markersize" : 2, 
        "linewidth"  : 2,
        "capsize"    : 0,
        "alpha"      : 1
    }

    all_sample_format_dicts = sample_format_dict_1, sample_format_dict_2, sample_format_dict_3, sample_format_dict_4, sample_format_dict_5, sample_format_dict_6
    
    writtings = {
        "titel"           : r"SIRVDZ Infectioncurve",
        "x_beschriftung"  : r"Time [days]",
        "y_beschriftung"  : r"Population [\%]"
    }
    
    general_format_dict = standard_format_dict
    zoom_parameters = no_zooming
    
    data_set_1 = time, False, 100* susceptible, False       
    data_set_2 = time, False, 100* infected, False   
    data_set_3 = time, False, 100* recovered, False 
    data_set_4 = time, False, 100* vaccinated, False 
    data_set_5 = time, False, 100* dead, False
    data_set_6 = time, False, 100* zombified, False   
    all_data = data_set_1 + data_set_2 + data_set_3 + data_set_4 + data_set_5 + data_set_6                       
                                                                    
    save_plot = True, "./plots/SIRVDZ_365_graph.jpg"                                   
        
    ultimate_plot_advanced( all_data, writtings, zoom_parameters, save_plot, all_sample_format_dicts, general_format_dict )

ultimate_plot_SIRVDZ_365()
