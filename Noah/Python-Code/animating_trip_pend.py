
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

    # importing the generated data
data_bsp        = np.loadtxt("../data/data_trip_pend.txt", skiprows=1 )
time            = data_bsp[ : , 0 ]
theta1          = data_bsp[ : , 1 ]
theta1_dot      = data_bsp[ : , 2 ]
theta2          = data_bsp[ : , 3 ]
theta2_dot      = data_bsp[ : , 4 ]
theta3          = data_bsp[ : , 5 ]
theta3_dot      = data_bsp[ : , 6 ]
params          = data_bsp[ : , 7 ]
energy          = data_bsp[ : , 8 ]

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


    # theta should always be within [-pi, +pi]
theta1   = (theta1 + np.pi) % (2*np.pi) - np.pi
theta2   = (theta2 + np.pi) % (2*np.pi) - np.pi
theta3   = (theta3 + np.pi) % (2*np.pi) - np.pi

    # Setup the the plot
fig, ax = plt.subplots()
ax.set_xlim(-1.1 * (length_1+length_2+length_3), 1.1 * (length_1+length_2+length_3) )
ax.set_ylim(-1.1 * (length_1+length_2+length_3), 1.1 * (length_1+length_2+length_3) )
ax.set_aspect('equal', adjustable='box')
ax.set_xticks([])
ax.set_yticks([])

# plt.rc ('text',   usetex    = True)
# plt.rc ('font',   family    = "computer modern")    

    # assigns each mass within [0,+inf] a bulksize within [min_size,max_size]
def bulk_size( min_size, max_size, mass):
    return min_size + (max_size-min_size)*(1 - np.exp(-mass/100))
min_size = 10
max_size = 25

    # Initialise dynamic objects to be drawn
bulk0,          = ax.plot([0,0], [0,0], 'o', markersize=10, color="black")
line1,          = ax.plot([], [],       '-', lw=3, color=sns.color_palette("dark")[0])
bulk1,          = ax.plot([], [],       'o', markersize=bulk_size(min_size, max_size, mass_1), color=sns.color_palette("dark")[0])
line2,          = ax.plot([], [],       '-', lw=3, color=sns.color_palette("dark")[9])
bulk2,          = ax.plot([], [],       'o', markersize=bulk_size(min_size, max_size, mass_2), color=sns.color_palette("dark")[9])
line3,          = ax.plot([], [],       '-', lw=3, color=sns.color_palette("dark")[6])
bulk3,          = ax.plot([], [],       'o', markersize=bulk_size(min_size, max_size, mass_3), color=sns.color_palette("dark")[6])
time_label      = ax.text(0.97, 0.97, " ", transform=ax.transAxes, fontsize=12, ha='right', va="top")
energy_label    = ax.text(0.97, 0.91, " ", transform=ax.transAxes, fontsize=12, ha='right', va="top")
theta1_label    = ax.text(0.03, 0.97, " ", transform=ax.transAxes, fontsize=12, ha='left', va="top")
theta2_label    = ax.text(0.03, 0.91, " ", transform=ax.transAxes, fontsize=12, ha='left', va="top")
theta3_label    = ax.text(0.03, 0.85, " ", transform=ax.transAxes, fontsize=12, ha='left', va="top")

    # ensure the lines are in the background
line1.set_zorder(0)
line2.set_zorder(0)
line3.set_zorder(0)

    # static parameters for the simulation
if( True ):
    g_grav_label    = ax.text(1.03, 0.97, f"$g = {g_grav:.2f}m/s^2 $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    mass_1_label    = ax.text(1.03, 0.91, f"$m_1 = {mass_1:.2f}kg $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    mass_2_label    = ax.text(1.03, 0.85, f"$m_2 = {mass_2:.2f}kg $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    mass_3_label    = ax.text(1.03, 0.79, f"$m_3 = {mass_3:.2f}kg $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    length_1_label  = ax.text(1.03, 0.73, f"$l_1 = {length_1:.2f}m $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    length_2_label  = ax.text(1.03, 0.67, f"$l_2 = {length_2:.2f}m $", transform=ax.transAxes, fontsize=12, ha='left', va="top")
    length_3_label  = ax.text(1.03, 0.61, f"$l_3 = {length_3:.2f}m $", transform=ax.transAxes, fontsize=12, ha='left', va="top")



    # Controls the length and framerate of the animation
def framerate_control( t_end, h, t_max, framerate_max):

    step = 1
    while(1/(h*step) > framerate_max):
        step = step + 1
    frame_count = int( min(t_end, t_max)/(h*step) )
    return frame_count, step

    # the animation last at most "t_max" seconds with at most "framerate_max" frames per second
t_max           = 20.0
framerate_max   = 25
frame_count, step = framerate_control(t_end, h, t_max, framerate_max)


    # update function for all dynamic objects
def update(frame):
    
        # First pendulum's position
    x1 = length_1 * np.sin(theta1[frame])
    y1 = -length_1 * np.cos(theta1[frame])
    line1.set_data([0, x1], [0, y1]) 
    bulk1.set_data([x1, x1], [y1, y1]) 
    
        # Second pendulum's position
    x2 = x1 + length_2 * np.sin(theta2[frame])
    y2 = y1 - length_2 * np.cos(theta2[frame])
    line2.set_data([x1, x2], [y1, y2])
    bulk2.set_data([x2, x2], [y2, y2]) 

        # Third pendulum's position
    x3 = x2 + length_3 * np.sin(theta3[frame])
    y3 = y2 - length_3 * np.cos(theta3[frame])
    line3.set_data([x2, x3], [y2, y3]) 
    bulk3.set_data([x3, x3], [y3, y3])  
    
        # dynamic labels
    time_label.set_text(f"$t = {time[frame]:.2f} s $")
    energy_label.set_text(f"$E = {energy[frame]:.2f} J $")
    theta1_label.set_text(f"$\\theta_1 = {theta1[frame]/np.pi:+.2f}\\pi $")
    theta2_label.set_text(f"$\\theta_2 = {theta2[frame]/np.pi:+.2f}\\pi $")
    theta3_label.set_text(f"$\\theta_3 = {theta3[frame]/np.pi:+.2f}\\pi $")

    return line1, bulk1, line2, bulk2, line3, bulk3, theta1_label, theta2_label, theta3_label, time_label, energy_label

    # Create the animation
ani = FuncAnimation(fig, update, frames=range(0, frame_count*step, step), interval= h * step * 1000, blit=True)
ani.save( "../plots/anim_trip_pend.mp4")

print("Triple Pendulum animated")
