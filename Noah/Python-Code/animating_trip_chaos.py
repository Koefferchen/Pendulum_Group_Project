
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

    # making sure an argument is passed
if len(sys.argv) != 2:
    print("Usage: python3 script.py <1 or 0>")
    sys.exit(1)


usage = int(sys.argv[1])
if( usage == 0 ):
        # plotting sensitivity on initial conditions 
    save_as = "../plots/anim_trip_chaos.gif"
    info    = "Triple Pendulum Chaos animated"
    file_a  = "../data/data_trip_chaos_a.txt"
    file_b  = "../data/data_trip_chaos_b.txt"
    file_c  = "../data/data_trip_chaos_c.txt"
elif( usage == 1 ):
        # plotting sensitivity on numerical solver
    save_as = "../plots/anim_trip_nums.gif"
    info    = "Triple Pendulum Numerical Solvers animated"
    file_a  = "../data/data_trip_RuKU2.txt"
    file_b  = "../data/data_trip_RuKu4.txt"
    file_c  = "../data/data_trip_RuKu6.txt"
else:
    print("Could not resolve argument in 'animating_triple_chaos.py'.")
    sys.exit(1)


    # importing the data generated
data_a          = np.loadtxt(file_a, skiprows=1 )
time            = data_a[ : , 0 ]
params          = data_a[ : , 7 ]

theta1a         = data_a[ : , 1 ]
theta2a         = data_a[ : , 3 ]
theta3a         = data_a[ : , 5 ]
theta1a_dot     = data_a[ : , 2 ]
theta2a_dot     = data_a[ : , 4 ]
theta3a_dot     = data_a[ : , 6 ]

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

    # Setup the the plot
fig, ax = plt.subplots()
ax.set_xlim(-1.1 * (length_1+length_2+length_3), 1.1 * (length_1+length_2+length_3) )
ax.set_ylim(-1.1 * (length_1+length_2+length_3), 1.1 * (length_1+length_2+length_3) )
ax.set_aspect('equal', adjustable='box')
ax.set_xticks([])
ax.set_yticks([])

# plt.rc ('text',   usetex    = True)
# plt.rc ('font',   family    = "computer modern")    

    # Initialise dynamic objects to be drawn
line1a,          = ax.plot([], [], 'o-', lw=2, color="blue")
line2a,          = ax.plot([], [], 'o-', lw=2, color="blue")
line3a,          = ax.plot([], [], 'o-', lw=2, color="blue")
line1b,          = ax.plot([], [], 'o-', lw=2, color="green")
line2b,          = ax.plot([], [], 'o-', lw=2, color="green")
line3b,          = ax.plot([], [], 'o-', lw=2, color="green")
line1c,          = ax.plot([], [], 'o-', lw=2, color="orange")
line2c,          = ax.plot([], [], 'o-', lw=2, color="orange")
line3c,          = ax.plot([], [], 'o-', lw=2, color="orange")
time_label      = ax.text(0.97, 0.97, " ", transform=ax.transAxes, fontsize=12, ha='right', va="top")


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

    multiplier = 1
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
    x1 = length_1 * np.sin(theta1a[frame])
    y1 = -length_1 * np.cos(theta1a[frame])
    line1a.set_data([0, x1], [0, y1]) 
    x2 = x1 + length_2 * np.sin(theta2a[frame])
    y2 = y1 - length_2 * np.cos(theta2a[frame])
    line2a.set_data([x1, x2], [y1, y2]) 
    x3 = x2 + length_3 * np.sin(theta3a[frame])
    y3 = y2 - length_3 * np.cos(theta3a[frame])
    line3a.set_data([x2, x3], [y2, y3]) 

        # Seconds pendulum's position
    x1 = length_1 * np.sin(theta1b[frame])
    y1 = -length_1 * np.cos(theta1b[frame])
    line1b.set_data([0, x1], [0, y1]) 
    x2 = x1 + length_2 * np.sin(theta2b[frame])
    y2 = y1 - length_2 * np.cos(theta2b[frame])
    line2b.set_data([x1, x2], [y1, y2]) 
    x3 = x2 + length_3 * np.sin(theta3b[frame])
    y3 = y2 - length_3 * np.cos(theta3b[frame])
    line3b.set_data([x2, x3], [y2, y3])  

        # Third pendulum's position
    x1 = length_1 * np.sin(theta1c[frame])
    y1 = -length_1 * np.cos(theta1c[frame])
    line1c.set_data([0, x1], [0, y1]) 
    x2 = x1 + length_2 * np.sin(theta2c[frame])
    y2 = y1 - length_2 * np.cos(theta2c[frame])
    line2c.set_data([x1, x2], [y1, y2]) 
    x3 = x2 + length_3 * np.sin(theta3c[frame])
    y3 = y2 - length_3 * np.cos(theta3c[frame])
    line3c.set_data([x2, x3], [y2, y3])  
    
        # dynamic labels
    time_label.set_text(f"$t = {time[frame]:.2f} s $")

    return line1a, line2a, line3a, line1b, line2b, line3b, line1c, line2c, line3c, time_label 

    # Create the animation
ani = FuncAnimation(fig, update, frames=range(0, frame_count*step, step), interval= h * step * 1000, blit=True)
ani.save( save_as )

print( info )
