
import numpy as np
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

    # Initialise dynamic objects to be drawn
line1,          = ax.plot([], [], 'o-', lw=2, color="blue")
line2,          = ax.plot([], [], 'o-', lw=2, color="green")
line3,          = ax.plot([], [], 'o-', lw=2, color="orange")
time_label      = ax.text(0.95, 0.9, " ", transform=ax.transAxes, fontsize=12, ha='right')
theta1_label    = ax.text(0.05, 0.9, " ", transform=ax.transAxes, fontsize=12, ha='left')
theta2_label    = ax.text(0.05, 0.84, " ", transform=ax.transAxes, fontsize=12, ha='left')
theta3_label    = ax.text(0.05, 0.78, " ", transform=ax.transAxes, fontsize=12, ha='left')
t_max           = 20.0
framerate_max   = 25


    # Controls the length and framerate of the animation
def framerate_control( t_end, h, t_max, framerate_max):

    multiplier = 1
    step = 1
    while(1/(h*step) > framerate_max):
        step = step + 1
    frame_count = int( min(t_end, t_max)/(h*step) )
    return frame_count, step
frame_count, step = framerate_control(t_end, h, t_max, framerate_max)


    # Update function
def update(frame):
    
        # First pendulum's position
    x1 = length_1 * np.sin(theta1[frame])
    y1 = -length_1 * np.cos(theta1[frame])
    line1.set_data([0, x1], [0, y1]) 
    
        # Second pendulum's position
    x2 = x1 + length_2 * np.sin(theta2[frame])
    y2 = y1 - length_2 * np.cos(theta2[frame])
    line2.set_data([x1, x2], [y1, y2]) 

        # Third pendulum's position
    x3 = x2 + length_3 * np.sin(theta3[frame])
    y3 = y2 - length_3 * np.cos(theta3[frame])
    line3.set_data([x2, x3], [y2, y3])  
    
        # dynamic labels
    time_label.set_text(f"$t = {time[frame]:.2f} $")
    theta1_label.set_text(f"$\\theta_1 = {theta1[frame]:.2f} $")
    theta2_label.set_text(f"$\\theta_2 = {theta2[frame]:.2f} $")
    theta3_label.set_text(f"$\\theta_3 = {theta3[frame]:.2f} $")

    return line1, line2, line3, theta1_label, theta2_label, theta3_label 

    # Create the animation
ani = FuncAnimation(fig, update, frames=range(0, frame_count*step, step), interval= h * step * 1000, blit=True)
ani.save( "../plots/anim_trip_pend.mp4")

print("Triple Pendulum animated")
