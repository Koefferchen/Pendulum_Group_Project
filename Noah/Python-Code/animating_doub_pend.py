
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

    # importing the data generated
data_bsp        = np.loadtxt("../data/data_doub_pend.txt", skiprows=1 )
time            = data_bsp[ : , 0 ]
theta1          = data_bsp[ : , 1 ]
theta1_dot      = data_bsp[ : , 2 ]
theta2          = data_bsp[ : , 3 ]
theta2_dot      = data_bsp[ : , 4 ]
params          = data_bsp[ : , 5 ]

t_end        = params[0]
h            = params[1]
g_grav       = params[2]
mass_up      = params[3]
mass_down    = params[4]
length_up    = params[5]
length_down  = params[6]


# Create the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(-1.5 * length_up, 1.5 * length_up)
ax.set_ylim(-1.5 * length_up, 1.5 * length_up)
ax.set_aspect('equal', adjustable='box')

# Pendulum line and bob
line1, = ax.plot([], [], 'o-', lw=2)
line2, = ax.plot([], [], 'o-', lw=2)


# Update function
def update(frame):
    
    theta_1 = theta1[frame]  # Angle of the first pendulum
    theta_2 = theta2[frame]  # Angle of the second pendulum
    
    # First pendulum's position
    x1 = length_up * np.sin(theta_1)
    y1 = -length_up * np.cos(theta_1)
    line1.set_data([0, x1], [0, y1])  # Line from pivot to first bob
    
    # Second pendulum's position
    x2 = x1 + length_down * np.sin(theta_2)
    y2 = y1 - length_down * np.cos(theta_2)
    line2.set_data([x1, x2], [y1, y2])  # Line from first bob to second bob
    
    return line1, line2  # Return both lines for animation


# Create the animation
ani = FuncAnimation(fig, update, frames=range(len(theta1)), interval=h * 1000, blit=True)

# Show the animation
plt.show()
