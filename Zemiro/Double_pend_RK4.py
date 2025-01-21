import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import math
from scipy import constants

def double_pend(y,L1,L2,m1,m2,g):
    
    y[0] = theta1_dot
    y[1] = theta2_dot
    y[2]= (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2)- 2 * np.sin(theta1 - theta2) * m2 * (theta2_dot**2 * L2 + w1**2 * L1 * np.cos(theta1 - theta2))
                       )/ (L1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))
    y[3]= (2 * np.sin(theta1 - theta2)* ( theta1_dot**2 * L1 * (m1 + m2)+ g * (m1 + m2) * np.cos(theta1)+ w2**2 * L2 * m2 * np.cos(theta1 - theta2) )
                        ) / L2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)) 

    return y

def rk4(y, t, dt, n,L1,L2,m1,m2,g,derivs):
    state=L1,L2,m1,m2,g
    for i in range(n):
        k1 = dt * derivs(state,t, y[i])
        k2 = dt * derivs(state,t + dt / 2, y[i] + k1 / 2)
        k3 = dt * derivs(state,t + dt / 2, y[i] + k2 / 2)
        k4 = dt * derivs(state,t + dt, y[i] + k3)
        y_next[i] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return y_next


def RK_4M(t1,theta1_0,theta2_0,theta1_dot_0,theta2_dot_0, g, L1,L2,m1,m2, h):
    """
    ODE Solver for ODEs of the form dx/dt = v, dv/dt = -(k/l)*sin(x)
    using Runge-Kutta 4th order Method.

    Args:
        t1: Time to solve the equation up to
        x0: Initial position at t = 0
        k: Spring constant or equivalent parameter
        l: Length of pendulum or equivalent parameter
        h: Step size of the method

    Returns:
        t: Array of time points
        x: Array of position values at each time point
        v: Array of velocity values at each time point
        E: Array of energy values at each time point
    """
    # Initialize arrays
    t = np.arange(0, t1 + h, h)
    n = len(t)
    x = np.zeros(n)
    v = np.zeros(n)
    E = np.zeros(n)

    # Set initial conditions
    x[0] = x0
    v[0] = 0
    
    # Main RK4 loop
    for i in range(1, n):
        # Compute RK4 steps for position
       
    x = np.fmod(x+np.pi, 2*np.pi)-np.pi  
    return t, x, v, E
