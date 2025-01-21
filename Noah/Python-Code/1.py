from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.pyplot as plt
import numpy as np

def add_colorbar(fig, ax, colormap='hls', label="Colorbar", cbar_range=(-np.pi, np.pi)):
    """
    Adds a standalone colorbar strip next to the plot with the specified colormap.

    Parameters:
        fig : matplotlib.figure.Figure
            The figure object.
        ax : matplotlib.axes.Axes
            The axes object to align the colorbar with.
        colormap : str or matplotlib.colors.Colormap
            The colormap to use for the colorbar (default is 'hls').
        label : str
            Label for the colorbar.
        cbar_range : tuple
            Range of the colorbar (min, max).
    """
    # Create a new axes for the colorbar
    cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])  # [left, bottom, width, height]

    # Normalize the colorbar data
    norm = Normalize(vmin=cbar_range[0], vmax=cbar_range[1])

    # Use the specified colormap (either a string or a Colormap instance)
    cmap = plt.get_cmap(colormap) if isinstance(colormap, str) else colormap

    # Create a ScalarMappable for the colorbar
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])  # Required for the colorbar

    # Add the colorbar
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(label)  # Set the colorbar label

# Example plot data
x = np.linspace(0, 10, 100)
y = np.sin(x)

# Create the main plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(x, y, label="Sine Wave")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Main Plot with Custom Colorbar")
ax.legend()

# Add a colorbar with the 'viridis' colormap
add_colorbar(fig, ax, colormap='viridis', label="Intensity", cbar_range=(-np.pi, np.pi))

# Show the plot
plt.show()
