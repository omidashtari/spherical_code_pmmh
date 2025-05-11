""" This script was writen by Juan Cruz Gonzalez Sembla.
This code will read data of the phi and theta coordinates
as well as of the variable data in the mid-gap.
A contour data plot will be produced."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rcParams
import matplotlib.cm as cm
from matplotlib.colors import LightSource
from scipy.interpolate import griddata
import os

rcParams['font.family'] = 'serif'

# We write the directory and filename
directory = r'path/to/directory' # Path to directory
filename = r'Data_mid_gap.dat' # Choose filename

# Load theta and phi coordinates
theta = np.loadtxt(directory + r'\theta.dat')
phi = np.loadtxt(directory + r'\phi.dat')

# Initialize empty lists to store the values
column1 = []
column2 = []
column3 = []
column4 = []

# Read data from file
with open(os.path.join(directory, filename), 'r') as file:
    next(file)
    for i, line in enumerate(file):
        # Split the line into two values
        values = line.split()
        # Convert the values to floats and append to the respective lists
        column1.append(float(values[0]))
        column2.append(float(values[1]))
        column3.append(float(values[2]))
        column4.append(float(values[3]))

# Convert the lists to NumPy arrays
Ur = np.array(column1)
Ut = np.array(column2)
Up = np.array(column3)
T = np.array(column4)

# Reshaping
Ur = np.reshape(Ur, newshape=(np.size(theta), np.size(phi)), order='C')
Ut = np.reshape(Ut, newshape=(np.size(theta), np.size(phi)), order='C')
Up = np.reshape(Up, newshape=(np.size(theta), np.size(phi)), order='C')
T = np.reshape(T, newshape=(np.size(theta), np.size(phi)), order='C')

# Get the value of mres
mres = np.rint(2 * np.pi / phi[-1]).astype(int)
length = np.size(phi)

i = 1
while i < mres:
    Ur = np.concatenate([Ur, Ur[:, :length]], axis=1)
    Ut = np.concatenate([Ut, Ut[:, :length]], axis=1)
    Up = np.concatenate([Up, Up[:, :length]], axis=1)
    T = np.concatenate([T, T[:, :length]], axis=1)
    i += 1

# If mres is different from 1 we have to redefine phi
if mres != 1:
    dphi = 2 * np.pi / mres / length
    phi = [dphi * k for k in range(length * mres)]

# We add the zero at the end to account for the space between the 2pi and 0
phi = np.append(phi, 0)
Ur = np.concatenate([Ur, np.reshape(Ur[:, 0], (np.size(theta), 1))], axis=1)
Ut = np.concatenate([Ut, np.reshape(Ut[:, 0], (np.size(theta), 1))], axis=1)
Up = np.concatenate([Up, np.reshape(Up[:, 0], (np.size(theta), 1))], axis=1)
T = np.concatenate([T, np.reshape(T[:, 0], (np.size(theta), 1))], axis=1)

# Create a meshgrid for data plotting
theta_grid, phi_grid = np.meshgrid(theta, phi)

# Interpolate data_mid_gap to a finer grid
theta_fine = np.linspace(theta.min(), theta.max(), 200)
phi_fine = np.linspace(phi.min(), phi.max(), 199)
phi_fine = np.append(phi_fine, 0)
theta_fine_grid, phi_fine_grid = np.meshgrid(theta_fine, phi_fine)

# Flatten grids and data for interpolation
points = np.array([theta_grid.flatten(), phi_grid.flatten()]).T

# Calculate Cartesian coordinates for the finer grid
X_fine = np.sin(theta_fine_grid) * np.cos(phi_fine_grid)
Y_fine = np.sin(theta_fine_grid) * np.sin(phi_fine_grid)
Z_fine = np.cos(theta_fine_grid)

# Interpolate data onto finer grid
Ur = np.transpose(Ur)
values = Ur.flatten()
Ur_fine = griddata(points, values, (theta_fine_grid, phi_fine_grid), method='cubic')

Ut = np.transpose(Ut)
values = Ut.flatten()
Ut_fine = griddata(points, values, (theta_fine_grid, phi_fine_grid), method='cubic')

Up = np.transpose(Up)
values = Up.flatten()
Up_fine = griddata(points, values, (theta_fine_grid, phi_fine_grid), method='cubic')

T = np.transpose(T)
values = T.flatten()
T_fine = griddata(points, values, (theta_fine_grid, phi_fine_grid), method='cubic')

""" Plot data of Ur """

# Prepare the sphere surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hide the axes and grid
ax.set_axis_off()

# Define normalization and light source
norm = plt.Normalize(vmin=-1, vmax=1)
ls = LightSource(270, 45)
rgb = ls.shade(Ur_fine, cmap=cm.RdBu, vert_exag=0.1, blend_mode='overlay')

# Plot surface with interpolated data
surf = ax.plot_surface(X_fine, Y_fine, Z_fine, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False)

# Add colorbar
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.RdBu), ax=ax, shrink=0.5, aspect=10, pad=0)

""" Plot data of Ut """

# Prepare the sphere surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hide the axes and grid
ax.set_axis_off()

# Define normalization and light source
norm = plt.Normalize(vmin=-1, vmax=1)
ls = LightSource(270, 45)
rgb = ls.shade(Ut_fine, cmap=cm.RdBu, vert_exag=0.1, blend_mode='overlay')

# Plot surface with interpolated data
surf = ax.plot_surface(X_fine, Y_fine, Z_fine, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False)

# Add colorbar
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.RdBu), ax=ax, shrink=0.5, aspect=10, pad=0)

""" Plot data of Up """

# Prepare the sphere surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hide the axes and grid
ax.set_axis_off()

# Define normalization and light source
norm = plt.Normalize(vmin=-1, vmax=1)
ls = LightSource(270, 45)
rgb = ls.shade(Up_fine, cmap=cm.RdBu, vert_exag=0.1, blend_mode='overlay')

# Plot surface with interpolated data
surf = ax.plot_surface(X_fine, Y_fine, Z_fine, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False)

# Add colorbar
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.RdBu), ax=ax, shrink=0.5, aspect=10, pad=0)

""" Plot data of T """

# Prepare the sphere surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Hide the axes and grid
ax.set_axis_off()

# Define normalization and light source
norm = plt.Normalize(vmin=-1, vmax=1)
ls = LightSource(270, 45)
rgb = ls.shade(T_fine, cmap=cm.RdBu, vert_exag=0.1, blend_mode='overlay')

# Plot surface with interpolated data
surf = ax.plot_surface(X_fine, Y_fine, Z_fine, rstride=1, cstride=1, facecolors=rgb, linewidth=0, antialiased=False)

# Add colorbar
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.RdBu), ax=ax, shrink=0.5, aspect=10, pad=0)

plt.show()
