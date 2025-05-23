""" This script was writen by Juan Cruz Gonzalez Sembla.
This code will read data of the r and theta coordinates
as well as of the variable data in the meridional plane.
A contour plot will be produced."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
import matplotlib as mpl
import os

# We write the directory and filename
directory = r'path/to/directory' # Path to directory
filename = r'Data_mer_plane.dat' # Choose filename

borders = True
normalized = False

# We load the r coordinate
Rin = 7 / 13
Rout = 20 / 13
r = np.loadtxt(directory + r'\r.dat')

# We load the phi coordinate
theta = np.loadtxt(directory + r'\theta.dat')

# Initialize empty lists to store the values
column1 = []
column2 = []
column3 = []
column4 = []

with open(os.path.join(directory, filename), 'r') as file:
    next(file)
    # Loop through each line in the file
    for line in file:
        # Split the line into two values
        values = line.split()
        # Convert the values to floats and append to the respective lists
        column1.append(float(values[0]))
        column2.append(float(values[1]))
        column3.append(float(values[2]))
        column4.append(float(values[3]))

        # temp = np.append(np.array(row_data), np.array(row_data)[0])
        # data_merid_plane[i, :] = np.append(temp, temp[-1])
        # i += 1

# Convert the lists to NumPy arrays
Ur = np.array(column1)
Ut = np.array(column2)
Up = np.array(column3)
T = np.array(column4)

# Reshaping
Ur = np.reshape(Ur, newshape=(np.size(r), np.size(theta)), order='C')
Ut = np.reshape(Ut, newshape=(np.size(r), np.size(theta)), order='C')
Up = np.reshape(Up, newshape=(np.size(r), np.size(theta)), order='C')
T = np.reshape(T, newshape=(np.size(r), np.size(theta)), order='C')

# We add the 0 and pi values
theta = np.append(0, theta)
theta = np.append(theta, np.pi)
Ur = np.concatenate([np.reshape(Ur[:, 0], (np.size(r), 1)), Ur, np.reshape(Ur[:, -1], (np.size(r), 1))], axis=1)
Ut = np.concatenate([np.reshape(Ut[:, 0], (np.size(r), 1)), Ut, np.reshape(Ut[:, -1], (np.size(r), 1))], axis=1)
Up = np.concatenate([np.reshape(Up[:, 0], (np.size(r), 1)), Up, np.reshape(Up[:, -1], (np.size(r), 1))], axis=1)
T = np.concatenate([np.reshape(T[:, 0], (np.size(r), 1)), T, np.reshape(T[:, -1], (np.size(r), 1))], axis=1)

# Create a meshgrid for data plotting
r_grid, theta_grid = np.meshgrid(r, theta)

# Convert polar coordinates to Cartesian coordinates
X = r_grid * np.sin(theta_grid)
Y = r_grid * np.cos(theta_grid)

# Now we plot

mpl.rc('text', usetex=True)

""" Plot data of Ur """

min_val = np.nanmin(Ur)
max_val = np.nanmax(Ur)

fig, ax = plt.subplots()

if normalized:
    levels = np.linspace(-1, 1, 1000)
    Ur = 2 * (Ur - min_val) / (max_val - min_val) - 1
else:
    levels = np.linspace(min_val, max_val, 1000)

contour = plt.contourf(X, Y, np.transpose(Ur), levels=levels, cmap='RdBu', extend='neither')

# Add the gray borderline based on the data grid edge
if borders:
    plt.plot(X[0, :], Y[0, :], color='gray', linewidth=0.5)  # Straight top border
    plt.plot(X[-1, :], Y[-1, :], color='gray', linewidth=0.5)  # Straight bottom border
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

plt.axis('equal')

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of $U_r$ in the meridional plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(theta, Ur[np.size(r) // 2, :], 'r+', lw=2, label="mid gap")
ax.plot(theta, Ur[np.size(r) // 3, :], 'b+', lw=2, label="1/3 gap")
ax.plot(theta, Ur[np.size(r) * 2 // 3, :], 'g+', lw=2, label="2/3 gap")
ax.set_xlabel(r'$\theta$', fontsize=25)
ax.set_ylabel(r'$U_r(\theta)$', fontsize=25)
ax.set_title(r'$U_r$ as a function of $\theta$', fontsize=25)
ax.grid()
plt.legend()

""" Plot data of Ut """

min_val = np.nanmin(Ut)
max_val = np.nanmax(Ut)

fig, ax = plt.subplots()

if normalized:
    levels = np.linspace(-1, 1, 1000)
    Ut = 2 * (Ut - min_val) / (max_val - min_val) - 1
else:
    levels = np.linspace(min_val, max_val, 1000)

contour = plt.contourf(X, Y, np.transpose(Ut), levels=levels, cmap='RdBu', extend='neither')

# Add the gray borderline based on the data grid edge
if borders:
    plt.plot(X[0, :], Y[0, :], color='gray', linewidth=0.5)  # Straight top border
    plt.plot(X[-1, :], Y[-1, :], color='gray', linewidth=0.5)  # Straight bottom border
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

plt.axis('equal')

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of $U_\theta$ in the meridional plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(theta, Ut[np.size(r) // 2, :], 'r+', lw=2, label="mid gap")
ax.plot(theta, Ut[np.size(r) // 3, :], 'b+', lw=2, label="1/3 gap")
ax.plot(theta, Ut[np.size(r) * 2 // 3, :], 'g+', lw=2, label="2/3 gap")
ax.set_xlabel(r'$\theta$', fontsize=25)
ax.set_ylabel(r'$U_\theta(\theta)$', fontsize=25)
ax.set_title(r'$U_\theta$ as a function of $\theta$', fontsize=25)
ax.grid()
plt.legend()

""" Plot data of Up """

min_val = np.nanmin(Up)
max_val = np.nanmax(Up)

fig, ax = plt.subplots()

if normalized:
    levels = np.linspace(-1, 1, 1000)
    Up = 2 * (Up - min_val) / (max_val - min_val) - 1
else:
    levels = np.linspace(min_val, max_val, 1000)

contour = plt.contourf(X, Y, np.transpose(Up), levels=levels, cmap='RdBu', extend='neither')

# Add the gray borderline based on the data grid edge
if borders:
    plt.plot(X[0, :], Y[0, :], color='gray', linewidth=0.5)  # Straight top border
    plt.plot(X[-1, :], Y[-1, :], color='gray', linewidth=0.5)  # Straight bottom border
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

plt.axis('equal')

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of $U_\phi$ in the meridional plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(theta, Up[np.size(r) // 2, :], 'r+', lw=2, label="mid gap")
ax.plot(theta, Up[np.size(r) // 3, :], 'b+', lw=2, label="1/3 gap")
ax.plot(theta, Up[np.size(r) * 2 // 3, :], 'g+', lw=2, label="2/3 gap")
ax.set_xlabel(r'$\theta$', fontsize=25)
ax.set_ylabel(r'$U_\phi(\theta)$', fontsize=25)
ax.set_title(r'$U_\phi$ as a function of $\theta$', fontsize=25)
ax.grid()
plt.legend()

""" Plot data of T """

min_val = np.nanmin(T)
max_val = np.nanmax(T)

fig, ax = plt.subplots()

if normalized:
    levels = np.linspace(-1, 1, 1000)
    T = 2 * (T - min_val) / (max_val - min_val) - 1
else:
    levels = np.linspace(min_val, max_val, 1000)

contour = plt.contourf(X, Y, np.transpose(T), levels=levels, cmap='RdBu', extend='neither')

# Add the gray borderline based on the data grid edge
if borders:
    plt.plot(X[0, :], Y[0, :], color='gray', linewidth=0.5)  # Straight top border
    plt.plot(X[-1, :], Y[-1, :], color='gray', linewidth=0.5)  # Straight bottom border
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

plt.axis('equal')

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of T in the meridional plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(theta, T[np.size(r) // 2, :], 'r+', lw=2, label="mid gap")
ax.plot(theta, T[np.size(r) // 3, :], 'b+', lw=2, label="1/3 gap")
ax.plot(theta, T[np.size(r) * 2 // 3, :], 'g+', lw=2, label="2/3 gap")
ax.set_xlabel(r'$\theta$', fontsize=25)
ax.set_ylabel(r'$T(\theta)$', fontsize=25)
ax.set_title(r'T as a function of $\theta$', fontsize=25)
ax.grid()
plt.legend()

plt.show()


