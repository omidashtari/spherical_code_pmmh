""" This script was writen by Juan Cruz Gonzalez Sembla.
It will read data of the r and phi coordinates
as well as of the variable data in the equatorial plane.
A contour plot will be produced and a plot of the variable
as a function or r."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

rcParams['font.family'] = 'serif'
import matplotlib as mpl

# We write the directory and filename
directory = r'path/to/directory' # Path to directory
filename = r'Data_eq_plane.dat' # Choose filename

borders = True
normalized = False

# We load the r coordinate
Rin = 7 / 13
Rout = 20 / 13
r = np.loadtxt(directory + r'\r.dat')

# We load the phi coordinate
phi = np.loadtxt(directory + r'\phi.dat')

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

# Convert the lists to NumPy arrays
Ur = np.array(column1)
Ut = np.array(column2)
Up = np.array(column3)
T = np.array(column4)

# Reshaping
Ur = np.reshape(Ur, newshape=(np.size(r), np.size(phi)), order='C')
Ut = np.reshape(Ut, newshape=(np.size(r), np.size(phi)), order='C')
Up = np.reshape(Up, newshape=(np.size(r), np.size(phi)), order='C')
T = np.reshape(T, newshape=(np.size(r), np.size(phi)), order='C')

# Get the value of mres
mres = np.floor(2 * np.pi / phi[-1]).astype(int)
length = np.size(phi)

print(f'This is mres = {mres}')

max_index = np.unravel_index(np.argmin(Ur), Ur.shape)

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
Ur = np.concatenate([Ur, np.reshape(Ur[:, 0], (np.size(r), 1))], axis=1)
Ut = np.concatenate([Ut, np.reshape(Ut[:, 0], (np.size(r), 1))], axis=1)
Up = np.concatenate([Up, np.reshape(Up[:, 0], (np.size(r), 1))], axis=1)
T = np.concatenate([T, np.reshape(T[:, 0], (np.size(r), 1))], axis=1)

# Create a meshgrid for data plotting
r_grid, phi_grid = np.meshgrid(r, phi)

# Convert polar coordinates to Cartesian coordinates
X = r_grid * np.cos(phi_grid)
Y = r_grid * np.sin(phi_grid)

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

if borders:
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title('Contour plot of $U_r$ in the equatorial plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(r, Ur[:, max_index[1]], 'r+', lw=2)
ax.set_xlabel(r'$r$', fontsize=25)
ax.set_ylabel(r'$U_r(r)$', fontsize=25)
ax.set_title('$U_r$ as a function of $r$', fontsize=25)
ax.grid()

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

if borders:
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of $U_\theta$ in the equatorial plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(r, Ut[:, 0], 'r+', lw=2)
ax.set_xlabel(r'$r$', fontsize=25)
ax.set_ylabel(r'$U_\theta(r)$', fontsize=25)
ax.set_title(r'$U_\theta$ as a function of $r$', fontsize=25)
ax.grid()

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

if borders:
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of $U_\phi$ in the equatorial plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(r, Up[:, 0], 'r+', lw=2)
ax.set_xlabel(r'$r$', fontsize=25)
ax.set_ylabel(r'$U_\phi(r)$', fontsize=25)
ax.set_title(r'$U_\phi$ as a function of $r$', fontsize=25)
ax.grid()

""" Plot data of T """

T -= np.transpose((Rout*Rin) / r_grid - Rin)

min_val = np.nanmin(T)
max_val = np.nanmax(T)

fig, ax = plt.subplots()

if normalized:
    levels = np.linspace(-1, 1, 1000)
    T = 2 * (T - min_val) / (max_val - min_val) - 1
else:
    levels = np.linspace(min_val, max_val, 1000)

contour = plt.contourf(X, Y, np.transpose(T), levels=levels, cmap='RdBu', extend='neither')

if borders:
    plt.plot(X[:, 0], Y[:, 0], color='gray', linewidth=0.5)  # Inner border
    plt.plot(X[:, -1], Y[:, -1], color='gray', linewidth=0.5)  # Outer border

colorbar = plt.colorbar(contour)
if normalized:
    colorbar.set_ticks(np.arange(-1, 1.1, 0.2))  # Set specific colorbar ticks from -1 to 1 with 0.1 steps
    colorbar.set_ticklabels([f'{x:.1f}' for x in np.arange(-1, 1.1, 0.2)])  # Format colorbar tick labels

# Remove the plot border
for spine in plt.gca().spines.values():
    spine.set_visible(False)

plt.title(r'Contour plot of T in the equatorial plane', fontsize=20)
plt.xticks([])
plt.yticks([])

fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=10)
ax.plot(r, T[:, 0], 'r+', lw=2)
ax.set_xlabel(r'$r$', fontsize=25)
ax.set_ylabel(r'T', fontsize=25)
ax.set_title(r'T as a function of $r$', fontsize=25)
ax.grid()

plt.show()