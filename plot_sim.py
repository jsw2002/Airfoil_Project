import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy.interpolate import griddata
import glob
import re

# Load NACA airfoil data files
code = input("Enter NACA code (e.g. 2412): ")
data_folder = f"data/NACA-{code}"

files = sorted(glob.glob(f"{data_folder}/output_*.csv"), 
               key=lambda x: int(re.search(r'output_(\d+).csv', x).group(1)))

# CHECK FILES FOUND
if not files:
    print("No files found!")
    exit()

# Help function to generate smoothNACA airfoil coordinates
def get_naca_coords(naca_code, n_points=200):
    # Parse code
    m = int(naca_code[0]) / 100.0
    p = int(naca_code[1]) / 10.0
    t = int(naca_code[2:]) / 100.0
    
    x = np.linspace(0, 1, n_points)
    
    # Thickness distribution
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 
                  0.2843 * x**3 - 0.1015 * x**4)
    
    # Camber line and gradients
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    
    for i, xi in enumerate(x):
        if xi < p:
            yc[i] = (m / p**2) * (2*p*xi - xi**2)
            dyc_dx[i] = (2*m / p**2) * (p - xi)
        else:
            yc[i] = (m / (1-p)**2) * ((1-2*p) + 2*p*xi - xi**2)
            dyc_dx[i] = (2*m / (1-p)**2) * (p - xi)
            
    theta = np.arctan(dyc_dx)
    
    # Upper and Lower surface
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    
    # Combine into a single polygon loop
    X_poly = np.concatenate((xu, xl[::-1]))
    Y_poly = np.concatenate((yu, yl[::-1]))
    
    return X_poly, Y_poly

# Pre-calculate the perfect shape
poly_x, poly_y = get_naca_coords(code)

# Scale and Shift to match simulation grid 
poly_x = poly_x * 1.0 + 0.5
poly_y = poly_y * 1.0 + 0.5

# Prepare grid for interpolation
df_init = pd.read_csv(files[0])
xi = np.linspace(df_init['x'].min(), df_init['x'].max(), 500)
yi = np.linspace(df_init['y'].min(), df_init['y'].max(), 200)
X, Y = np.meshgrid(xi, yi)

fig, ax = plt.subplots(figsize=(10, 5))

# Animation update function
def update(frame_idx):
    ax.clear()
    df = pd.read_csv(files[frame_idx])
    step = int(files[frame_idx].split('_')[-1].split('.')[0])

    # Interpolate Flow
    Z_rho = griddata((df['x'], df['y']), df['rho'], (X, Y), method='linear')
    
    # Plot Flow
    ax.imshow(Z_rho, extent=[0, 2.5, 0, 1.0], origin='lower', 
              cmap='jet', alpha=0.95, aspect='auto', vmin=0.5, vmax=2.5)
    
    # Plot The Perfect Wing
    ax.fill(poly_x, poly_y, 'black', zorder=10) 
    
    ax.set_title(f"NACA {code} - Step: {step}")
    return ax,

ani = animation.FuncAnimation(fig, update, frames=len(files), interval=50)
plt.show()
