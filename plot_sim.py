import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from scipy.interpolate import griddata
import glob
import re
import sys
import os

def generate_naca_coordinates(naca_code, n_points=200):
    # Parse parameter values from code
    m = int(naca_code[0]) / 100.0  # Max camber
    p = int(naca_code[1]) / 10.0   # Camber position
    t = int(naca_code[2:]) / 100.0 # Thickness

    x = np.linspace(0, 1, n_points)

    # Thickness distribution formula for NACA 4-digit series
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 
                  0.2843 * x**3 - 0.1015 * x**4)

    # Calculate camber line and gradient
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    for i, xi in enumerate(x):
        if xi < p:
            # Forward section equations
            yc[i] = (m / p**2) * (2*p*xi - xi**2)
            dyc_dx[i] = (2*m / p**2) * (p - xi)
        else:
            # Aft section equations
            yc[i] = (m / (1-p)**2) * ((1-2*p) + 2*p*xi - xi**2)
            dyc_dx[i] = (2*m / (1-p)**2) * (p - xi)

    theta = np.arctan(dyc_dx)

    # Apply thickness perpendicular to camber line
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    # Concatenate to form closed loop polygon
    X_poly = np.concatenate((xu, xl[::-1]))
    Y_poly = np.concatenate((yu, yl[::-1]))

    return X_poly, Y_poly

def main():
    # Prompt for configuration to locate data
    if len(sys.argv) > 1:
        code = sys.argv[1]
    else:
        code = input("Enter NACA code (e.g. 2412): ")
    
    data_folder = f"data/NACA-{code}"
    
    if not os.path.exists(data_folder):
        print(f"Error: Directory {data_folder} not found.")
        return

    # Sort files by numerical step index
    files = sorted(glob.glob(f"{data_folder}/output_*.csv"), 
                   key=lambda x: int(re.search(r'output_(\d+).csv', x).group(1)))

    if not files:
        print("No simulation output files found.")
        return

    # Derive airfoil geometry for visualization
    poly_x, poly_y = generate_naca_coordinates(code)
    
    # Transform geometry to match simulation domain (Chord=1.0, LE at x=0.5)
    poly_x = poly_x * 1.0 + 0.5
    poly_y = poly_y * 1.0 + 0.5

    # Setup interpolation grid based on first frame
    df_init = pd.read_csv(files[0])
    xi = np.linspace(df_init['x'].min(), df_init['x'].max(), 500)
    yi = np.linspace(df_init['y'].min(), df_init['y'].max(), 200)
    X, Y = np.meshgrid(xi, yi)

    fig, ax = plt.subplots(figsize=(10, 5))

    def update_frame(frame_idx):
        ax.clear()
        
        # Load and process snapshot
        file_path = files[frame_idx]
        df = pd.read_csv(file_path)
        step = int(re.search(r'output_(\d+).csv', file_path).group(1))

        # Griddata interpolation for contour plotting
        Z_rho = griddata((df['x'], df['y']), df['rho'], (X, Y), method='linear')

        # Visualization
        ax.imshow(Z_rho, extent=[0, 2.5, 0, 1.0], origin='lower', 
                  cmap='viridis', alpha=0.95, aspect='auto', vmin=0.5, vmax=2.5)
        
        # Overlay airfoil geometry
        ax.fill(poly_x, poly_y, 'black', zorder=10)
        
        ax.set_title(f"NACA {code} - Step: {step}")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        return ax,

    # Animate
    ani = animation.FuncAnimation(fig, update_frame, frames=len(files), interval=50)
    plt.show()

if __name__ == "__main__":
    main()
