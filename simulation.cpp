//
// Created by Wilson, James S on 19/12/2025.
//

#include "simulation.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>

// Define global variable
std::string GlobalAirfoilName = "2412";

// Default airfoil NACA 2412
double NACA_M = 0.02; // Max camber
double NACA_P = 0.40; // Camber position
double NACA_T = 0.12; // Thickness

// Creating airfoil from input
void setup_airfoil() {
    std::string code;
    std::cout << "Enter NACA 4-digit code: ";
    std::cin >> code;

    GlobalAirfoilName = code;

    // Safety to default to NACA 0012
    if (code.length() != 4) {
        std::cout << "Invalid NACA code length! Defaulting to NACA 0012." << std::endl;
        NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
        return;
    }

    // Safety to check characters are ints
    for (char c : code) {
        if (!std::isdigit(c)) {
            std::cout << "NACA code contains letters or symbols! Defaulting to NACA 0012." <<std::endl;
            NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
            return;
        }
    }

    // Parse thickness
    int thick_int = std::stoi(code.substr(2,2));

    // Ensure thickness >1% and <40%
    if (thick_int < 1 || thick_int > 40) {
        std::cout << "Thickness must be between 1 and 40! Defaulting to NACA 0012." <<std::endl;
        NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
        return;
    }

    // Create digits 1 and 2
    // - '0' for ASCII
    double digit1 = code[0] - '0';
    double digit2 = code[1] - '0';

    // Converting to percentage and assigning
    NACA_M = digit1 / 100.0;
    NACA_P = digit2 / 10.0;
    NACA_T = thick_int / 100.0;
}

// Create physics state
STATE create_state(double rho, double u, double v, double p) {
    STATE U;
    U.rho = rho;
    U.rho_u = rho * u;
    U.rho_v = rho * v;
    U.E = (p / (GAMMA - 1.0)) + 0.5 * rho * (u * u + v * v);
    return U;
}

// Function to get state
void get_state(const STATE& U, double& rho, double& u, double& v, double& p) {
    rho = U.rho;
    // Protect against divide-by-zero if rho is too small
    if (std::abs(rho) < 1e-6) rho = 1e-6;

    u = U.rho_u / rho;
    v = U.rho_v / rho;
    p = (GAMMA - 1.0) * (U.E - 0.5 * rho * (u * u + v * v));
}

// Function to compute the flux state
STATE compute_physical_flux(const STATE& U, double nx, double ny) {
    STATE flux;
    double rho, u, v, p;
    get_state(U, rho, u, v , p);

    // Normal Velocity to face
    double vn = u * nx + v * ny;

    flux.rho   = rho * vn;
    flux.rho_u = rho * u * vn + p * nx;
    flux.rho_v = rho * v * vn + p * ny;
    flux.E     = (U.E + p) * vn;

    return flux;
}

// Helper function to calculate NACA params at specific zeta
void get_params(double zeta, double& m, double& theta, double& t_half) {
    m = 0.0; // Camber height
    double dm_dx = 0.0; // Slope

    // Check if airfoil is non symmetrical
    if (NACA_M > 0.0 && NACA_P > 0.0) {
        // Front of airfoil
        if (zeta <= NACA_P) {
            m = (NACA_M / (NACA_P * NACA_P)) * (2 * NACA_P * zeta - zeta * zeta);
            dm_dx = (2 * NACA_M / (NACA_P * NACA_P)) * (NACA_P - zeta); // Derivative for theta
        }
        // Back of airfoil
        else {
            m = (NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (1.0 - zeta) * (1.0 + zeta - 2.0 * NACA_P);
            dm_dx = (2.0 * NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (NACA_P - zeta); // Derivative for theta
        }
        // Calculate theta
        theta = std::atan(dm_dx);
    } else {
        theta = 0.0;
    }

    // Thickness terms
    double term1 = 0.2969 * std::sqrt(zeta);
    double term2 = -0.1260 * zeta;
    double term3 = -0.3516 * std::pow(zeta, 2);
    double term4 = 0.2843 * std::pow(zeta, 3);
    double term5 = -0.1015 * std::pow(zeta, 4);

    // Scale standard polynomial to our thickness
    t_half = (NACA_T / 0.2) * (term1 + term2 + term3 + term4 + term5);

}

// Boolean function to check if coordinate is inside grid
bool is_inside(double x, double y) {
    // Setting up geometry
    double chord = 1.0; // Length
    double x_LE = 0.5; // Start of airfoil in x
    double centre_y = 0.5; // Centre of airfoil in y

    // Normalise X to grid scale
    double x_local = (x - x_LE) / chord;
    double y_local = y - centre_y;

    // Find upper zeta
    double zeta_U = x_local; // Initial guess

    // loop to guess zeta for upper
    for (int i = 0; i < 5; ++i) {
        // Quick exit if not in airfoil length
        if (zeta_U < 0.0 || zeta_U > 1.0) break;

        double m, theta, t;
        get_params(zeta_U, m, theta, t);

        // x-location of the surface at this guess
        double x_loc = zeta_U - t * std::sin(theta);

        //change zeta to fix error
        double error = x_local - x_loc;
        if (std::abs(error) < 1e-5) break;
        zeta_U += error;
    }

    // Find lower zeta
    double zeta_L = x_local; // Initial guess

    // loop to guess zeta for upper
    for (int i = 0; i < 5; ++i) {
        // Quick exit if not in airfoil length
        if (zeta_L < 0.0 || zeta_L > 1.0) break;

        double m, theta, t;
        get_params(zeta_L, m, theta, t);

        // x-location of the surface at this guess
        double x_loc = zeta_L + t * std::sin(theta);

        //change zeta to fix error
        double error = x_local - x_loc;
        if (std::abs(error) < 1e-5) break;
        zeta_L += error;
    }

    // Check if search has gone out of bounds
    if (zeta_U < 0.0 || zeta_U > 1.0 || zeta_L < 0.0 || zeta_L > 1.0) return false;

    // Get Y-heights at corrected zetas
    double m_U, theta_U, t_U;
    get_params(zeta_U, m_U, theta_U, t_U);
    double y_upper = m_U + t_U * std::cos(theta_U);

    double m_L, theta_L, t_L;
    get_params(zeta_L, m_L, theta_L, t_L);
    double y_lower = m_L - t_L * std::cos(theta_L);

    // Check if inside airfoil
    return (y_local <= y_upper && y_local >= y_lower);
}

// Helper function to get the normal vector for the airfoil
void get_normal(double x, double y, double &nx, double &ny) {
    // Normalise X
    double chord = 1.0; // Length
    double x_LE = 0.5; // Start of airfoil in x
    double x_local = (x - x_LE) / chord;

    // Safety for nose/tail
    if (x_local < 0.001) x_local = 0.001;
    if (x_local > 0.999) x_local = 0.999;

    // Get slope
    double m, theta, t;
    get_params(x_local, m, theta, t);

    // Assign m_local for Upper/Lower surface check
    double m_local = m;

    // Convert slope to normal vector
    // Slope vector is (cos(theta) sin(theta))
    // Orthogonal vector is (-sin(theta) cos(theta)
    if (y > m_local + 0.5) {
        // Upper surface
        nx = -std::sin(theta);
        ny = std::cos(theta);
    } else {
        // Lower surface
        nx = -std::sin(theta);
        ny = -std::cos(theta);
    }

}

// Initialising grid
void initialise_grid(std::vector<STATE> &grid) {
    // Defining normalised freestream
    double rho_inf = 1.0;
    double p_inf = 1.0 / GAMMA;
    double u_inf = 2.0; // Mach 2
    double v_inf = 0.0; // Horizontal flow only

    // Creating states for freestream and boundary condition
    STATE freestream = create_state(rho_inf, u_inf, v_inf, p_inf);
    STATE wall = create_state(1.0, 0.0, 0.0, p_inf);

    // For loop to initialise grid
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = j * NX + i;

            // Converting index to physical coordinate
            double x_phys = (double)i * DX;
            double y_phys = (double)j * DY;

            // Implementing airfoil check and assigning state
            if (is_inside(x_phys, y_phys)) {
                grid[idx] = wall;       // Initialise wall cell
                grid[idx].is_solid = true; // Mark for solver
                get_normal(x_phys, y_phys, grid[idx].nx, grid[idx].ny); // Call normal vector
            } else {
                grid[idx] = freestream; // Initialise freestream cell
                grid[idx].is_solid = false;
                grid[idx].nx = 0.0;
                grid[idx].ny = 0.0;
            }
        }
    }
}

// Calculate dt using the CFL condition
double calculate_dt(const std::vector<STATE>& grid, double CFL_number) {
    // Initialise variable
    double max_speed = 0.0;

    // Loop through grid
    for (int idx = 0; idx < NX * NY; ++idx) {
        // Skip inside airfoil
        if (grid[idx].is_solid) continue;

        // Get state and assign to variables
        double rho = grid[idx].rho;
        double u = grid[idx].rho_u / rho;
        double v = grid[idx].rho_v / rho;
        double p = (GAMMA - 1.0) * (grid[idx].E - 0.5 * rho * (u*u + v*v));

        // Speed of sound
        double c = std::sqrt(GAMMA * p / rho);

        // Speed of the fastest wave crossing the cell
        double speed_x = std::abs(u) + c;
        double speed_y = std::abs(v) + c;
        double local_speed = speed_x / DX + speed_y / DY;

        // Check if max speed needs updating
        if (local_speed > max_speed) {
            max_speed = local_speed;
        }
    }

    // Return dt using CFL formula: dt = CFL / (max_speed / dx)
    // Note already divided by dx/dy
    return CFL_number / max_speed;
}

void save_to_csv(std::vector<STATE> &grid, int step) {
    // Create filename
    std::string filename = "data/NACA-" + GlobalAirfoilName + "/output_" + std::to_string(step) + ".csv";
    std::ofstream file(filename);

    // Write header
    file << "x,y,rho,u,v,p,is_solid\n";

    // Loop through grid and write data
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = j * NX + i;
            STATE& S = grid[idx];

            // Convert to physical position
            double x = i * DX;
            double y = j * DY;

            // Calculate primitive variables
            double rho = S.rho;
            double u = 0.0, v = 0.0, p = 0.0;

            // Values for wall
            if (S.is_solid) {
                u = 0.0; v = 0.0;
                p = 1.0 / GAMMA;
            } else {
                // Values for freestream
                if (rho < 1e-6) rho = 1e-6; // Safety
                u = S.rho_u / rho;
                v = S.rho_v / rho;
                p = (GAMMA - 1.0) * (S.E - 0.5 * rho * (u*u + v*v));
            }
            // Write the line
            file << std::fixed << std::setprecision(5)
                 << x << "," << y << ","
                 << rho << "," << u << "," << v << "," << p << ","
                 << (S.is_solid ? 1 : 0) << "\n";
        }
    }
    file.close();
    std::cout << "Saved " << filename << std::endl;
}

// State to compute Rusanov Flux for shockwaves
STATE compute_rusanov_flux(STATE& UL, STATE& UR, double nx, double ny) {
    // Calculate physical fluxes
    STATE FL = compute_physical_flux(UL, nx, ny);
    STATE FR = compute_physical_flux(UR, nx, ny);

    // Calculate Wave Speed
    double rhoL, uL, vL, pL;
    get_state(UL, rhoL, uL, vL, pL);
    double rhoR, uR, vR, pR;
    get_state(UR, rhoR, uR, vR, pR);

    double cL = std::sqrt(GAMMA * pL / rhoL);
    double cR = std::sqrt(GAMMA * pR / rhoR);
    double vnL = uL * nx + vL * ny;
    double vnR = uR * nx + vR * ny;

    // Max speed
    double max_speed = std::max(std::abs(vnL) + cL, std::abs(vnR) + cR);

    // Create rusanov flux state
    STATE flux;
    flux.rho = 0.5 * (FL.rho   + FR.rho) - 0.5 * max_speed * (UR.rho - UL.rho);
    flux.rho_u = 0.5 * (FL.rho_u + FR.rho_u) - 0.5 * max_speed * (UR.rho_u - UL.rho_u);
    flux.rho_v = 0.5 * (FL.rho_v + FR.rho_v) - 0.5 * max_speed * (UR.rho_v - UL.rho_v);
    flux.E = 0.5 * (FL.E + FR.E) - 0.5 * max_speed * (UR.E - UL.E);

    return flux;
}

// CFD solver
double update_grid(std::vector<STATE> &grid, double dt) {
    // Temporary new grid to store updates
    std::vector<STATE> new_grid = grid;

    // Variable to track max change of density for stopping logic
    double max_rho_change = 0.0;

    // Loop to iterate over internal cells
    for (int j = 1; j < NY - 1; ++j) {
        for (int i = 1; i < NX - 1; ++i) {
            int idx = j * NX + i;

            // Check if cell is inside airfoil
            if (grid[idx].is_solid) continue;

            STATE U_Centre = grid[idx];
            STATE NetFlux = {0, 0, 0, 0};

            // Neighbouring indices
            int idx_L = j * NX + (i -1);
            int idx_R = j * NX + (i + 1);
            int idx_D = (j - 1) * NX + i;
            int idx_U = (j + 1) * NX + i;

            STATE U_L = grid[idx_L];
            STATE U_R = grid[idx_R];
            STATE U_D = grid[idx_D];
            STATE U_U = grid[idx_U];

            // Flux from left
            if (U_L.is_solid) {
                double rho, u, v, p;
                get_state(U_Centre, rho, u, v, p);
                NetFlux.rho_u += p / DX; // Pressure pushes to the right
            } else {
                STATE F = compute_rusanov_flux(U_L, U_Centre, 1.0, 0.0);
                NetFlux.rho += F.rho / DX;
                NetFlux.rho_u += F.rho_u / DX;
                NetFlux.rho_v += F.rho_v / DX;
                NetFlux.E += F.E / DX;
            }
            // Flux from right
            if (U_R.is_solid) {
                double rho, u, v, p;
                get_state(U_Centre, rho, u, v, p);
                NetFlux.rho_u -= p / DX; // Pressure pushes to the left
            } else {
                STATE F = compute_rusanov_flux(U_Centre, U_R, 1.0, 0.0);
                NetFlux.rho -= F.rho / DX;
                NetFlux.rho_u -= F.rho_u / DX;
                NetFlux.rho_v -= F.rho_v / DX;
                NetFlux.E -= F.E / DX;
            }
            // Flux from below
            if (U_D.is_solid) {
                double rho, u, v, p;
                get_state(U_Centre, rho, u, v, p);
                NetFlux.rho_v += p / DY; // Pressure pushes up
            } else {
                STATE G = compute_rusanov_flux(U_D, U_Centre, 0.0, 1.0);
                NetFlux.rho += G.rho / DY;
                NetFlux.rho_u += G.rho_u / DY;
                NetFlux.rho_v += G.rho_v / DY;
                NetFlux.E += G.E / DY;
            }
            // Flux from above
            if (U_U.is_solid) {
                double rho, u, v, p;
                get_state(U_Centre, rho, u, v, p);
                NetFlux.rho_v -= p / DY; // Pressure pushes down
            } else {
                STATE G = compute_rusanov_flux(U_Centre, U_U, 0.0, 1.0);
                NetFlux.rho -= G.rho / DY;
                NetFlux.rho_u -= G.rho_u / DY;
                NetFlux.rho_v -= G.rho_v / DY;
                NetFlux.E -= G.E / DY;
            }

            // Calculate desnity change
            double rho_change = NetFlux.rho * dt;

            // Update state of cell
            new_grid[idx].rho += rho_change;
            new_grid[idx].rho_u += NetFlux.rho_u * dt;
            new_grid[idx].rho_v += NetFlux.rho_v * dt;
            new_grid[idx].E += NetFlux.E * dt;

            // Check if density change is greater than maximum
            if (std::abs(rho_change) > max_rho_change) {
                max_rho_change = std::abs(rho_change);
            }
        }
    }
    // Overwrite grid
    grid = new_grid;

    return max_rho_change;
}

// Function to update ghost cells
void update_ghost_cells(std::vector<STATE>& grid) {
    // Loop through grid starts at
    // j = (0.3 / DY) - 10 and i = (0.5 / DX) - 10
    // ends at
    // j = (0.7 / DY) + 10 and i = (1.5 / DX) + 10
    for (int j = 110; j < 280; ++j) {
        for (int i = 190; i < 610; ++i) {
            int idx = j * NX + i;

            // Check to only update solid cells
            if (!grid[idx].is_solid) continue;

            // Neighbouring indices
            int idx_L = j * NX + (i -1);
            int idx_R = j * NX + (i + 1);
            int idx_D = (j - 1) * NX + i;
            int idx_U = (j + 1) * NX + i;

            // Initialise pointer
            STATE* fluid_neighbour = nullptr;

            // Due to geometry of airfoil first check vertically
            if (!grid[idx_U].is_solid) fluid_neighbour = &grid[idx_U];
            else if (!grid[idx_D].is_solid) fluid_neighbour = &grid[idx_D];
            else if (!grid[idx_L].is_solid) fluid_neighbour = &grid[idx_L];
            else if (!grid[idx_R].is_solid) fluid_neighbour = &grid[idx_R];

            // if cell fully inside airfoil continue
            if (fluid_neighbour == nullptr) continue;

            // Copy state of cell
            double rho, u, v, p;
            get_state(*fluid_neighbour, rho, u, v, p);

            // Reflect Velocity Vector
            double nx = grid[idx].nx;
            double ny = grid[idx].ny;

            double dot_product = u * nx + v * ny;

            double u_ghost = u - 2.0 * dot_product * nx;
            double v_ghost = v - 2.0 * dot_product * ny;

            // Write back to Ghost Cell
            grid[idx].rho = rho;
            grid[idx].rho_u = rho * u_ghost;
            grid[idx].rho_v = rho * v_ghost;
            grid[idx].E  = (p / (GAMMA - 1.0)) + 0.5 * rho * (u_ghost * u_ghost + v_ghost * v_ghost);
        }
    }
}