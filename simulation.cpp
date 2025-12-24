//
// Created by Wilson, James S on 19/12/2025.
//

#include "simulation.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

// Default airfoil NACA 2412
double NACA_M = 0.02; // Max camber
double NACA_P = 0.40; // Camber position
double NACA_T = 0.12; // Thickness

// Creating airfoil from input
void setup_airfoil() {
    std::string code;
    std::cout << "Enter NACA 4-digit code: ";
    std::cin >> code;

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
STATE compute_flux(const STATE& U, double nx, double ny) {
    STATE flux;
    double rho, u, v, p;
    get_state(U, rho, u, v , p);

    // Normal Velocity to face
    double vn = u * nx + v + ny;

    flux.rho   = rho * vn;
    flux.rho_u = rho * u * vn + p * nx;
    flux.rho_v = rho * v * vn + p * ny;
    flux.E     = (U.E + p) * vn;

    return flux;
}

// Boolean function to check if coordinate is inside grid
bool is_inside(double x, double y) {
    // Setting up geometry
    double chord = 1.0; // Length
    double x_LE = 0.5; // Start of airfoil in x
    double centre_y = 0.5; // Centre of airfoil in y

    // Normalise x into a percentage of chord length
    double zeta = (x - x_LE) / chord;

    // Quick exit if not in airfoil length
    if (zeta < 0.0 || zeta > chord) return false;

    // Initialising variables for airfoil structure
    double m = 0.0; // Camber height
    double dm_dx = 0.0; // Slope
    double theta = 0.0; // Angle

    // Check if airfoil is non symmetrical
    if (NACA_M > 0.0 && NACA_P > 0.0) {
        // Front of airfoil
        if (zeta <= NACA_P * chord) {
            m = (NACA_M / (NACA_P * NACA_P)) * (2 * NACA_P * zeta - zeta * zeta);
            dm_dx = (2 * NACA_M / (NACA_P * NACA_P)) * (NACA_P - zeta); // Derivative for theta
        }
        // Back of airfoil
        if (NACA_P * chord <= zeta) {
            m = (NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (1.0 - zeta) * (1.0 + zeta - 2.0 * NACA_P);
            dm_dx = (2.0 * NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (NACA_P - zeta); // Derivative for theta
        }
        // Calculate theta
        theta = std::atan(dm_dx);

        // Thickness terms
        double term1 = 0.2969 * std::sqrt(zeta);
        double term2 = -0.1260 * zeta;
        double term3 = -0.3516 * std::pow(zeta, 2);
        double term4 = 0.2843 * std::pow(zeta, 3);
        double term5 = -0.1015 * std::pow(zeta, 4);

        // Scale standard polynomial to our thickness
        double t_half = (NACA_T / 0.2) * (term1 + term2 + term3 + term4 + term5);
    }

}

// Initialising grid
void initialise_grid(std::vector<STATE> &grid) {
    // Defining normalised freestream
    double rho_inf = 1.0;
    double p_inf = 1.0 / GAMMA;
    double u_inf = 2.0; // Mach 2
    double v_inf = 0.0; // Horizontal flow only

    // Creating states for freestream and boundary conditions
    STATE freestream = create_state(rho_inf, u_inf, v_inf, p_inf);
    STATE wall = create_state(1.0, 0.0, 0.0, p_inf);

    // For loop to initialise grid
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = j * NX + i; // Convert to 1D index

            // Initialising cells
            grid[idx].rho = rho_inf;
            grid[idx].rho_u = u_inf * rho_inf;
            grid[idx].rho_v = v_inf * rho_inf;
            grid[idx].E = E_inf * rho_inf;
        }
    }
}

