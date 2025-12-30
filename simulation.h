//
// Created by Wilson, James S on 19/12/2025.

#ifndef AIRFOILCFD_SIMULATION_H
#define AIRFOILCFD_SIMULATION_H

#include <vector>
#include <string>

// Global Variable to be used
extern std::string GlobalAirfoilName;

// Defining the state structure of a cell
struct STATE {
    double rho; // Density
    double rho_u; // Momentum X
    double rho_v; // Momentum Y
    double E; // Energy
    bool is_solid; // True if inside airfoil structure
};

// Simulation constants
const int NX = 500;
const int NY = 200;
const double DX = 0.005;
const double DY = 0.005;
const double GAMMA = 1.4;

// Function to setup airfoil params from user input
void setup_airfoil();

// Function to initialise grid
void initialise_grid(std::vector<STATE>& grid);

// Function to calculate dt
double calculate_dt(const std::vector<STATE>& grid, double CFL_number);

// Function to save to csv for python plotting
void save_to_csv(std::vector<STATE>& grid, int step);

// Function to update grid
double update_grid(std::vector<STATE>& grid, double dt);

#endif //AIRFOILCFD_SIMULATION_H