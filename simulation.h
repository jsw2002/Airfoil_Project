//
// Created by Wilson, James S on 19/12/2025.

#ifndef AIRFOILCFD_SIMULATION_H
#define AIRFOILCFD_SIMULATION_H

#include <vector>

// Defining the state structure
struct STATE {
    double rho; // Density
    double rho_u; // Momentum X
    double rho_v; // Momentum Y
    double E; // Energy
};

// Simulation constants
const int NX = 100;
const int NY = 50;
const double DX = 0.01;
const double DY = 0.01;
const double GAMMA = 1.4;

// Function to initialise grid
void initialise_grid(std::vector<STATE>& grid);

// Function to save to csv for python plotting
void save_to_csv(std::vector<STATE>& grid, int step);

// Function to update grid
void update_grid(std::vector<STATE>& grid, double dt);

#endif //AIRFOILCFD_SIMULATION_H