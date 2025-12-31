//
// Created by Wilson, James S on 19/12/2025.

#ifndef AIRFOILCFD_SIMULATION_H
#define AIRFOILCFD_SIMULATION_H

#include <vector>
#include <string>

// Sim Global: Airfoil name used for directory and file naming
extern std::string GlobalAirfoilName;

// Sim Global: NACA airfoil code and Mach number to be used in testing
extern double NACA_M;
extern double NACA_P;
extern double NACA_T;
extern double MachNumber;

// Represents the fluid properties at a single grid point
struct FluidState {
    double rho;      // Density
    double rho_u;    // Momentum X
    double rho_v;    // Momentum Y
    double E;        // Total Energy
    bool is_solid;   // Boundary condition flag
    double nx;       // Surface normal X (for solid boundaries)
    double ny;       // Surface normal Y (for solid boundaries)
};

// Domain Constants
const int NX = 500;
const int NY = 200;
const double DX = 0.005;
const double DY = 0.005;
const double GAMMA = 1.4;

void promptUserForAirfoilParams();
void initializeGrid(std::vector<FluidState>& grid, double mach);
double calculateTimeStep(const std::vector<FluidState>& grid, double CFL_number);
void exportSnapshot(std::vector<FluidState>& grid, int step);
double advanceTimeStep(std::vector<FluidState>& grid, double dt);
void updateGhostCells(std::vector<FluidState>& grid);

#endif //AIRFOILCFD_SIMULATION_H