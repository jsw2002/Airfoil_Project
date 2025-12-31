#include <iostream>
#include <vector>
#include <filesystem>
#include "simulation.h"

namespace fs = std::filesystem;

int main() {
    // Get simulation configuration from user
    promptUserForAirfoilParams();

    std::cout << "Enter Mach Number (e.g. 0.8, 1.5, 2.0): ";
    std::cin >> Mach_Number;

    // Prepare output directory
    std::string folder_path = "data/NACA-" + GlobalAirfoilName;
    if (fs::exists(folder_path)) {
        fs::remove_all(folder_path);
    }
    fs::create_directories(folder_path);

    // Initialize simulation domain
    std::vector<FluidState> grid(NX * NY);
    initializeGrid(grid, Mach_Number);

    // Save initial state
    exportSnapshot(grid, 0);

    // Time integration parameters
    double total_time = 4.0;
    double current_time = 0.0;
    int step = 0;
    double CFL_target = 0.5;
    double tolerance = 1e-4; // Convergence threshold
    int save_interval = 50;

    // Main time-stepping loop
    while (current_time < total_time) {
        updateGhostCells(grid);

        double dt = calculateTimeStep(grid, CFL_target);
        double max_change = advanceTimeStep(grid, dt);

        current_time += dt;
        step++;

        // Check for steady-state convergence
        if (step > 100 && max_change < tolerance) {
            exportSnapshot(grid, step);
            std::cout << "Converged at step " << step << std::endl;
            break;
        }

        // Periodic validaton output
        if (step % save_interval == 0) {
            std::cout << "Step: " << step
                      << " | Time: " << current_time << "s" << std::endl;
            exportSnapshot(grid, step);
        }
    }

    std::cout << "Simulation Complete!" << std::endl;
    return 0;
}
