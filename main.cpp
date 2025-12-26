//
// Created by Wilson, James S on 19/12/2025.
//

#include <iostream>
#include <vector>
#include <string>
#include "simulation.h"

int main() {
    // -------------------------------------------------
    // 1. CONFIGURATION
    // -------------------------------------------------
    std::cout << "--- NACA Airfoil CFD Solver ---" << std::endl;

    // Ask user for the shape (Sets the global variables)
    setup_airfoil();

    // -------------------------------------------------
    // 2. INITIALIZATION
    // -------------------------------------------------
    // Create the memory for the grid (Size defined in simulation.h)
    // NX=250, NY=100 -> 25,000 cells
    std::vector<STATE> grid(NX * NY);

    std::cout << "Initializing Grid..." << std::endl;
    initialise_grid(grid);

    // Save the initial state (t=0) to verify geometry
    save_to_csv(grid, 0);

    // -------------------------------------------------
    // 3. SIMULATION LOOP
    // -------------------------------------------------
    double total_time = 4.0;  // Run for 4 simulation seconds
    double current_time = 0.0;
    int step = 0;

    // CFL Condition: dt must be small enough so info doesn't skip cells.
    // Max Speed ~ Mach 2 (u=2.0) + Sound Speed (c=1.0) = 3.0
    // dx = 0.01
    // Safe dt < 0.01 / 3.0 = 0.0033
    double dt = 0.0005;

    int save_interval = 50;   // Save data every 50 steps

    std::cout << "Starting Simulation Loop..." << std::endl;

    while (current_time < total_time) {

        // A. Run the Physics Engine
        // TODO: Uncomment this once update_grid is implemented!
        // update_grid(grid, dt);

        // B. Update Time
        current_time += dt;
        step++;

        // C. Save / Print Progress
        if (step % save_interval == 0) {
            std::cout << "Step: " << step
                      << " | Time: " << current_time << "s" << std::endl;

            save_to_csv(grid, step);
        }
    }

    std::cout << "Simulation Complete!" << std::endl;
    return 0;
}
