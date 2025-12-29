//
// Created by Wilson, James S on 19/12/2025.
//

#include <iostream>
#include <vector>
#include <cstdlib>
#include "simulation.h"

int main() {

    // Ask user for input
    setup_airfoil();

    // Setup directory
    std::string folder = "data/NACA-" + GlobalAirfoilName;

    // If the folder exists delete it
    std::string clean_cmd = "rm -rf " + folder;
    system(clean_cmd.c_str());

    // Recreate the folder
    std::string create_cmd = "mkdir -p " + folder;
    system(create_cmd.c_str());

    // Initialising grid
    std::vector<STATE> grid(NX * NY);
    initialise_grid(grid);

    // Save the initial state (t=0)
    save_to_csv(grid, 0);

    // Simulation params
    double total_time = 4.0;
    double current_time = 0.0;
    int step = 0;
    double CFL_number = 0.5;

    int save_interval = 50;   // Save data every 50 steps

    while (current_time < total_time) {

        // Calculate dt
        double dt = calculate_dt(grid, CFL_number);

        // Run simulation
        update_grid(grid, dt);

        // Update Time
        current_time += dt;
        step++;

        // Save / Print Progress
        if (step % save_interval == 0) {
            std::cout << "Step: " << step
                      << " | Time: " << current_time << "s" << std::endl;

            save_to_csv(grid, step);
        }
    }

    std::cout << "Simulation Complete!" << std::endl;
    return 0;
}
