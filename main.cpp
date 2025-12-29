//
// Created by Wilson, James S on 19/12/2025.
//

#include <iostream>
#include <vector>
#include "simulation.h"

int main() {

    // Ask user for input
    setup_airfoil();

    // Initialising grid
    std::vector<STATE> grid(NX * NY);
    initialise_grid(grid);

    // Save the initial state (t=0)
    save_to_csv(grid, 0);

    // Simulation params
    double total_time = 4.0;
    double current_time = 0.0;
    int step = 0;

    // Safe dt will later implement CFL calculation
    // TODO: CFL calculation
    double dt = 0.0005;

    int save_interval = 50;   // Save data every 50 steps

    while (current_time < total_time) {

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
