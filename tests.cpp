#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "simulation.h"
#include <cmath>

// Test Suite for CFD Simulation Logic
// Focuses on geometry checks, time-step calculations, and flux computations.

TEST_CASE("Airfoil Geometry Check (isPointInsideAirfoil)", "[geometry]") {
    // Manually configure global variables for NACA 2412
    NACA_M = 0.02;
    NACA_P = 0.40;
    NACA_T = 0.12;

    SECTION("Point definitely inside the airfoil") {
        // Airfoil is placed at x_LE = 0.5, centre_y = 0.5, Chord = 1.0
        // Mid-chord x = 0.5 + 0.5 = 1.0
        // Center-line y = 0.5
        // Point (1.0, 0.5) should be inside
        REQUIRE(isPointInsideAirfoil(1.0, 0.5) == true);
    }

    SECTION("Point definitely outside the airfoil") {
        // Point (1.0, 1.0) is far above the airfoil (y_local = 0.5)
        REQUIRE(isPointInsideAirfoil(1.0, 1.0) == false);
    }

    SECTION("Leading Edge (0.5, 0.5)") {
        // The leading edge is at (0.5, 0.5) in physical space
        REQUIRE(isPointInsideAirfoil(0.5, 0.5) == true);
    }

    SECTION("Trailing Edge (Near 1.5, 0.5)") {
        // The trailing edge is at (1.5, 0.5). Testing slightly inside to avoid singularity precision issues.
        REQUIRE(isPointInsideAirfoil(1.499, 0.5) == true);
    }
}

TEST_CASE("Time Step Calculation (calculateTimeStep)", "[numerics]") {
    // Setup 10x10 grid with uniform flow
    int test_NX = 10;
    int test_NY = 10;
    // Note: calculateTimeStep iterates up to NX*NY defined in header (500*200).
    // However, the function signature takes a vector.
    // Ideally we should mock NX/NY or ensure the vector is large enough.
    // Since NX/NY are const globals, we must provide a vector of size NX*NY
    // or rely on the loop bounds.
    // The loop in calculateTimeStep uses NX*NY constant.
    // We must respect that size to avoid segfaults if the loop isn't size-aware via vector.size().
    // Checking implementation: loop goes to NX*NY.
    // So we must resize to NX*NY even if we only care about values.

    std::vector<FluidState> grid(NX * NY);

    // Fill with Uniform Flow (Mach 0.8)
    double rho = 1.0;
    double u = 0.8 * std::sqrt(1.4); // Mach 0.8 * c (c=1.0 approx if p=1/1.4)
    double v = 0.0;
    double p = 1.0 / GAMMA;

    // Create state manually
    FluidState uniform;
    uniform.rho = rho;
    uniform.rho_u = rho * u;
    uniform.rho_v = rho * v;
    uniform.E = (p / (GAMMA - 1.0)) + 0.5 * rho * (u * u + v * v);
    uniform.is_solid = false;

    std::fill(grid.begin(), grid.end(), uniform);

    SECTION("Positive and reasonable dt") {
        // Verify CFL condition ensures stability
        double dt = calculateTimeStep(grid, 0.5);
        REQUIRE(dt > 0.0);
        REQUIRE(dt < 1.0);
    }
}

TEST_CASE("Inviscid Flux Computation (computeInviscidFlux)", "[physics]") {
    // Stationary Fluid State
    FluidState U;
    U.rho = 1.0;
    U.rho_u = 0.0;
    U.rho_v = 0.0;
    U.E = 2.5; // p = 1.0, E = p/(gamma-1) => 1.0/0.4 = 2.5

    double nx = 1.0;
    double ny = 0.0;

    SECTION("Zero velocity results in pressure-only flux") {
        // Access flux components
        FluidState flux = computeInviscidFlux(U, nx, ny);

        // Mass flux should be zero
        REQUIRE(flux.rho == Approx(0.0));

        // Momentum flux (x) should be Pressure * nx
        // Pressure p = (1.4 - 1.0) * (2.5 - 0) = 1.0
        REQUIRE(flux.rho_u == Approx(1.0));

        // Momentum flux (y) should be Pressure * ny (0)
        REQUIRE(flux.rho_v == Approx(0.0));

        // Energy flux should be 0 (enthalpy * vn, vn=0)
        REQUIRE(flux.E == Approx(0.0));
    }
}
