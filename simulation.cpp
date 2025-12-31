#include "simulation.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>

// Global variable definition
std::string GlobalAirfoilName = "2412";

// Default airfoil parameters (NACA 2412)
double NACA_M = 0.02; // Max camber
double NACA_P = 0.40; // Camber position
double NACA_T = 0.12; // Thickness
double Mach_Number = 2.0;

void promptUserForAirfoilParams() {
    std::string code;
    std::cout << "Enter NACA 4-digit code: ";
    std::cin >> code;

    GlobalAirfoilName = code;

    // Default to NACA 0012 on invalid input length
    if (code.length() != 4) {
        std::cout << "Invalid NACA code length! Defaulting to NACA 0012." << std::endl;
        NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
        return;
    }

    // Validate numeric input
    for (char c : code) {
        if (!std::isdigit(c)) {
            std::cout << "NACA code contains non-digits! Defaulting to NACA 0012." << std::endl;
            NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
            return;
        }
    }

    int thick_int = std::stoi(code.substr(2, 2));

    // Enforce thickness constraints
    if (thick_int < 1 || thick_int > 40) {
        std::cout << "Thickness must be between 1 and 40! Defaulting to NACA 0012." << std::endl;
        NACA_M = 0.00; NACA_P = 0.00; NACA_T = 0.12;
        return;
    }

    // Parse digits for camber params
    double digit1 = code[0] - '0';
    double digit2 = code[1] - '0';

    NACA_M = digit1 / 100.0;
    NACA_P = digit2 / 10.0;
    NACA_T = thick_int / 100.0;
}

FluidState createFluidState(double rho, double u, double v, double p) {
    FluidState U;
    U.rho = rho;
    U.rho_u = rho * u;
    U.rho_v = rho * v;
    U.E = (p / (GAMMA - 1.0)) + 0.5 * rho * (u * u + v * v);
    return U;
}

void extractPrimitives(const FluidState& U, double& rho, double& u, double& v, double& p) {
    rho = U.rho;
    // Prevent divide-by-zero errors in near-vacuum regions
    if (std::abs(rho) < 1e-6) rho = 1e-6;

    u = U.rho_u / rho;
    v = U.rho_v / rho;
    p = (GAMMA - 1.0) * (U.E - 0.5 * rho * (u * u + v * v));
}

FluidState computeInviscidFlux(const FluidState& U, double nx, double ny) {
    FluidState flux;
    double rho, u, v, p;
    extractPrimitives(U, rho, u, v, p);

    // Calculate normal velocity component
    double vn = u * nx + v * ny;

    flux.rho   = rho * vn;
    flux.rho_u = rho * u * vn + p * nx;
    flux.rho_v = rho * v * vn + p * ny;
    flux.E     = (U.E + p) * vn;

    return flux;
}

void calculateNACA4Geometry(double zeta, double& m, double& theta, double& t_half) {
    m = 0.0;
    double dm_dx = 0.0;

    // Calculate camber line height and slope if non-symmetric
    if (NACA_M > 0.0 && NACA_P > 0.0) {
        if (zeta <= NACA_P) {
            // Forward camber equations
            m = (NACA_M / (NACA_P * NACA_P)) * (2 * NACA_P * zeta - zeta * zeta);
            dm_dx = (2 * NACA_M / (NACA_P * NACA_P)) * (NACA_P - zeta);
        } else {
            // Aft camber equations
            m = (NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (1.0 - zeta) * (1.0 + zeta - 2.0 * NACA_P);
            dm_dx = (2.0 * NACA_M / ((1 - NACA_P) * (1 - NACA_P))) * (NACA_P - zeta);
        }
        theta = std::atan(dm_dx);
    } else {
        theta = 0.0;
    }

    // Standard NACA 4-series thickness distribution
    double term1 =  0.2969 * std::sqrt(zeta);
    double term2 = -0.1260 * zeta;
    double term3 = -0.3516 * std::pow(zeta, 2);
    double term4 =  0.2843 * std::pow(zeta, 3);
    double term5 = -0.1015 * std::pow(zeta, 4);

    t_half = (NACA_T / 0.2) * (term1 + term2 + term3 + term4 + term5);
}

bool isPointInsideAirfoil(double x, double y) {
    double chord = 1.0;
    double x_LE = 0.5;
    double centre_y = 0.5;

    double x_local = (x - x_LE) / chord;
    double y_local = y - centre_y;

    // Iterative Newton-Raphson approximation for upper surface zeta
    double zeta_U = x_local;
    for (int i = 0; i < 5; ++i) {
        if (zeta_U < 0.0 || zeta_U > 1.0) break;

        double m, theta, t;
        calculateNACA4Geometry(zeta_U, m, theta, t);
        
        double x_loc = zeta_U - t * std::sin(theta);
        double error = x_local - x_loc;
        if (std::abs(error) < 1e-5) break;
        zeta_U += error;
    }

    // Iterative Newton-Raphson approximation for lower surface zeta
    double zeta_L = x_local;
    for (int i = 0; i < 5; ++i) {
        if (zeta_L < 0.0 || zeta_L > 1.0) break;

        double m, theta, t;
        calculateNACA4Geometry(zeta_L, m, theta, t);

        double x_loc = zeta_L + t * std::sin(theta);
        double error = x_local - x_loc;
        if (std::abs(error) < 1e-5) break;
        zeta_L += error;
    }

    if (zeta_U < 0.0 || zeta_U > 1.0 || zeta_L < 0.0 || zeta_L > 1.0) return false;

    // Compute bounding surfaces at corrected coordinates
    double m_U, theta_U, t_U;
    calculateNACA4Geometry(zeta_U, m_U, theta_U, t_U);
    double y_upper = m_U + t_U * std::cos(theta_U);

    double m_L, theta_L, t_L;
    calculateNACA4Geometry(zeta_L, m_L, theta_L, t_L);
    double y_lower = m_L - t_L * std::cos(theta_L);

    return (y_local <= y_upper && y_local >= y_lower);
}

void computeSurfaceNormal(double x, double y, double &nx, double &ny) {
    double chord = 1.0;
    double x_LE = 0.5;
    double x_local = (x - x_LE) / chord;

    // Clamp values to prevent singularity at leading/trailing edges
    if (x_local < 0.001) x_local = 0.001;
    if (x_local > 0.999) x_local = 0.999;

    double m, theta, t;
    calculateNACA4Geometry(x_local, m, theta, t);

    // Determine normal vector based on upper/lower surface location
    // Upper surface normal points upwards/backwards, lower points downwards/backwards
    if (y > m + 0.5) {
        nx = -std::sin(theta);
        ny =  std::cos(theta);
    } else {
        nx = -std::sin(theta);
        ny = -std::cos(theta);
    }
}

void initializeGrid(std::vector<FluidState> &grid, double mach) {
    double rho_inf = 1.0;
    double p_inf = 1.0 / GAMMA;
    double c_inf = std::sqrt(GAMMA * p_inf / rho_inf);
    double u_inf = mach * c_inf;
    double v_inf = 0.0;

    FluidState freestream = createFluidState(rho_inf, u_inf, v_inf, p_inf);
    FluidState wall = createFluidState(1.0, 0.0, 0.0, p_inf);

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = j * NX + i;
            double x_phys = (double)i * DX;
            double y_phys = (double)j * DY;

            if (isPointInsideAirfoil(x_phys, y_phys)) {
                grid[idx] = wall; 
                grid[idx].is_solid = true;
                computeSurfaceNormal(x_phys, y_phys, grid[idx].nx, grid[idx].ny);
            } else {
                grid[idx] = freestream;
                grid[idx].is_solid = false;
                grid[idx].nx = 0.0;
                grid[idx].ny = 0.0;
            }
        }
    }
}

double calculateTimeStep(const std::vector<FluidState>& grid, double CFL_number) {
    double max_signal_speed = 0.0;

    for (int idx = 0; idx < NX * NY; ++idx) {
        if (grid[idx].is_solid) continue;

        double rho = grid[idx].rho;
        double u = grid[idx].rho_u / rho;
        double v = grid[idx].rho_v / rho;
        double p = (GAMMA - 1.0) * (grid[idx].E - 0.5 * rho * (u * u + v * v));
        double c = std::sqrt(GAMMA * p / rho);

        // Calculate max wave propagation speed in both directions
        double speed_x = std::abs(u) + c;
        double speed_y = std::abs(v) + c;
        double local_speed = speed_x / DX + speed_y / DY;

        if (local_speed > max_signal_speed) {
            max_signal_speed = local_speed;
        }
    }

    return CFL_number / max_signal_speed;
}

void exportSnapshot(std::vector<FluidState> &grid, int step) {
    std::string filename = "data/NACA-" + GlobalAirfoilName + "/output_" + std::to_string(step) + ".csv";
    std::ofstream file(filename);

    file << "x,y,rho,u,v,p,is_solid\n";

    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int idx = j * NX + i;
            FluidState& S = grid[idx];

            double x = i * DX;
            double y = j * DY;

            double rho = S.rho;
            double u = 0.0, v = 0.0, p = 0.0;

            if (S.is_solid) {
                p = 1.0 / GAMMA;
            } else {
                if (rho < 1e-6) rho = 1e-6;
                u = S.rho_u / rho;
                v = S.rho_v / rho;
                p = (GAMMA - 1.0) * (S.E - 0.5 * rho * (u * u + v * v));
            }

            file << std::fixed << std::setprecision(5)
                 << x << "," << y << ","
                 << rho << "," << u << "," << v << "," << p << ","
                 << (S.is_solid ? 1 : 0) << "\n";
        }
    }
    file.close();
    std::cout << "Saved " << filename << std::endl;
}

FluidState computeRusanovFlux(FluidState& UL, FluidState& UR, double nx, double ny) {
    FluidState FL = computeInviscidFlux(UL, nx, ny);
    FluidState FR = computeInviscidFlux(UR, nx, ny);

    double rhoL, uL, vL, pL;
    extractPrimitives(UL, rhoL, uL, vL, pL);
    double rhoR, uR, vR, pR;
    extractPrimitives(UR, rhoR, uR, vR, pR);

    double cL = std::sqrt(GAMMA * pL / rhoL);
    double cR = std::sqrt(GAMMA * pR / rhoR);
    double vnL = uL * nx + vL * ny;
    double vnR = uR * nx + vR * ny;

    // Rusanov wave speed (maximum eigenvalue)
    double max_speed = std::max(std::abs(vnL) + cL, std::abs(vnR) + cR);

    // Compute flux with dissipation term
    FluidState flux;
    flux.rho   = 0.5 * (FL.rho   + FR.rho)   - 0.5 * max_speed * (UR.rho   - UL.rho);
    flux.rho_u = 0.5 * (FL.rho_u + FR.rho_u) - 0.5 * max_speed * (UR.rho_u - UL.rho_u);
    flux.rho_v = 0.5 * (FL.rho_v + FR.rho_v) - 0.5 * max_speed * (UR.rho_v - UL.rho_v);
    flux.E     = 0.5 * (FL.E     + FR.E)     - 0.5 * max_speed * (UR.E     - UL.E);

    return flux;
}

double advanceTimeStep(std::vector<FluidState> &grid, double dt) {
    std::vector<FluidState> new_grid = grid;
    double max_rho_change = 0.0;

    // Parallelize loop for performance
    #pragma omp parallel for reduction(max:max_rho_change)
    for (int j = 1; j < NY - 1; ++j) {
        for (int i = 1; i < NX - 1; ++i) {
            int idx = j * NX + i;

            if (grid[idx].is_solid) continue;

            FluidState U_Centre = grid[idx];
            FluidState NetFlux = {0, 0, 0, 0};

            int idx_L = j * NX + (i - 1);
            int idx_R = j * NX + (i + 1);
            int idx_D = (j - 1) * NX + i;
            int idx_U = (j + 1) * NX + i;

            FluidState U_L = grid[idx_L];
            FluidState U_R = grid[idx_R];
            FluidState U_D = grid[idx_D];
            FluidState U_U = grid[idx_U];

            // Compute fluxes from all four directions
            // Handle solid boundaries via pressure extrapolation
            if (U_L.is_solid) {
                double rho, u, v, p;
                extractPrimitives(U_Centre, rho, u, v, p);
                NetFlux.rho_u += p / DX; 
            } else {
                FluidState F = computeRusanovFlux(U_L, U_Centre, 1.0, 0.0);
                NetFlux.rho += F.rho / DX;
                NetFlux.rho_u += F.rho_u / DX;
                NetFlux.rho_v += F.rho_v / DX;
                NetFlux.E += F.E / DX;
            }

            if (U_R.is_solid) {
                double rho, u, v, p;
                extractPrimitives(U_Centre, rho, u, v, p);
                NetFlux.rho_u -= p / DX;
            } else {
                FluidState F = computeRusanovFlux(U_Centre, U_R, 1.0, 0.0);
                NetFlux.rho -= F.rho / DX;
                NetFlux.rho_u -= F.rho_u / DX;
                NetFlux.rho_v -= F.rho_v / DX;
                NetFlux.E -= F.E / DX;
            }

            if (U_D.is_solid) {
                double rho, u, v, p;
                extractPrimitives(U_Centre, rho, u, v, p);
                NetFlux.rho_v += p / DY;
            } else {
                FluidState G = computeRusanovFlux(U_D, U_Centre, 0.0, 1.0);
                NetFlux.rho += G.rho / DY;
                NetFlux.rho_u += G.rho_u / DY;
                NetFlux.rho_v += G.rho_v / DY;
                NetFlux.E += G.E / DY;
            }

            if (U_U.is_solid) {
                double rho, u, v, p;
                extractPrimitives(U_Centre, rho, u, v, p);
                NetFlux.rho_v -= p / DY;
            } else {
                FluidState G = computeRusanovFlux(U_Centre, U_U, 0.0, 1.0);
                NetFlux.rho -= G.rho / DY;
                NetFlux.rho_u -= G.rho_u / DY;
                NetFlux.rho_v -= G.rho_v / DY;
                NetFlux.E -= G.E / DY;
            }

            double rho_change = NetFlux.rho * dt;

            new_grid[idx].rho   += rho_change;
            new_grid[idx].rho_u += NetFlux.rho_u * dt;
            new_grid[idx].rho_v += NetFlux.rho_v * dt;
            new_grid[idx].E     += NetFlux.E * dt;

            if (std::abs(rho_change) > max_rho_change) {
                max_rho_change = std::abs(rho_change);
            }
        }
    }
    grid = new_grid;

    return max_rho_change;
}

void updateGhostCells(std::vector<FluidState>& grid) {
    // Restrict update area for efficiency
    for (int j = 55; j < 140; ++j) {
        for (int i = 95; i < 305; ++i) {
            int idx = j * NX + i;

            if (!grid[idx].is_solid) continue;

            int idx_L = j * NX + (i - 1);
            int idx_R = j * NX + (i + 1);
            int idx_D = (j - 1) * NX + i;
            int idx_U = (j + 1) * NX + i;

            FluidState* fluid_neighbour = nullptr;

            if (!grid[idx_U].is_solid) fluid_neighbour = &grid[idx_U];
            else if (!grid[idx_D].is_solid) fluid_neighbour = &grid[idx_D];
            else if (!grid[idx_L].is_solid) fluid_neighbour = &grid[idx_L];
            else if (!grid[idx_R].is_solid) fluid_neighbour = &grid[idx_R];

            if (fluid_neighbour == nullptr) continue;

            double rho, u, v, p;
            extractPrimitives(*fluid_neighbour, rho, u, v, p);

            // Reflect velocity vector across the surface normal
            double nx = grid[idx].nx;
            double ny = grid[idx].ny;
            double dot_product = u * nx + v * ny;

            double u_ghost = u - 2.0 * dot_product * nx;
            double v_ghost = v - 2.0 * dot_product * ny;

            grid[idx].rho = rho;
            grid[idx].rho_u = rho * u_ghost;
            grid[idx].rho_v = rho * v_ghost;
            grid[idx].E  = (p / (GAMMA - 1.0)) + 0.5 * rho * (u_ghost * u_ghost + v_ghost * v_ghost);
        }
    }
}