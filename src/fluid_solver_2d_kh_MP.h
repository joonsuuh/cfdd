#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

float square(float x) { return x * x; }

class fluid_solver_2d {
public:
    fluid_solver_2d(float Lx0, float Ly0, int Nx0, int Ny0, float gamma0,
            float output_dt0)
            : gamma(gamma0), Lx(Lx0), Ly(Ly0), Nx(Nx0), Ny(Ny0) {
        output_dt = output_dt0;
        dx = Lx / (Nx - 2);
        dy = Ly / (Ny - 2);

        rho.resize(Nx * Ny);
        vx.resize(Nx * Ny);
        vy.resize(Nx * Ny);
        P.resize(Nx * Ny);

        mass.resize(Nx * Ny);
        mom_x.resize(Nx * Ny);
        mom_y.resize(Nx * Ny);
        energy.resize(Nx * Ny);

        rho_tmp.resize(Nx * Ny);
        vx_tmp.resize(Nx * Ny);
        vy_tmp.resize(Nx * Ny);
        P_tmp.resize(Nx * Ny);

        rho_Lx.resize(Nx * Ny);
        rho_Rx.resize(Nx * Ny);
        rho_Ly.resize(Nx * Ny);
        rho_Ry.resize(Nx * Ny);

        vx_Lx.resize(Nx * Ny);
        vx_Rx.resize(Nx * Ny);
        vx_Ly.resize(Nx * Ny);
        vx_Ry.resize(Nx * Ny);
        vy_Lx.resize(Nx * Ny);
        vy_Rx.resize(Nx * Ny);
        vy_Ly.resize(Nx * Ny);
        vy_Ry.resize(Nx * Ny);
        P_Lx.resize(Nx * Ny);
        P_Rx.resize(Nx * Ny);
        P_Ly.resize(Nx * Ny);
        P_Ry.resize(Nx * Ny);

        mass_flux_x.resize(Nx * Ny);
        mass_flux_y.resize(Nx * Ny);
        momx_flux_x.resize(Nx * Ny);
        momx_flux_y.resize(Nx * Ny);
        momy_flux_x.resize(Nx * Ny);
        momy_flux_y.resize(Nx * Ny);
        energy_flux_x.resize(Nx * Ny);
        energy_flux_y.resize(Nx * Ny);
    }

    ~fluid_solver_2d() {}

    void primitive_to_conserved() {
        // Compute conserved variables from primitive ones (use safe rho and explicit volume)
        const float vol = dx * dy;
        const float rho_eps = 1e-16;
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                int idx = i + j * Nx;

                float rho_safe = std::fmax(rho[idx], rho_eps);
                mass[idx] = rho_safe * vol;
                mom_x[idx] = mass[idx] * vx[idx];
                mom_y[idx] = mass[idx] * vy[idx];
        
                float KE = 0.5 * (square(vx[idx]) + square(vy[idx]));
                float PE = P[idx] / (gamma - 1) / rho_safe;
                energy[idx] = mass[idx] * (PE + KE);
            }
        }   
    }

    void conserved_to_primitive() {
        // Compute primitive variables from conserved ones (use safe mass -> safe rho)
        const float vol = dx * dy;
        const float mass_eps = 1e-16;
        const float P_floor = 1e-12;
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                int idx = i + j * Nx;

                float mass_safe = std::fmax(mass[idx], mass_eps);
                rho[idx] = mass_safe / vol;
                vx[idx] = mom_x[idx] / mass_safe;
                vy[idx] = mom_y[idx] / mass_safe;
                float e_spec = energy[idx] / mass_safe;
                float KE = 0.5 * (square(vx[idx]) + square(vy[idx]));
                float Pcalc = (e_spec - KE) * rho[idx] * (gamma - 1);                P[idx] = std::fmax(Pcalc, P_floor);
            }
        }   
    }

    void init(const std::vector<float> &rho0, const std::vector<float> &vx0,
            const std::vector<float> &vy0, const std::vector<float> &P0) {
        // Initialize the primitive variables using the given initial condition
        rho = rho0;
        vx = vx0;
        vy = vy0;
        P = P0; 

        // ensure ghost cells are initialized before forming conserved variables
        periodic_boundary(rho);
        periodic_boundary(vx);
        periodic_boundary(vy);
        periodic_boundary(P);

        primitive_to_conserved(); // compute conserved from primitives once
    }

    float find_dt() {
        // Find dt satisfying CFL; compute only on interior cells
        float C_cfl = 0.20;
        float min_dxdy = std::min(dx, dy);
        float max_speed = 0.0;
        const float rho_eps = 1e-16;

        // interior only: skip ghost layers
        for (int j = 1; j < Ny - 1; j++){
            for (int i = 1; i < Nx - 1; i++){
                int idx = i + j * Nx;
                float rho_safe = std::fmax(rho[idx], rho_eps);
                float sound = std::sqrt(std::fmax(0.0, gamma * P[idx] / rho_safe));
                float vel = std::sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]);
                float max_temp = sound + vel;
                if (max_temp > max_speed) max_speed = max_temp;
            }
        }
        float safe_max = std::fmax(max_speed, 1e-12);
        float dt = C_cfl * min_dxdy / safe_max;
		// std::cout << "CFL dt: " << dt << std::endl;
        return dt;
    }

    void solve(float t0, float t_end) {
        // Solve the fluid equations, starting from t0 and stoping at t_end
        float t = t0;
        int n = 0; // n labels the output file
        while (t < t_end) {
            if (t >= output_dt * n) {
                output(n);
                n += 1;
            }
            float dt = find_dt();
            step(dt);
            t += dt;
        }
    }

    void step(float dt) {
        // extrapolate a half step in time using primitive equations
        primitive_update(0.5 * dt);

        // compute fluxes
        compute_fluxes();

        // update solultion from fluxes
        update_conserved(dt);

        // update primitive variables
        conserved_to_primitive();
    }

    void periodic_boundary(std::vector<float> &f) {
        // apply periodic boundary conditions to an array f

        // f(0, j) = f(Nx - 2, j)
        // f(Nx - 1, j) = f(1 , j)
        for (int j = 0; j < Ny; j++) {
            f[j * Nx] = f[(Nx - 2) + j * Nx];
            f[(Nx - 1) + j * Nx] = f[1 + j * Nx];
        }

        // f(i, 0) = f(i, Ny - 2)
        // f(i, Ny - 1) = f(i, 1)
        for (int i = 0; i < Nx; i++) {
            f[i] = f[i + (Ny - 2) * Nx];
            f[i + (Nx * (Ny - 1))] = f[i + Nx];
        }
    }

    void primitive_update(float dt) {
        // update the primitive variables using Euler equations in primitive
        // form using an FTCS scheme; use safe rho when dividing
        const float rho_eps = 1e-16;
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++){
            for (int i = 1; i < Nx - 1; i++){
                int ij = i + j * Nx;
                int i_plus1_j = i + 1 + j * Nx;
                int i_minus1_j = i - 1 + j * Nx;
                int i_j_minus1 = i + (j - 1) * Nx;
                int i_j_plus1 = i + (j + 1) * Nx;

                float rho_c = std::fmax(rho[ij], rho_eps);

                // continuity
                float rho_inner = (vx[ij] * (rho[i_plus1_j] - rho[i_minus1_j]) / 2 / dx) + 
                        (vy[ij] * (rho[i_j_plus1] - rho[i_j_minus1]) / 2 / dy) +
                        (rho[ij] * (((vx[i_plus1_j] - vx[i_minus1_j]) / 2 / dx) + ((vy[i_j_plus1] - vy[i_j_minus1]) / 2 / dy)));
                rho_tmp[ij] = rho[ij] - 0.5 * dt * rho_inner;

                // momentum in x: advective and pressure-gradient terms
                float vx_adv = (vx[ij] * (vx[i_plus1_j] - vx[i_minus1_j]) / 2 / dx) + 
                        (vy[ij] * (vx[i_j_plus1] - vx[i_j_minus1]) / 2 / dy);
                float dPdx = (P[i_plus1_j] - P[i_minus1_j]) / 2 / dx;
                vx_tmp[ij] = vx[ij] - 0.5 * dt * (vx_adv + dPdx / rho_c);

                // momentum in y
                float vy_adv = (vx[ij] * (vy[i_plus1_j] - vy[i_minus1_j]) / 2 / dx) + 
                        (vy[ij] * (vy[i_j_plus1] - vy[i_j_minus1]) / 2 / dy);
                float dPdy = (P[i_j_plus1] - P[i_j_minus1]) / 2 / dy;
                vy_tmp[ij] = vy[ij] - 0.5 * dt * (vy_adv + dPdy / rho_c);

                // pressure evolution
                float P_inner = (vx[ij] * (P[i_plus1_j] - P[i_minus1_j]) / 2 / dx) + 
                        (vy[ij] * (P[i_j_plus1] - P[i_j_minus1]) / 2 / dy) +
                        (gamma * P[ij] * (((vx[i_plus1_j] - vx[i_minus1_j]) / 2 / dx) + ((vy[i_j_plus1] - vy[i_j_minus1]) / 2 / dy)));
                P_tmp[ij] = P[ij] - 0.5 * dt * P_inner;

            }
        }

        periodic_boundary(rho_tmp);
        periodic_boundary(vx_tmp);
        periodic_boundary(vy_tmp);
        periodic_boundary(P_tmp);

        rho = rho_tmp;
        vx = vx_tmp;
        vy = vy_tmp;
        P = P_tmp;

    }

    void extrapolate_to_interface() {
        // compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++){
            for (int i = 1; i < Nx - 1; i++){
                int ij = i + j * Nx;
                int i_plus1_j = i + 1 + j * Nx;
                int i_minus1_j = i - 1 + j * Nx;
                int i_j_minus1 = i + (j - 1) * Nx;
                int i_j_plus1 = i + (j + 1) * Nx;

                rho_Lx[ij] = rho[ij] - 0.25 * (rho[i_plus1_j] - rho[i_minus1_j]);
                rho_Rx[ij] = rho[ij] + 0.25 * (rho[i_plus1_j] - rho[i_minus1_j]);
                rho_Ly[ij] = rho[ij] - 0.25 * (rho[i_j_plus1] - rho[i_j_minus1]);
                rho_Ry[ij] = rho[ij] + 0.25 * (rho[i_j_plus1] - rho[i_j_minus1]);

                vx_Lx[ij] = vx[ij] - 0.25 * (vx[i_plus1_j] - vx[i_minus1_j]);
                vx_Rx[ij] = vx[ij] + 0.25 * (vx[i_plus1_j] - vx[i_minus1_j]);
                vx_Ly[ij] = vx[ij] - 0.25 * (vx[i_j_plus1] - vx[i_j_minus1]);
                vx_Ry[ij] = vx[ij] + 0.25 * (vx[i_j_plus1] - vx[i_j_minus1]);

                vy_Lx[ij] = vy[ij] - 0.25 * (vy[i_plus1_j] - vy[i_minus1_j]);
                vy_Rx[ij] = vy[ij] + 0.25 * (vy[i_plus1_j] - vy[i_minus1_j]);
                vy_Ly[ij] = vy[ij] - 0.25 * (vy[i_j_plus1] - vy[i_j_minus1]);
                vy_Ry[ij] = vy[ij] + 0.25 * (vy[i_j_plus1] - vy[i_j_minus1]);

                P_Lx[ij] = P[ij] - 0.25 * (P[i_plus1_j] - P[i_minus1_j]);
                P_Rx[ij] = P[ij] + 0.25 * (P[i_plus1_j] - P[i_minus1_j]);
                P_Ly[ij] = P[ij] - 0.25 * (P[i_j_plus1] - P[i_j_minus1]);
                P_Ry[ij] = P[ij] + 0.25 * (P[i_j_plus1] - P[i_j_minus1]);
            }
        }

        periodic_boundary(rho_Lx);
        periodic_boundary(rho_Rx);
        periodic_boundary(rho_Ly);
        periodic_boundary(rho_Ry);
        periodic_boundary(vx_Lx);
        periodic_boundary(vx_Rx);
        periodic_boundary(vx_Ly);
        periodic_boundary(vx_Ry);
        periodic_boundary(vy_Lx);
        periodic_boundary(vy_Rx);
        periodic_boundary(vy_Ly);
        periodic_boundary(vy_Ry);
        periodic_boundary(P_Lx);
        periodic_boundary(P_Rx);
        periodic_boundary(P_Ly);
        periodic_boundary(P_Ry);
    }

    void compute_fluxes() {
        // compute the fluxes (use safe rho when dividing / sqrt)
        float vL, vR, vL_2, vR_2, vmax_x, vmax_y;
        float F_mass_L, F_mass_R, F_momx_L, F_momx_R, F_momy_L, F_momy_R, F_energy_L, F_energy_R;
        float Q_mass_L, Q_mass_R, Q_momx_L, Q_momx_R, Q_momy_L, Q_momy_R, Q_energy_L, Q_energy_R;

        extrapolate_to_interface();
        const float rho_eps = 1e-16;
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int ij = i + j * Nx;
                int ip1j = i + 1 + j * Nx;
                int ijp1 = i + (j + 1) * Nx;
                
                // x-direction fluxes (use safe densities)
                float rhoR_x = std::fmax(rho_Rx[ij], rho_eps);
                float rhoL_xp = std::fmax(rho_Lx[ip1j], rho_eps);
                vL = sqrt(std::fmax(0.0, gamma * P_Rx[ij] / rhoR_x));
                // USE normal velocity (vx) for x-interface signal speed
                // vL_2 = std::fabs(vx_Rx[ij]);                        // changed
				vL_2 = sqrt(vx_Rx[ij]*vx_Rx[ij] + vy_Rx[ij]*vy_Rx[ij]); // changed
                vR = sqrt(std::fmax(0.0, gamma * P_Lx[ip1j] / rhoL_xp));
                // vR_2 = std::fabs(vx_Lx[ip1j]);                      // changed
				vR_2 = sqrt(vx_Lx[ip1j]*vx_Lx[ip1j] + vy_Lx[ip1j]*vy_Lx[ip1j]); // changed
                vmax_x = std::max(vL_2 + vL, vR_2 + vR);

                F_mass_L = rhoR_x * vx_Rx[ij];
                F_mass_R = rhoL_xp * vx_Lx[ip1j];
                F_momx_L = rhoR_x * vx_Rx[ij] * vx_Rx[ij] + P_Rx[ij];
                F_momx_R = rhoL_xp * vx_Lx[ip1j] * vx_Lx[ip1j] + P_Lx[ip1j];
                F_momy_L = rhoR_x * vx_Rx[ij] * vy_Rx[ij];
                F_momy_R = rhoL_xp * vx_Lx[ip1j] * vy_Lx[ip1j];
                F_energy_L = (P_Rx[ij] / (gamma - 1) + 0.5 * rhoR_x * (square(vx_Rx[ij]) + square(vy_Rx[ij])) + P_Rx[ij]) * vx_Rx[ij];
                F_energy_R = (P_Lx[ip1j] / (gamma - 1) + 0.5 * rhoL_xp * (square(vx_Lx[ip1j]) + square(vy_Lx[ip1j])) + P_Lx[ip1j]) * vx_Lx[ip1j];
                Q_mass_L = rhoR_x;
                Q_mass_R = rhoL_xp;
                Q_momx_L = rhoR_x * vx_Rx[ij];
                Q_momx_R = rhoL_xp * vx_Lx[ip1j];
                Q_momy_L = rhoR_x * vy_Rx[ij];
                Q_momy_R = rhoL_xp * vy_Lx[ip1j];
                Q_energy_L = (P_Rx[ij] / (gamma - 1) + 0.5 * rhoR_x * (square(vx_Rx[ij]) + square(vy_Rx[ij])));
                Q_energy_R = (P_Lx[ip1j] / (gamma - 1) + 0.5 * rhoL_xp * (square(vx_Lx[ip1j]) + square(vy_Lx[ip1j])));

                mass_flux_x[ij] = 0.5 * (F_mass_L + F_mass_R - vmax_x * (Q_mass_R - Q_mass_L));
                momx_flux_x[ij] = 0.5 * (F_momx_L + F_momx_R - vmax_x * (Q_momx_R - Q_momx_L));
                momy_flux_x[ij] = 0.5 * (F_momy_L + F_momy_R - vmax_x * (Q_momy_R - Q_momy_L));
                energy_flux_x[ij] = 0.5 * (F_energy_L + F_energy_R - vmax_x * (Q_energy_R - Q_energy_L));

                // y-direction fluxes (use safe densities)
                float rhoR_y = std::fmax(rho_Ry[ij], rho_eps);
                float rhoL_yp = std::fmax(rho_Ly[ijp1], rho_eps);
                vL = sqrt(std::fmax(0.0, gamma * P_Ry[ij] / rhoR_y));
                // USE normal velocity (vy) for y-interface signal speed
                // vL_2 = std::fabs(vy_Ry[ij]);                        // changed
				vL_2 = sqrt(vx_Ry[ij]*vx_Ry[ij] + vy_Ry[ij]*vy_Ry[ij]); // changed
                vR = sqrt(std::fmax(0.0, gamma * P_Ly[ijp1] / rhoL_yp));
                // vR_2 = std::fabs(vy_Ly[ijp1]);                      // changed
				vR_2 = sqrt(vx_Ly[ijp1]*vx_Ly[ijp1] + vy_Ly[ijp1]*vy_Ly[ijp1]); // changed
                vmax_y = std::max(vL_2 + vL, vR_2 + vR);

                F_mass_L = rhoR_y * vy_Ry[ij];
                F_mass_R = rhoL_yp * vy_Ly[ijp1];
                F_momx_L = rhoR_y * vy_Ry[ij] * vx_Ry[ij];
                F_momx_R = rhoL_yp * vy_Ly[ijp1] * vx_Ly[ijp1];
                F_momy_L = rhoR_y * vy_Ry[ij] * vy_Ry[ij] + P_Ry[ij];
                F_momy_R = rhoL_yp * vy_Ly[ijp1] * vy_Ly[ijp1] + P_Ly[ijp1];
                F_energy_L = (P_Ry[ij] / (gamma - 1) + 0.5 * rhoR_y * (square(vx_Ry[ij]) + square(vy_Ry[ij])) + P_Ry[ij]) * vy_Ry[ij];
                F_energy_R = (P_Ly[ijp1] / (gamma - 1) + 0.5 * rhoL_yp * (square(vx_Ly[ijp1]) + square(vy_Ly[ijp1])) + P_Ly[ijp1]) * vy_Ly[ijp1];
                Q_mass_L = rhoR_y;
                Q_mass_R = rhoL_yp;
                Q_momx_L = rhoR_y * vx_Ry[ij];
                Q_momx_R = rhoL_yp * vx_Ly[ijp1];
                Q_momy_L = rhoR_y * vy_Ry[ij];
                Q_momy_R = rhoL_yp * vy_Ly[ijp1];
                Q_energy_L = (P_Ry[ij] / (gamma - 1) + 0.5 * rhoR_y * (square(vx_Ry[ij]) + square(vy_Ry[ij])));
                Q_energy_R = (P_Ly[ijp1] / (gamma - 1) + 0.5 * rhoL_yp * (square(vx_Ly[ijp1]) + square(vy_Ly[ijp1])));

                mass_flux_y[ij] = 0.5 * (F_mass_L + F_mass_R - vmax_y * (Q_mass_R - Q_mass_L));			
                momx_flux_y[ij] = 0.5 * (F_momx_L + F_momx_R - vmax_y * (Q_momx_R - Q_momx_L));
                momy_flux_y[ij] = 0.5 * (F_momy_L + F_momy_R - vmax_y * (Q_momy_R - Q_momy_L));
                energy_flux_y[ij] = 0.5 * (F_energy_L + F_energy_R - vmax_y * (Q_energy_R - Q_energy_L));
            }
        }

        periodic_boundary(mass_flux_x);
        periodic_boundary(mass_flux_y);
        periodic_boundary(momx_flux_x);
        periodic_boundary(momx_flux_y);
        periodic_boundary(momy_flux_x);
        periodic_boundary(momy_flux_y);
        periodic_boundary(energy_flux_x);
        periodic_boundary(energy_flux_y);
    }

    void update_conserved(float dt) {
        // update the conserved variables using the fluxes
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int ij = i + j * Nx;
                int i_minus1_j = i - 1 + j * Nx;
                int i_j_minus1 = i + (j - 1) * Nx;

                mass[ij] = mass[ij] - ((mass_flux_x[ij] - mass_flux_x[i_minus1_j]) * dy * dt) - ((mass_flux_y[ij] - mass_flux_y[i_j_minus1]) * dx * dt);
                mom_x[ij] = mom_x[ij] - ((momx_flux_x[ij] - momx_flux_x[i_minus1_j]) * dy * dt) - ((momx_flux_y[ij] - momx_flux_y[i_j_minus1]) * dx * dt);
                mom_y[ij] = mom_y[ij] - ((momy_flux_x[ij] - momy_flux_x[i_minus1_j]) * dy * dt) - ((momy_flux_y[ij] - momy_flux_y[i_j_minus1]) * dx * dt);
                energy[ij] = energy[ij] - ((energy_flux_x[ij] - energy_flux_x[i_minus1_j]) * dy * dt) - ((energy_flux_y[ij] - energy_flux_y[i_j_minus1]) * dx * dt);
                
            }
        }
        // update the boundary
        periodic_boundary(mass);
        periodic_boundary(mom_x);
        periodic_boundary(mom_y);
        periodic_boundary(energy);
    }

    void output(int n) {
        std::ofstream outfile("../data/kelvin_helmholtz/output_rho_" + std::to_string(n) + ".csv");
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                int idx = i + j * Nx;
                outfile << rho[idx];
                if (i != Nx - 2)
                    outfile << ", ";
                else
                    outfile << std::endl;
            }
        }
        outfile.close();
    }

    int Nx, Ny;
    float Lx, Ly;
    float dx, dy;
    float gamma, output_dt;
    std::vector<float> rho, vx, vy, P;			// primitive variables
    std::vector<float> mass, mom_x, mom_y, energy;		// conserved variables
    // arrays to hold the results during primitive_update
    std::vector<float> rho_tmp, vx_tmp, vy_tmp, P_tmp;
    // arrays of fluxes for each conserved variable:
    std::vector<float> mass_flux_x, mass_flux_y;
    std::vector<float> momx_flux_x, momx_flux_y;
    std::vector<float> momy_flux_x, momy_flux_y;
    std::vector<float> energy_flux_x, energy_flux_y;
    // arrays for extrapolating to cell interfaces:
    std::vector<float> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
    std::vector<float> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
    std::vector<float> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
    std::vector<float> P_Lx, P_Ly, P_Rx, P_Ry;
};