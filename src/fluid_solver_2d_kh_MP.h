#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

double square(double x) { return x * x; }

class fluid_solver_2d {
public:
	fluid_solver_2d(double Lx0, double Ly0, int Nx0, int Ny0, double gamma0,
			double output_dt0)
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
		// DONE: Compute conserved variables from primitive ones
		#pragma omp parallel for collapse(2)
		for (int j = 0; j < Ny; j++){
			for (int i = 0; i < Nx; i++){
				int idx = i + j * Nx;

				mass[idx] = rho[idx] * dx * dy;
				mom_x[idx] = mass[idx] * vx[idx];
				mom_y[idx] = mass[idx] * vy[idx];
		
				double KE = 0.5 * (square(vx[idx]) + square(vy[idx]));
				double PE = P[idx] / (gamma - 1) / rho[idx];
				energy[idx] = mass[idx] * (PE + KE);
			}
		}   
	}

	void conserved_to_primitive() {
		// DONE: Compute primitive variables from conserved ones
		#pragma omp parallel for collapse(2)
		for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                int idx = i + j * Nx;

                rho[idx] = mass[idx] / dx / dy;
                vx[idx] = mom_x[idx] / mass[idx];
                vy[idx] = mom_y[idx] / mass[idx];
                P[idx] = (energy[idx] / mass[idx] - (0.5 * (square(vx[idx]) + square(vy[idx])))) * rho[idx] * (gamma - 1);
            }
        }   
	}

	void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
			const std::vector<double> &vy0, const std::vector<double> &P0) {
		// DONE??: Initialize the primitive variables using the given initial 
		// condition
		rho = rho0;
		vx = vx0;
		vy = vy0;
		P = P0; 

		primitive_to_conserved(); // i think we need to do this once before time step
	}

	double find_dt() {
		// DONE??: Find the optimal dt that satisfies the CFL condition, and return
		// its value
		double C_cfl = 0.4; // change to find stability limit
		double min = std::min(dx, dy);
		double max_temp;
		double max = 0.0;
		for (int j = 0; j < Ny; j++){  // interior points only
			for (int i = 0; i < Nx; i++){
				int idx = i + j * Nx;
				max_temp = std::sqrt(gamma * P[idx] / rho[idx]) + std::sqrt(vx[idx] * vx[idx] + vy[idx] + vy[idx]);
				if (max_temp > max) {
					max = max_temp;
				}
			}
		}
		// periodic_boundary(max_vec);
		// max_vec.clear();

		return C_cfl * min / max;
	}

	void solve(double t0, double t_end) {
		// Solve the fluid equations, starting from t0 and stoping at t_end
		double t = t0;
		int n = 0; // n labels the output file
		while (t < t_end) {
			if (t >= output_dt * n) {
				output(n);
				n += 1;
			}
			double dt = find_dt();
			step(dt);
			t += dt;
		}
	}

	void step(double dt) {
		// extrapolate a half step in time using primitive equations
		primitive_update(0.5 * dt);

		// compute fluxes
		compute_fluxes();

		// update solultion from fluxes
		update_conserved(dt);

		// update primitive variables
		conserved_to_primitive();
	}

	void periodic_boundary(std::vector<double> &f) {
		// DONE: apply periodic boundary conditions to an array f

		// f(0, j) = f(Nx - 2, j)
		// 		=> idx = i + j * Nx = 0 + j * Nx
		// f(Nx - 1, j) = f(1 , j)
		//		=> idx = (Nx - 1) + j * Nx
		for (int j = 0; j < Ny; j++) {
			f[j * Nx] = f[(Nx - 2) + j * Nx];
			f[(Nx - 1) + j * Nx] = f[1 + j * Nx];
		}

		// f(i, 0) = f(i, Ny - 2)
		// f(i, Ny - 1) = f(i, 1)
		for (int i = 0; i < Nx; i++) {
			f[i] = f[i + (Ny - 2) * Nx];
			f[i + (Ny - 1) * Nx] = f[i + Nx];
		}
	}

	void primitive_update(double dt) {
		// DONE: update the primitive variables using Euler equations in primitive
		// form using an FTCS scheme
		#pragma omp parallel for collapse(2)
		for (int j = 1; j < Ny - 1; j++){
            for (int i = 1; i < Nx - 1; i++){
                int ij = i + j * Nx;
				int i_plus1_j = i + 1 + j * Nx;
				int i_minus1_j = i - 1 + j * Nx;
				int i_j_minus1 = i + (j - 1) * Nx;
				int i_j_plus1 = i + (j + 1) * Nx;

				double rho_inner = (vx[ij] * (rho[i_plus1_j] - rho[i_minus1_j]) / 2 / dx) + 
								   (vy[ij] * (rho[i_j_plus1] - rho[i_j_minus1]) / 2 / dy) +
								   (rho[ij] * (((vx[i_plus1_j] - vx[i_minus1_j]) / 2 / dx) + ((vy[i_j_plus1] - vy[i_j_minus1]) / 2 / dy)));
				rho_tmp[ij] = rho[ij] - 0.5 * dt * rho_inner;

				double vx_inner = (vx[ij] * (vx[i_plus1_j] - vx[i_minus1_j]) / 2 / dx) + 
								   (vy[ij] * (vx[i_j_plus1] - vx[i_j_minus1]) / 2 / dy) +
								   ((P[i_plus1_j] - P[i_minus1_j]) / 2 / dx);
				vx_tmp[ij] = vx[ij] - 0.5 * dt * vx_inner / rho[ij];

				double vy_inner = (vx[ij] * (vy[i_plus1_j] - vy[i_minus1_j]) / 2 / dx) + 
								   (vy[ij] * (vy[i_j_plus1] - vy[i_j_minus1]) / 2 / dy) +
								   ((P[i_j_plus1] - P[i_j_minus1]) / 2 / dy);
				vy_tmp[ij] = vy[ij] - 0.5 * dt * vy_inner / rho[ij];

				double P_inner = (vx[ij] * (P[i_plus1_j] - P[i_minus1_j]) / 2 / dx) + 
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
		// DONE: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
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
		// DONE: compute the fluxes
		double vL, vR, vL_2, vR_2, vmax_x, vmax_y;
		double F_mass_L, F_mass_R, F_momx_L, F_momx_R, F_momy_L, F_momy_R, F_energy_L, F_energy_R;
		double Q_mass_L, Q_mass_R, Q_momx_L, Q_momx_R, Q_momy_L, Q_momy_R, Q_energy_L, Q_energy_R;

		extrapolate_to_interface();
		#pragma omp parallel for collapse(2)
		for (int j = 1; j < Ny - 1; j++) {
			for (int i = 1; i < Nx - 1; i++) {
				int ij = i + j * Nx;
				int ip1j = i + 1 + j * Nx;
				int ijp1 = i + (j + 1) * Nx;
				
				// x-direction fluxes
				vL = sqrt(gamma * P_Rx[ij] / rho_Rx[ij]);
				vL_2 = sqrt(vx_Rx[ij]*vx_Rx[ij] + vy_Rx[ij]*vy_Rx[ij]);
				vR = sqrt(gamma * P_Lx[ip1j] / rho_Lx[ip1j]);
				vR_2 = sqrt(vx_Lx[ip1j]*vx_Lx[ip1j] + vy_Lx[ip1j]*vy_Lx[ip1j]);
				vmax_x = std::max(vL_2 + vL, vR_2 + vR);

				F_mass_L = rho_Rx[ij] * vx_Rx[ij];
				F_mass_R = rho_Lx[ip1j] * vx_Lx[ip1j];
				F_momx_L = rho_Rx[ij] * vx_Rx[ij] * vx_Rx[ij] + P_Rx[ij];
				F_momx_R = rho_Lx[ip1j] * vx_Lx[ip1j] * vx_Lx[ip1j] + P_Lx[ip1j];
				F_momy_L = rho_Rx[ij] * vx_Rx[ij] * vy_Rx[ij];
				F_momy_R = rho_Lx[ip1j] * vx_Lx[ip1j] * vy_Lx[ip1j];
				F_energy_L = (P_Rx[ij] / (gamma - 1) + 0.5 * rho_Rx[ij] * (square(vx_Rx[ij]) + square(vy_Rx[ij])) + P_Rx[ij]) * vx_Rx[ij];
				F_energy_R = (P_Lx[ip1j] / (gamma - 1) + 0.5 * rho_Lx[ip1j] * (square(vx_Lx[ip1j]) + square(vy_Lx[ip1j])) + P_Lx[ip1j]) * vx_Lx[ip1j];
				Q_mass_L = rho_Rx[ij];
				Q_mass_R = rho_Lx[ip1j];
				Q_momx_L = rho_Rx[ij] * vx_Rx[ij];
				Q_momx_R = rho_Lx[ip1j] * vx_Lx[ip1j];
				Q_momy_L = rho_Rx[ij] * vy_Rx[ij];
				Q_momy_R = rho_Lx[ip1j] * vy_Lx[ip1j];
				Q_energy_L = (P_Rx[ij] / (gamma - 1) + 0.5 * rho_Rx[ij] * (square(vx_Rx[ij]) + square(vy_Rx[ij])));
				Q_energy_R = (P_Lx[ip1j] / (gamma - 1) + 0.5 * rho_Lx[ip1j] * (square(vx_Lx[ip1j]) + square(vy_Lx[ip1j])));

				mass_flux_x[ij] = 0.5 * (F_mass_L + F_mass_R - vmax_x * (Q_mass_R - Q_mass_L));
				momx_flux_x[ij] = 0.5 * (F_momx_L + F_momx_R - vmax_x * (Q_momx_R - Q_momx_L));
				momy_flux_x[ij] = 0.5 * (F_momy_L + F_momy_R - vmax_x * (Q_momy_R - Q_momy_L));
				energy_flux_x[ij] = 0.5 * (F_energy_L + F_energy_R - vmax_x * (Q_energy_R - Q_energy_L));

				// y-direction fluxes
				vL = sqrt(gamma * P_Ry[ij] / rho_Ry[ij]);
				vL_2 = sqrt(vx_Ry[ij]*vx_Ry[ij] + vy_Ry[ij]*vy_Ry[ij]);
				vR = sqrt(gamma * P_Ly[ijp1] / rho_Ly[ijp1]);
				vR_2 = sqrt(vx_Ly[ijp1]*vx_Ly[ijp1] + vy_Ly[ijp1]*vy_Ly[ijp1]);
				vmax_y = std::max(vL_2 + vL, vR_2 + vR);

				F_mass_L = rho_Ry[ij] * vy_Ry[ij];
				F_mass_R = rho_Ly[ijp1] * vy_Ly[ijp1];
				F_momx_L = rho_Ry[ij] * vy_Ry[ij] * vx_Ry[ij];
				F_momx_R = rho_Ly[ijp1] * vy_Ly[ijp1] * vx_Ly[ijp1];
				F_momy_L = rho_Ry[ij] * vy_Ry[ij] * vy_Ry[ij] + P_Ry[ij];
				F_momy_R = rho_Ly[ijp1] * vy_Ly[ijp1] * vy_Ly[ijp1] + P_Ly[ijp1];
				F_energy_L = (P_Ry[ij] / (gamma - 1) + 0.5 * rho_Ry[ij] * (square(vx_Ry[ij]) + square(vy_Ry[ij])) + P_Ry[ij]) * vy_Ry[ij];
				F_energy_R = (P_Ly[ijp1] / (gamma - 1) + 0.5 * rho_Ly[ijp1] * (square(vx_Ly[ijp1]) + square(vy_Ly[ijp1])) + P_Ly[ijp1]) * vy_Ly[ijp1];
				Q_mass_L = rho_Ry[ij];
				Q_mass_R = rho_Ly[ijp1];
				Q_momx_L = rho_Ry[ij] * vx_Ry[ij];
				Q_momx_R = rho_Ly[ijp1] * vx_Ly[ijp1];
				Q_momy_L = rho_Ry[ij] * vy_Ry[ij];
				Q_momy_R = rho_Ly[ijp1] * vy_Ly[ijp1];
				Q_energy_L = (P_Ry[ij] / (gamma - 1) + 0.5 * rho_Ry[ij] * (square(vx_Ry[ij]) + square(vy_Ry[ij])));
				Q_energy_R = (P_Ly[ijp1] / (gamma - 1) + 0.5 * rho_Ly[ijp1] * (square(vx_Ly[ijp1]) + square(vy_Ly[ijp1])));

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

	void update_conserved(double dt) {
		// DONE: update the conserved variables using the fluxes
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
	double Lx, Ly;
	double dx, dy;
	double gamma, output_dt;
	std::vector<double> rho, vx, vy, P;			// primitive variables
	std::vector<double> mass, mom_x, mom_y, energy;		// conserved variables
	// arrays to hold the results during primitive_update
	std::vector<double> rho_tmp, vx_tmp, vy_tmp, P_tmp;
	// arrays of fluxes for each conserved variable:
	std::vector<double> mass_flux_x, mass_flux_y;
	std::vector<double> momx_flux_x, momx_flux_y;
	std::vector<double> momy_flux_x, momy_flux_y;
	std::vector<double> energy_flux_x, energy_flux_y;
	// arrays for extrapolating to cell interfaces:
	std::vector<double> rho_Lx, rho_Ly, rho_Rx, rho_Ry;
	std::vector<double> vx_Lx, vx_Ly, vx_Rx, vx_Ry;
	std::vector<double> vy_Lx, vy_Ly, vy_Rx, vy_Ry;
	std::vector<double> P_Lx, P_Ly, P_Rx, P_Ry;
};
