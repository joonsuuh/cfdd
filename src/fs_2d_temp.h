#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

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
    // TODO: Compute conserved variables from primitive ones
    for(int j = 1; j < Ny - 1; j++){
      for(int i = 1; i < Nx - 1; i++){
        int pos = i + j * Nx;
        mass[pos] = rho[pos] * dx * dy;
        mom_x[pos] = rho[pos] * vx[pos] * dx * dy;
        mom_y[pos] = rho[pos] * vy[pos] * dx * dy;
        energy[pos] = rho[pos] * (0.5 * (vx[pos] * vx[pos] + vy[pos] * vy[pos]) + P[pos]/((gamma - 1) * rho[pos])) * dx * dy;
      }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }

  void conserved_to_primitive() {
    // TODO: Compute primitive variables from conserved ones
    for(int j = 1; j < Ny - 1; j++){
      for(int i = 1; i < Nx - 1; i++){
        int pos = i + j * Nx;
        rho[pos] = mass[pos] / (dx * dy);
        vx[pos] = mom_x[pos] / mass[pos];
        vy[pos] = mom_y[pos] / mass[pos];
        P[pos] = ((energy[pos]/(dx * dy))/rho[pos] - 0.5 * (vx[pos] * vx[pos] + vy[pos] * vy[pos])) * (gamma - 1) * rho[pos];
      }
    }
    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);
  }

  void init(const std::vector<double> &rho0, const std::vector<double> &vx0,
            const std::vector<double> &vy0, const std::vector<double> &P0) {
    // TODO: Initialize the primitive variables using the given initial
    // condition
    rho = rho0;
    vx = vx0;
    vy = vy0;
    P = P0;
    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);
    primitive_to_conserved();
  }

  double find_dt() {
    double C_CFL = 0.4;
    double max_val = 0;
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            int pos = i + j * Nx;

            if (rho[pos] <= 0.0 || P[pos] <= 0.0) {
                std::cerr << "Invalid rho or P at pos " << pos
                          << ": rho=" << rho[pos] << ", P=" << P[pos] << std::endl;
                exit(1);
            }

            double curr_val = std::sqrt(gamma * P[pos] / rho[pos]) +
                              std::sqrt(vx[pos] * vx[pos] + vy[pos] * vy[pos]);
            //std::cout << P[pos] << std::endl; //P is bad :(
            max_val = std::max(max_val, curr_val);
        }
    }

    if (max_val <= 0.0) {
        std::cerr << "Invalid max_val: " << max_val << std::endl;
        exit(1);
    }

    return C_CFL * std::min(dx, dy) / max_val;
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
      std::cout << "t: " << t << ", dt: " << dt << ", n: " << n << ", output_dt * n: " << output_dt * n << std::endl;
      if (dt <= 0.0 || std::isnan(dt) || std::isinf(dt)) {
        std::cerr << "Invalid dt detected: " << dt << std::endl;
        exit(1); // Stop the program if dt is invalid
      }

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
    for (int j = 1; j < Ny - 1; j++) {
        f[j * Nx] = f[j * Nx + Nx - 2];
        f[j * Nx + Nx - 1] = f[j * Nx + 1];
    }
    for (int i = 1; i < Nx - 1; i++) {
        f[i] = f[(Ny - 2) * Nx + i];
        f[(Ny - 1) * Nx + i] = f[i + Nx];
    }
  }


  void primitive_update(double dt) {
  // TODO: update the primitive variables using Euler equations in primitive
  // form using an FTCS scheme
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            int curr_pos = i + j * Nx;
            int i_plus = (i + 1) + j * Nx;
            int i_minus = (i - 1) + j * Nx;
            int j_plus = i + (j + 1) * Nx;
            int j_minus = i + (j - 1) * Nx;

            // Calculate derivatives using centered differences
            double rho_x_diff = (rho[i_plus] - rho[i_minus]) / (2.0 * dx);
            double rho_y_diff = (rho[j_plus] - rho[j_minus]) / (2.0 * dy);
            double vx_x_diff = (vx[i_plus] - vx[i_minus]) / (2.0 * dx);
            double vx_y_diff = (vx[j_plus] - vx[j_minus]) / (2.0 * dy);
            double vy_x_diff = (vy[i_plus] - vy[i_minus]) / (2.0 * dx);
            double vy_y_diff = (vy[j_plus] - vy[j_minus]) / (2.0 * dy);

            // Update rho using FTCS
            rho_tmp[curr_pos] = rho[curr_pos] - 0.5 * dt * (
                vx[curr_pos] * rho_x_diff + vy[curr_pos] * rho_y_diff +
                rho[curr_pos] * (vx_x_diff + vy_y_diff)
            );

            // Update vx
            vx_tmp[curr_pos] = vx[curr_pos] - 0.5 * dt * (
                vx[curr_pos] * vx_x_diff + vy[curr_pos] * vx_y_diff +
                (1.0 / rho[curr_pos]) * (P[i_plus] - P[i_minus]) / (2.0 * dx)
            );

            // Update vy
            vy_tmp[curr_pos] = vy[curr_pos] - 0.5 * dt * (
                vx[curr_pos] * vy_x_diff + vy[curr_pos] * vy_y_diff +
                (1.0 / rho[curr_pos]) * (P[j_plus] - P[j_minus]) / (2.0 * dy)
            );

            // Update P
            P_tmp[curr_pos] = P[curr_pos] - 0.5 * dt * (
                gamma * P[curr_pos] * (vx_x_diff + vy_y_diff) +
                vx[curr_pos] * (P[i_plus] - P[i_minus]) / (2.0 * dx) +
                vy[curr_pos] * (P[j_plus] - P[j_minus]) / (2.0 * dy)
            );
        }
    }

    std::swap(rho, rho_tmp);
    std::swap(vx, vx_tmp);
    std::swap(vy, vy_tmp);
    std::swap(P, P_tmp);
    periodic_boundary(rho);
    periodic_boundary(vx);
    periodic_boundary(vy);
    periodic_boundary(P);
  }
  void extrapolate_to_interface() {
    // TODO: compute rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, and P_R here
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            int curr_pos = i + j * Nx;
            int i_plus = (i + 1) + j * Nx;
            int i_minus = (i - 1) + j * Nx;
            int j_plus = i + (j + 1) * Nx;
            int j_minus = i + (j - 1) * Nx;

            // Extrapolate rho
            rho_Lx[curr_pos] = rho[curr_pos] - 0.25 * (rho[i_plus] - rho[i_minus]);
            rho_Rx[curr_pos] = rho[curr_pos] + 0.25 * (rho[i_plus] - rho[i_minus]);
            rho_Ly[curr_pos] = rho[curr_pos] - 0.25 * (rho[j_plus] - rho[j_minus]);
            rho_Ry[curr_pos] = rho[curr_pos] + 0.25 * (rho[j_plus] - rho[j_minus]);

            // Extrapolate vx
            vx_Lx[curr_pos] = vx[curr_pos] - 0.25 * (vx[i_plus] - vx[i_minus]);
            vx_Rx[curr_pos] = vx[curr_pos] + 0.25 * (vx[i_plus] - vx[i_minus]);
            vx_Ly[curr_pos] = vx[curr_pos] - 0.25 * (vx[j_plus] - vx[j_minus]);
            vx_Ry[curr_pos] = vx[curr_pos] + 0.25 * (vx[j_plus] - vx[j_minus]);

            // Extrapolate vy
            vy_Lx[curr_pos] = vy[curr_pos] - 0.25 * (vy[i_plus] - vy[i_minus]);
            vy_Rx[curr_pos] = vy[curr_pos] + 0.25 * (vy[i_plus] - vy[i_minus]);
            vy_Ly[curr_pos] = vy[curr_pos] - 0.25 * (vy[j_plus] - vy[j_minus]);
            vy_Ry[curr_pos] = vy[curr_pos] + 0.25 * (vy[j_plus] - vy[j_minus]);

            // Extrapolate P
            P_Lx[curr_pos] = P[curr_pos] - 0.25 * (P[i_plus] - P[i_minus]);
            P_Rx[curr_pos] = P[curr_pos] + 0.25 * (P[i_plus] - P[i_minus]);
            P_Ly[curr_pos] = P[curr_pos] - 0.25 * (P[j_plus] - P[j_minus]);
            P_Ry[curr_pos] = P[curr_pos] + 0.25 * (P[j_plus] - P[j_minus]);
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
    // TODO: compute the fluxes
    extrapolate_to_interface();
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            int curr_pos = i + j * Nx;

            //Calculate v_max for x and y
            double v_max_x = std::max(
              std::sqrt(gamma * P_Rx[curr_pos] / rho_Rx[curr_pos]) + std::sqrt(vx_Rx[curr_pos] * vx_Rx[curr_pos] + vy_Rx[curr_pos] * vy_Rx[curr_pos]),
              std::sqrt(gamma * P_Lx[curr_pos + 1] / rho_Lx[curr_pos + 1]) + std::sqrt(vx_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1] + vy_Lx[curr_pos + 1] * vy_Lx[curr_pos + 1])
            );
            double v_max_y = std::max(
              std::sqrt(gamma * P_Ry[curr_pos] / rho_Ry[curr_pos]) + std::sqrt(vx_Ry[curr_pos] * vx_Ry[curr_pos] + vy_Ry[curr_pos] * vy_Ry[curr_pos]),
              std::sqrt(gamma * P_Ly[curr_pos + Nx] / rho_Ly[curr_pos + Nx]) + std::sqrt(vx_Ly[curr_pos + Nx] * vx_Ly[curr_pos + Nx] + vy_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx])
            );

            //Calculate rho flux stuff
            double rho_F_L_x = rho_Rx[curr_pos] * vx_Rx[curr_pos];
            double rho_F_R_x = rho_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1];
            double rho_F_L_y = rho_Ry[curr_pos] * vy_Ry[curr_pos];
            double rho_F_R_y = rho_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx];
            double rho_Q_L_x = rho_Rx[curr_pos];
            double rho_Q_R_x = rho_Lx[curr_pos + 1];
            double rho_Q_L_y = rho_Ry[curr_pos];
            double rho_Q_R_y = rho_Ly[curr_pos + Nx];

            double rho_flux_x = 0.5 * (rho_F_L_x + rho_F_R_x) - v_max_x / 2.0 * (rho_Q_R_x - rho_Q_L_x);
            double rho_flux_y = 0.5 * (rho_F_L_y + rho_F_R_y) - v_max_y / 2.0 * (rho_Q_R_y - rho_Q_L_y);

            mass_flux_x[curr_pos] = rho_flux_x * dx * dy;
            mass_flux_y[curr_pos] = rho_flux_y * dx * dy;

            //Calculate rhovx flux stuff
            double rhovx_F_L_x = rho_Rx[curr_pos] * vx_Rx[curr_pos] * vx_Rx[curr_pos] + P_Rx[curr_pos];
            double rhovx_F_R_x = rho_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1] + P_Lx[curr_pos + 1];
            double rhovx_F_L_y = rho_Ry[curr_pos] * vx_Ry[curr_pos] * vy_Ry[curr_pos];
            double rhovx_F_R_y = rho_Ly[curr_pos + Nx] * vx_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx];
            double rhovx_Q_L_x = rho_Rx[curr_pos] * vx_Rx[curr_pos];
            double rhovx_Q_R_x = rho_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1];
            double rhovx_Q_L_y = rho_Ry[curr_pos] * vx_Ry[curr_pos];
            double rhovx_Q_R_y = rho_Ly[curr_pos + Nx] * vx_Ly[curr_pos + Nx];

            double rhovx_flux_x = 0.5 * (rhovx_F_L_x + rhovx_F_R_x) - v_max_x / 2.0 * (rhovx_Q_R_x - rhovx_Q_L_x);
            double rhovx_flux_y = 0.5 * (rhovx_F_L_y + rhovx_F_R_y) - v_max_y / 2.0 * (rhovx_Q_R_y - rhovx_Q_L_y);

            momx_flux_x[curr_pos] = rhovx_flux_x * dx * dy;
            momx_flux_y[curr_pos] = rhovx_flux_y * dx * dy;

            //Calculate rhovy flux stuff
            double rhovy_F_L_x = rho_Rx[curr_pos] * vx_Rx[curr_pos] * vy_Rx[curr_pos];
            double rhovy_F_R_x = rho_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1] * vy_Lx[curr_pos + 1];
            double rhovy_F_L_y = rho_Ry[curr_pos] * vy_Ry[curr_pos] * vy_Ry[curr_pos] + P_Ry[curr_pos];
            double rhovy_F_R_y = rho_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx] + P_Ly[curr_pos + Nx];
            double rhovy_Q_L_x = rho_Rx[curr_pos] * vy_Rx[curr_pos];
            double rhovy_Q_R_x = rho_Lx[curr_pos + 1] * vy_Lx[curr_pos + 1];
            double rhovy_Q_L_y = rho_Ry[curr_pos] * vy_Ry[curr_pos];
            double rhovy_Q_R_y = rho_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx];

            double rhovy_flux_x = 0.5 * (rhovy_F_L_x + rhovy_F_R_x) - v_max_x / 2.0 * (rhovy_Q_R_x - rhovy_Q_L_x);
            double rhovy_flux_y = 0.5 * (rhovy_F_L_y + rhovy_F_R_y) - v_max_y / 2.0 * (rhovy_Q_R_y - rhovy_Q_L_y);

            momy_flux_x[curr_pos] = rhovy_flux_x * dx * dy;
            momy_flux_y[curr_pos] = rhovy_flux_y * dx * dy;

            //Calculate U flux stuff
            //These first 4 variables are based off direction to the interface, not our current position.
            double U_Lx = rho_Lx[curr_pos + 1] * (0.5 * vx_Lx[curr_pos + 1] * vx_Lx[curr_pos + 1] + vy_Lx[curr_pos + 1] * vy_Lx[curr_pos + 1] + P_Lx[curr_pos + 1] / (gamma - 1) / rho_Lx[curr_pos + 1]);
            double U_Rx = rho_Rx[curr_pos] * (0.5 * vx_Rx[curr_pos] * vx_Rx[curr_pos] + vy_Rx[curr_pos] * vy_Rx[curr_pos] + P_Rx[curr_pos] / (gamma - 1) / rho_Rx[curr_pos]);
            double U_Ly = rho_Ly[curr_pos + Nx] * (0.5 * vx_Ly[curr_pos + Nx] * vx_Ly[curr_pos + Nx] + vy_Ly[curr_pos + Nx] * vy_Ly[curr_pos + Nx] + P_Ly[curr_pos + Nx] / (gamma - 1) / rho_Ly[curr_pos + Nx]);
            double U_Ry = rho_Ry[curr_pos] * (0.5 * vx_Ry[curr_pos] * vx_Ry[curr_pos] + vy_Ry[curr_pos] * vy_Ry[curr_pos] + P_Ry[curr_pos] / (gamma - 1) / rho_Ry[curr_pos]);

            double U_F_L_x = (U_Rx + P_Rx[curr_pos]) * vx_Rx[curr_pos];
            double U_F_R_x = (U_Lx + P_Lx[curr_pos + 1]) * vx_Lx[curr_pos + 1];
            double U_F_L_y = (U_Ry + P_Ry[curr_pos]) * vy_Ry[curr_pos];
            double U_F_R_y = (U_Ly + P_Ly[curr_pos + Nx]) * vy_Ly[curr_pos + Nx];
            double U_Q_L_x = U_Rx;
            double U_Q_R_x = U_Lx;
            double U_Q_L_y = U_Ry;
            double U_Q_R_y = U_Ly;

            double U_flux_x = 0.5 * (U_F_L_x + U_F_R_x) - v_max_x / 2.0 * (U_Q_R_x - U_Q_L_x);
            double U_flux_y = 0.5 * (U_F_L_y + U_F_R_y) - v_max_y / 2.0 * (U_Q_R_y - U_Q_L_y);

            energy_flux_x[curr_pos] = U_flux_x * dx * dy;
            energy_flux_y[curr_pos] = U_flux_y * dx * dy;
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
    // TODO: update the conserved variables using the fluxes
    for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
            int curr_pos = i + j * Nx;
            int i_minus = (i - 1) + j * Nx;
            int j_minus = i + (j - 1) * Nx;

            mass[curr_pos] -= dt * ((mass_flux_x[curr_pos] - mass_flux_x[i_minus]) / dx + (mass_flux_y[curr_pos] - mass_flux_y[j_minus]) / dy);
            mom_x[curr_pos] -= dt * ((momx_flux_x[curr_pos] - momx_flux_x[i_minus]) / dx + (momx_flux_y[curr_pos] - momx_flux_y[j_minus]) / dy);
            mom_y[curr_pos] -= dt * ((momy_flux_x[curr_pos] - momy_flux_x[i_minus]) / dx + (momy_flux_y[curr_pos] - momy_flux_y[j_minus]) / dy);
            energy[curr_pos] -= dt * ((energy_flux_x[curr_pos] - energy_flux_x[i_minus]) / dx + (energy_flux_y[curr_pos] - energy_flux_y[j_minus]) / dy);
        }
    }
    periodic_boundary(mass);
    periodic_boundary(mom_x);
    periodic_boundary(mom_y);
    periodic_boundary(energy);
  }

  void output(int n) {
    std::ofstream outfile("output_rho_" + std::to_string(n) + ".csv");
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
  std::vector<double> rho, vx, vy, P;                 // primitive variables
  std::vector<double> mass, mom_x, mom_y, energy;     // conserved variables
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