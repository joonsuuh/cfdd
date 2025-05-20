#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "fs_2d_temp.h"
#include <chrono>
using namespace std;

const double a = 0.25;
const double b = 0.75;
const double rho_d = 4;
const double rho_0 = 1;
const double v0 = 0.5;
const double P0 = 2.5;
const double k = 6 * M_PI;
const double v_small = 0.01;
const double sigma = 0.05;

int main(){
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // simulation parameters
    int Nx = 512;
    int Ny = 512;
    double Lx = 1.0;
    double Ly = 1.0;
    double gamma = 5.0 / 3.0;
    double output_dt = 0.01;

    fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, output_dt);

    std::vector<double> rho(Nx * Ny);
    std::vector<double> vx(Nx * Ny);
    std::vector<double> vy(Nx * Ny);
    std::vector<double> P(Nx * Ny);

    double dx = Lx / (Nx - 2);
    double dy = Ly / (Ny - 2);
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int idx = i + j * Nx;

            double x = i * dx;
            double y = j * dy;
            P[idx] = P0;

            if (a <= y && y <= b){
                rho[idx] = rho_d;
                vx[idx] = v0;
            }
            else if (y < a || y > b){
                rho[idx] = rho_0;
                vx[idx] = -v0; 
            }
            vy[idx] = v_small * sin(k*x) * (exp(-(y-a)*(y-a) / sigma / sigma) + exp(-(y-b)*(y-b) / sigma / sigma));
        }
    }

    solver.init(rho, vx, vy, P);
    solver.solve(0.0, 2.5);

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << " s" << endl;

    return 0;
}
