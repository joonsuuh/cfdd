#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include "fs_2d_mp_fast.h"
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
const double cfl = 0.4;

int main(){
    // start timer  
    auto start = std::chrono::high_resolution_clock::now();

    // set seed
    srand(1337);

    // simulation parameters
    int Nx = 514;
    int Ny = 514;
    double Lx = 1.0;
    double Ly = 1.0;
    double gamma = 5.0 / 3.0;
    double output_dt = 0.01;

    fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, output_dt, cfl);

    std::vector<double> rho(Nx * Ny);
    std::vector<double> vx(Nx * Ny);
    std::vector<double> vy(Nx * Ny);
    std::vector<double> P(Nx * Ny);

    int max_threads = omp_get_max_threads();
    cout << "Threads Used: " << max_threads << endl;

    // output vy to file
    std::ofstream outfile("../data/kelvin_helmholtz/init_vy.csv");
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
            // replace sin(kx) with random number generator to roll from 0 to 12 pi
            double k_rand = ((double)rand() / RAND_MAX) * 12 * M_PI;
            vy[idx] = v_small * sin(k_rand * x) * (exp(-(y-a)*(y-a) / sigma / sigma) + exp(-(y-b)*(y-b) / sigma / sigma));

            // output vy to file when y = 0.75
            // std::ofstream outfile("../data/kelvin_helmholtz/init_vy.csv", std::ios_base::app);
            if (j == 385){
                outfile << vy[idx];
                if (i != Nx - 1)
                    outfile << ", ";
                else
                    outfile << "\n";
            }
        }
    }
    outfile.close();

    solver.init(rho, vx, vy, P);
    solver.solve(0.0, 1.0);

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " s" << endl;
    return 0; 
}
