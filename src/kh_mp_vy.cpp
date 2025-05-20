#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include "kh_mp_vy.h"
#include <chrono>
using namespace std;

const double a = 0.25;
const double b = 0.75;
const double rho_d = 4;
const double rho_0 = 1;
const double v0 = 0.5;
const double P0 = 2.5;
const double k = 4.0 * M_PI;
const double v_small = 0.01;
const double sigma = 0.05;
const int N = 514;
const double end_t = 2.0;

int main(){
    // start timer  
    auto start = std::chrono::high_resolution_clock::now();
    // delete file content
    std::ofstream outfile("../data/kelvin_helmholtz/kh_output_vy.csv", std::ios_base::trunc);
    outfile.close();
    // simulation parameters
    int Nx = N;
    int Ny = N;
    double Lx = 1.0;
    double Ly = 1.0;
    double gamma = 5.0 / 3.0;
    double output_dt = 0.01;

    fluid_solver_2d solver(Lx, Ly, Nx, Ny, gamma, output_dt);

    std::vector<double> rho(Nx * Ny);
    std::vector<double> vx(Nx * Ny);
    std::vector<double> vy(Nx * Ny);
    std::vector<double> P(Nx * Ny);

    int max_threads = omp_get_max_threads();
    cout << "Threads Used: " << max_threads << endl;

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
    solver.solve(0.0, end_t);

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " s" << endl;
    return 0; 
}
