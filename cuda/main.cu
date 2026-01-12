#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "fs2d.cuh"
#include <chrono>

const float a = 0.25f;
const float b = 0.75f;
const float rho_d = 4.0f;
const float rho_0 = 1.0f;
const float v0 = 0.5f;
const float P0 = 2.5f;
const float k = 2.0f * M_PI;
const float v_small = 0.01f;
const float sigma = 0.05f;

int main(){
    // start timer  
    auto start = std::chrono::high_resolution_clock::now();

    // simulation parameters
    int Nx = 514;
    int Ny = 514;
    float Lx = 1.0f;
    float Ly = 1.0f;
    float gamma = 5.0 / 3.0f;
    float output_dt = 0.01f;

    FluidSolver2D solver(Lx, Ly, Nx, Ny, gamma, output_dt);

    std::vector<float> rho(Nx * Ny);
    std::vector<float> vx(Nx * Ny);
    std::vector<float> vy(Nx * Ny);
    std::vector<float> P(Nx * Ny);

    float dx = Lx / (Nx - 2);
    float dy = Ly / (Ny - 2);
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int idx = i + j * Nx;

            float x = (i - 1) * dx;
            float y = (j - 1) * dy;
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
    solver.solve(0.0f, 5.0f);
    cudaDeviceSynchronize();

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed = end - start;
    std::cout << "Time: " << elapsed.count() << " s" << std::endl;
    return 0; 
}