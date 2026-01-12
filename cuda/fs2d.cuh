#include <cassert>
// #include <cuda_runtime.h>

#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) 
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
#pragma once

__device__ __forceinline__ float square(float x) { return x * x; }

__global__ void primitive_to_conservative(float* rho, float* mass, float* vx, float* vy,
                                          float* momx, float* momy, float* P, float* energy,
                                          float dx, float dy, float gamma, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int idx = i + j * nx;

    if (i < 0 || i >= nx || j < 0 || j >= ny) { return; }

    float const vol = dx * dy;
    float const rho_c = fmaxf(rho[idx], 1e-12f);
    float const mass_c = rho_c * vol;
    float const vx_c = vx[idx];
    float const vy_c = vy[idx];

    mass[idx] = mass_c;
    momx[idx] = mass_c * vx_c;
    momy[idx] = mass_c * vy_c;

    float const KE = 0.5f * (vx_c * vx_c + vy_c * vy_c);
    float const PE = P[idx] / (gamma - 1.0f) / rho_c; 
    energy[idx] = mass_c * (PE + KE);
}

__global__ void conserved_to_primitive_kernel(float* rho, float* mass, float* vx, float* vy,
                                       float* momx, float* momy, float* P, float* energy,
                                       float dx, float dy, float gamma, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int idx = i + j * nx;

    if (i < 0 || i > nx - 1 || j < 0 || j > ny - 1) { return; }

    float const vol = dx * dy;
    float const mass_c = fmaxf(mass[idx], 1e-12f);

    float const rho_c = mass_c / vol;
    float const vx_c = momx[idx] / mass_c;
    float const vy_c = momy[idx] / mass_c;

    rho[idx] = rho_c;
    vx[idx] = vx_c;
    vy[idx] = vy_c;
    float const KE = 0.5f * (vx_c * vx_c + vy_c * vy_c);
    P[idx] = fmaxf((energy[idx] / mass_c - KE) * rho_c * (gamma - 1.0f), 1e-12f);
}

__global__ void update_conserved_kernel(
    float* mass, float* momx, float* momy, float* energy,
    float* mass_flux_x, float* momx_flux_x, float* momy_flux_x, float* energy_flux_x,
    float* mass_flux_y, float* momx_flux_y, float* momy_flux_y, float* energy_flux_y,
    int nx, int ny, float dx, float dy, float dt)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < 1 || i > nx - 2 || j < 1 || j > ny - 2) { return; }

    int ij = i + j * nx;
    int i_minus1_j = i - 1 + j * nx;
    int i_j_minus1 = i + (j - 1) * nx;

    mass[ij] = mass[ij] - ((mass_flux_x[ij] - mass_flux_x[i_minus1_j]) * dy * dt)
        - ((mass_flux_y[ij] - mass_flux_y[i_j_minus1]) * dx * dt);
    momx[ij] = momx[ij] - ((momx_flux_x[ij] - momx_flux_x[i_minus1_j]) * dy * dt)
        - ((momx_flux_y[ij] - momx_flux_y[i_j_minus1]) * dx * dt);
    momy[ij] = momy[ij] - ((momy_flux_x[ij] - momy_flux_x[i_minus1_j]) * dy * dt)
        - ((momy_flux_y[ij] - momy_flux_y[i_j_minus1]) * dx * dt);
    energy[ij] = energy[ij] - ((energy_flux_x[ij] - energy_flux_x[i_minus1_j]) * dy * dt)
        - ((energy_flux_y[ij] - energy_flux_y[i_j_minus1]) * dx * dt);
}

__global__ void compute_fluxes_kernel(
    float* mass_flux_x, float* momx_flux_x, float* momy_flux_x, float* energy_flux_x,
    float* mass_flux_y, float* momx_flux_y, float* momy_flux_y, float* energy_flux_y,
    float* rho_lx, float* rho_ly, float* rho_rx, float* rho_ry,
    float* vx_lx, float* vx_ly, float* vy_lx, float* vy_ly,
    float* vx_rx, float* vx_ry, float* vy_rx, float* vy_ry,
    float* P_lx, float* P_rx, float* P_ly, float* P_ry, float gamma, int nx, int ny)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
    {
        int ij = i + j * nx;
        int ip1j = i + 1 + j * nx;
        int ijp1 = i + (j + 1) * nx;

        float vl, vr, vl_2, vr_2, vmax_x, vmax_y;
        float F_mass_l, F_mass_r, F_momx_l, F_momx_r, F_momy_l, F_momy_r, F_energy_l, F_energy_r;
        float Q_mass_l, Q_mass_r, Q_momx_l, Q_momx_r, Q_momy_l, Q_momy_r, Q_energy_l, Q_energy_r;
        // x-direction fluxes
        float rhor_x = rho_rx[ij];
        float rhol_xp = rho_lx[ip1j];
        vl = sqrtf(gamma * P_rx[ij] / rhor_x);
        vl_2 = sqrtf(vx_rx[ij]*vx_rx[ij] + vy_rx[ij]*vy_rx[ij]);
        vr = sqrtf(gamma * P_lx[ip1j] / rhol_xp);
        vr_2 = sqrtf(vx_lx[ip1j]*vx_lx[ip1j] + vy_lx[ip1j]*vy_lx[ip1j]);
        vmax_x = fmaxf(vl_2 + vl, vr_2 + vr);

        F_mass_l = rhor_x * vx_rx[ij];
        F_mass_r = rhol_xp * vx_lx[ip1j];
        F_momx_l = rhor_x * vx_rx[ij] * vx_rx[ij] + P_rx[ij];
        F_momx_r = rhol_xp * vx_lx[ip1j] * vx_lx[ip1j] + P_lx[ip1j];
        F_momy_l = rhor_x * vx_rx[ij] * vy_rx[ij];
        F_momy_r = rhol_xp * vx_lx[ip1j] * vy_lx[ip1j];
        F_energy_l = (P_rx[ij] / (gamma - 1) + 0.5f * rhor_x * (square(vx_rx[ij]) + square(vy_rx[ij])) + P_rx[ij]) * vx_rx[ij];
        F_energy_r = (P_lx[ip1j] / (gamma - 1) + 0.5f * rhol_xp * (square(vx_lx[ip1j]) + square(vy_lx[ip1j])) + P_lx[ip1j]) * vx_lx[ip1j];
        Q_mass_l = rhor_x;
        Q_mass_r = rhol_xp;
        Q_momx_l = rhor_x * vx_rx[ij];
        Q_momx_r = rhol_xp * vx_lx[ip1j];
        Q_momy_l = rhor_x * vy_rx[ij];
        Q_momy_r = rhol_xp * vy_lx[ip1j];
        Q_energy_l = (P_rx[ij] / (gamma - 1) + 0.5f * rhor_x * (square(vx_rx[ij]) + square(vy_rx[ij])));
        Q_energy_r = (P_lx[ip1j] / (gamma - 1) + 0.5f * rhol_xp * (square(vx_lx[ip1j]) + square(vy_lx[ip1j])));

        mass_flux_x[ij] = 0.5f * (F_mass_l + F_mass_r - vmax_x * (Q_mass_r - Q_mass_l));
        momx_flux_x[ij] = 0.5f * (F_momx_l + F_momx_r - vmax_x * (Q_momx_r - Q_momx_l));
        momy_flux_x[ij] = 0.5f * (F_momy_l + F_momy_r - vmax_x * (Q_momy_r - Q_momy_l));
        energy_flux_x[ij] = 0.5f * (F_energy_l + F_energy_r - vmax_x * (Q_energy_r - Q_energy_l));

        // y-direction fluxes
        float rhor_y = rho_ry[ij];
        float rhol_yp = rho_ly[ijp1];
        vl = sqrtf(gamma * P_ry[ij] / rhor_y);
        vl_2 = sqrtf(vx_ry[ij]*vx_ry[ij] + vy_ry[ij]*vy_ry[ij]);
        vr = sqrtf(gamma * P_ly[ijp1] / rhol_yp);
        vr_2 = sqrtf(vx_ly[ijp1]*vx_ly[ijp1] + vy_ly[ijp1]*vy_ly[ijp1]);
        vmax_y = fmaxf(vl_2 + vl, vr_2 + vr);

        F_mass_l = rhor_y * vy_ry[ij];
        F_mass_r = rhol_yp * vy_ly[ijp1];
        F_momx_l = rhor_y * vy_ry[ij] * vx_ry[ij];
        F_momx_r = rhol_yp * vy_ly[ijp1] * vx_ly[ijp1];
        F_momy_l = rhor_y * vy_ry[ij] * vy_ry[ij] + P_ry[ij];
        F_momy_r = rhol_yp * vy_ly[ijp1] * vy_ly[ijp1] + P_ly[ijp1];
        F_energy_l = (P_ry[ij] / (gamma - 1) + 0.5f * rhor_y * (square(vx_ry[ij]) + square(vy_ry[ij])) + P_ry[ij]) * vy_ry[ij];
        F_energy_r = (P_ly[ijp1] / (gamma - 1) + 0.5f * rhol_yp * (square(vx_ly[ijp1]) + square(vy_ly[ijp1])) + P_ly[ijp1]) * vy_ly[ijp1];
        Q_mass_l = rhor_y;
        Q_mass_r = rhol_yp;
        Q_momx_l = rhor_y * vx_ry[ij];
        Q_momx_r = rhol_yp * vx_ly[ijp1];
        Q_momy_l = rhor_y * vy_ry[ij];
        Q_momy_r = rhol_yp * vy_ly[ijp1];
        Q_energy_l = (P_ry[ij] / (gamma - 1) + 0.5f * rhor_y * (square(vx_ry[ij]) + square(vy_ry[ij])));
        Q_energy_r = (P_ly[ijp1] / (gamma - 1) + 0.5f * rhol_yp * (square(vx_ly[ijp1]) + square(vy_ly[ijp1])));

        mass_flux_y[ij] = 0.5f * (F_mass_l + F_mass_r - vmax_y * (Q_mass_r - Q_mass_l));
        momx_flux_y[ij] = 0.5f * (F_momx_l + F_momx_r - vmax_y * (Q_momx_r - Q_momx_l));
        momy_flux_y[ij] = 0.5f * (F_momy_l + F_momy_r - vmax_y * (Q_momy_r - Q_momy_l));
        energy_flux_y[ij] = 0.5f * (F_energy_l + F_energy_r - vmax_y * (Q_energy_r - Q_energy_l));
    }
}

__global__ void extrapolate_to_interface_kernel(
    float* rho, float* vx, float* vy, float* P, float* rho_lx, float* rho_ly, float* rho_rx, float* rho_ry,
    float* vx_lx, float* vx_ly, float* vy_lx, float* vy_ly, float* vx_rx, float* vx_ry, float* vy_rx, float* vy_ry,
    float* P_lx, float* P_rx, float* P_ly, float* P_ry, int nx, int ny)

{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= 1 && i < nx - 1 && j >= 1 && j < ny - 1)
    {
        int ij = i + j * nx;
        int i_plus1_j = i + 1 + j * nx;
        int i_minus1_j = i - 1 + j * nx;
        int i_j_minus1 = i + (j - 1) * nx;
        int i_j_plus1 = i + (j + 1) * nx;

        rho_lx[ij] = rho[ij] - 0.25f * (rho[i_plus1_j] - rho[i_minus1_j]);
        rho_rx[ij] = rho[ij] + 0.25f * (rho[i_plus1_j] - rho[i_minus1_j]);
        rho_ly[ij] = rho[ij] - 0.25f * (rho[i_j_plus1] - rho[i_j_minus1]);
        rho_ry[ij] = rho[ij] + 0.25f * (rho[i_j_plus1] - rho[i_j_minus1]);

        vx_lx[ij] = vx[ij] - 0.25f * (vx[i_plus1_j] - vx[i_minus1_j]);
        vx_rx[ij] = vx[ij] + 0.25f * (vx[i_plus1_j] - vx[i_minus1_j]);
        vx_ly[ij] = vx[ij] - 0.25f * (vx[i_j_plus1] - vx[i_j_minus1]);
        vx_ry[ij] = vx[ij] + 0.25f * (vx[i_j_plus1] - vx[i_j_minus1]);

        vy_lx[ij] = vy[ij] - 0.25f * (vy[i_plus1_j] - vy[i_minus1_j]);
        vy_rx[ij] = vy[ij] + 0.25f * (vy[i_plus1_j] - vy[i_minus1_j]);
        vy_ly[ij] = vy[ij] - 0.25f * (vy[i_j_plus1] - vy[i_j_minus1]);
        vy_ry[ij] = vy[ij] + 0.25f * (vy[i_j_plus1] - vy[i_j_minus1]);

        P_lx[ij] = P[ij] - 0.25f * (P[i_plus1_j] - P[i_minus1_j]);
        P_rx[ij] = P[ij] + 0.25f * (P[i_plus1_j] - P[i_minus1_j]);
        P_ly[ij] = P[ij] - 0.25f * (P[i_j_plus1] - P[i_j_minus1]);
        P_ry[ij] = P[ij] + 0.25f * (P[i_j_plus1] - P[i_j_minus1]);
    }
}

__global__ void primitive_update_kernel(
    float* rho, float* vx, float* vy, float* P,
    float* rho_tmp, float* vx_tmp, float* vy_tmp, float* P_tmp,
    int Nx, int Ny, float dx, float dy, float dt, float gamma) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= 1 && i < Nx - 1 && j >= 1 && j < Ny - 1) {
        int ij = i + j * Nx;
        int i_plus1 = i + 1;
        int i_minus1 = i - 1;
        int j_plus1 = j + 1;
        int j_minus1 = j - 1;
        int i_plus1_j = i_plus1 + j * Nx;
        int i_minus1_j = i_minus1 + j * Nx;
        int i_j_minus1 = i + j_minus1 * Nx;
        int i_j_plus1 = i + j_plus1 * Nx;

        float rho_c = fmaxf(rho[ij], 1e-6f);

        // Continuity equation
        float rho_inner = (vx[ij] * (rho[i_plus1_j] - rho[i_minus1_j]) / (2.0f * dx)) +
            (vy[ij] * (rho[i_j_plus1] - rho[i_j_minus1]) / (2.0f * dy)) +
            (rho[ij] * (((vx[i_plus1_j] - vx[i_minus1_j]) / (2.0f * dx)) +
            ((vy[i_j_plus1] - vy[i_j_minus1]) / (2.0f * dy))));
        rho_tmp[ij] = rho[ij] - 0.5f * dt * rho_inner;

        // Momentum in x
        float vx_adv = (vx[ij] * (vx[i_plus1_j] - vx[i_minus1_j]) / (2.0f * dx)) +
            (vy[ij] * (vx[i_j_plus1] - vx[i_j_minus1]) / (2.0f * dy));
        float dPdx = (P[i_plus1_j] - P[i_minus1_j]) / (2.0f * dx);
        vx_tmp[ij] = vx[ij] - 0.5f * dt * (vx_adv + dPdx / rho_c);

        // Momentum in y
        float vy_adv = (vx[ij] * (vy[i_plus1_j] - vy[i_minus1_j]) / (2.0f * dx)) +
            (vy[ij] * (vy[i_j_plus1] - vy[i_j_minus1]) / (2.0f * dy));
        float dPdy = (P[i_j_plus1] - P[i_j_minus1]) / (2.0f * dy);
        vy_tmp[ij] = vy[ij] - 0.5f * dt * (vy_adv + dPdy / rho_c);

        // Pressure evolution
        float P_inner = (vx[ij] * (P[i_plus1_j] - P[i_minus1_j]) / (2.0f * dx)) +
            (vy[ij] * (P[i_j_plus1] - P[i_j_minus1]) / (2.0f * dy)) +
            (gamma * P[ij] * (((vx[i_plus1_j] - vx[i_minus1_j]) / (2.0f * dx)) +
            ((vy[i_j_plus1] - vy[i_j_minus1]) / (2.0f * dy))));
        P_tmp[ij] = P[ij] - 0.5f * dt * P_inner;
    }
}

__global__ void compute_speed(float* rho, float* vx, float* vy, float* P, float gamma, int nx, int ny, float* speeds)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (i < 1 || i > nx - 2 || j < 1 || j > ny - 2) return;

    int idx = i + j * nx;

    float const rho_c = fmaxf(rho[idx], 1e-6f); // Ensure non-negative density
    float const vx_c = vx[idx];
    float const vy_c = vy[idx];
    float const P_c = fmaxf(P[idx], 1e-6f); // Ensure non-negative pressure

    speeds[i - 1 + (j - 1) * (nx - 2)] = sqrtf(gamma * P_c / rho_c) + sqrtf(vx_c * vx_c + vy_c * vy_c);
}

__device__ float atomicMaxf(float* address, float val) {
    int* address_as_i = (int*)address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
                        __float_as_int(fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

// Modified maxReduceAtomic kernel with boundary checks and error handling
__global__ void maxReduceAtomic(const float* input, float* output, int size,
                                float identity) {
    extern __shared__ float s_input[];
    unsigned int tid = threadIdx.x;
    unsigned int segment = blockIdx.x * blockDim.x * 2;
    unsigned int i = segment + tid;

    if (i + blockDim.x < size) {
        s_input[tid] = fmaxf(input[i], input[i + blockDim.x]);
    } else if (i < size) {
        s_input[tid] = input[i];
    } else {
        s_input[tid] = identity;
    }

    // Reduction in shared memory
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        __syncthreads();

        if (tid < s) {
            s_input[tid] = fmaxf(s_input[tid], s_input[tid + s]);
        }
    }

    // Multi-block reduction
    if (tid == 0) {
        atomicMaxf(output, s_input[0]);
    }
}

__global__ void update_left_right(float* f, int nx, int ny)
{
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j > 0 && j < ny - 1)
    {
        f[j * nx] = f[(nx - 2) + j * nx];
        f[(nx - 1) + j * nx] = f[1 + j * nx];
    }
}

__global__ void update_top_bottom(float* f, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i > 0 && i < nx - 1)
    {
        f[i] = f[i + (ny - 2) * nx];
        f[i + nx * (ny - 1)] = f[i + nx];
    }
}
class FluidSolver2D
{
public:
    FluidSolver2D(float lx, float ly, int nx, int ny, float gamma, float output_dt)
        : lx_ { lx }, ly_ { ly }
        , gamma_ { gamma }, output_dt_ { output_dt }
    {
        dx_ = lx / (nx - 2.0f);
        dy_ = ly / (ny - 2.0f);
        nx_ = nx;
        ny_ = ny;
        n_ = static_cast<size_t>(nx_) * ny_;
        bytes_ = n_ * sizeof(float);
        rho_h.resize(nx * ny);
        printf("[FluidSolver2D] nx: %d, ny: %d, n: %zu, bytes: %zu\n", nx_, ny_, n_, bytes_);
        cudaMallocFluidSolver();
    }

    ~FluidSolver2D()
    {
        cudaFreeFluidSolver();
    }


    void init(auto const& rho, auto const& vx, auto const& vy, auto const& P)
    {
        printf("[init] rho.size(): %zu, vx.size(): %zu, vy.size(): %zu, P.size(): %zu\n", rho.size(), vx.size(), vy.size(), P.size());

        cudaMemcpy(rho_d, rho.data(), bytes_, cudaMemcpyHostToDevice);
        cudaMemcpy(vx_d, vx.data(), bytes_, cudaMemcpyHostToDevice);
        cudaMemcpy(vy_d, vy.data(), bytes_, cudaMemcpyHostToDevice);
        cudaMemcpy(P_d, P.data(), bytes_, cudaMemcpyHostToDevice);

        periodic_boundary(rho_d);
        periodic_boundary(vx_d);
        periodic_boundary(vy_d);
        periodic_boundary(P_d);

        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        primitive_to_conservative<<<grid_dim, block_dim>>>(rho_d, mass_d, vx_d, vy_d,
                                                           momx_d, momy_d, P_d, energy_d,
                                                           dx_, dy_, gamma_, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    float find_dt()
    {
        float C_cfl { 0.20f };
        float min_dxdy = std::min(dx_, dy_);

        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        compute_speed<<<grid_dim, block_dim>>>(rho_d, vx_d, vy_d, P_d, gamma_, nx_, ny_, speeds_d);
        cudaDeviceSynchronize();

        float max_speed = 0.0f;
        // float* max_speed_d;
        // cudaMalloc(&max_speed_d, sizeof(float));
        cudaMemcpy(max_speed_d, &max_speed, sizeof(float), cudaMemcpyHostToDevice);

        int blockSize = 128;
        int numBlocks = (n_ + blockSize * 2 - 1) / (blockSize * 2);
        maxReduceAtomic<<<numBlocks, blockSize, blockSize * sizeof(float)>>>(speeds_d, max_speed_d, n_, 1e-12f);
        cudaDeviceSynchronize();

        cudaMemcpy(&max_speed, max_speed_d, sizeof(float), cudaMemcpyDeviceToHost);

        // Free allocated memory
        // cudaFree(speeds_d);
        // cudaFree(max_speed_d);

        // Ensure max_speed is not zero to avoid division by zero
        max_speed = fmaxf(max_speed, 1e-6f);

        return C_cfl * min_dxdy / max_speed;
    }

    void solve(float t_start, float t_end)
    {
        float t { t_start };
        int n { 0 };
        float dt { 0.0f };

        while (t < t_end)
        {
            if (t >= output_dt_ * n)
            {
                output(n);
                n++;
            }

            cudaMemcpy(rho_h.data(), rho_d, bytes_, cudaMemcpyDeviceToHost);
            dt = find_dt();
            t += dt;

            step(dt);
        }
    }

    void step(float dt)
    {
        primitive_update(0.5f * dt);
        compute_fluxes();
        update_conserved(dt);
        conserved_to_primitive();
    }

    void periodic_boundary(float* field)
    {
        dim3 block_dim(256);
        dim3 grid_dim_x((nx_ + block_dim.x - 1) / block_dim.x);
        dim3 grid_dim_y((ny_ + block_dim.x - 1) / block_dim.x);
        update_left_right<<<grid_dim_y, block_dim>>>(field, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();
        update_top_bottom<<<grid_dim_x, block_dim>>>(field, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    void conserved_to_primitive()
    {
        dim3 block_dim(16, 16); // Define block dimensions explicitly
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);

        conserved_to_primitive_kernel<<<grid_dim, block_dim>>>(rho_d, mass_d, vx_d, vy_d,
                                                        momx_d, momy_d, P_d, energy_d,
                                                        dx_, dy_, gamma_, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    void primitive_update(float dt)
    {
        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        primitive_update_kernel<<<grid_dim, block_dim>>>(
            rho_d, vx_d, vy_d, P_d, rho_tmp_d, vx_tmp_d, vy_tmp_d, P_tmp_d,
            nx_, ny_, dx_, dy_, dt, gamma_);
        cudaGetLastError();
        cudaDeviceSynchronize();

        periodic_boundary(rho_tmp_d);
        periodic_boundary(vx_tmp_d);
        periodic_boundary(vy_tmp_d);
        periodic_boundary(P_tmp_d);

        cudaMemcpy(rho_d, rho_tmp_d, bytes_, cudaMemcpyDeviceToDevice);
        cudaMemcpy(vx_d, vx_tmp_d, bytes_, cudaMemcpyDeviceToDevice);
        cudaMemcpy(vy_d, vy_tmp_d, bytes_, cudaMemcpyDeviceToDevice);
        cudaMemcpy(P_d, P_tmp_d, bytes_, cudaMemcpyDeviceToDevice);
        cudaGetLastError();
        cudaDeviceSynchronize();
    }

    void extrapolate_to_interface()
    {
        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        extrapolate_to_interface_kernel<<<grid_dim, block_dim>>>(
            rho_d, vx_d, vy_d, P_d, rho_lx_d, rho_ly_d, rho_rx_d, rho_ry_d,
            vx_lx_d, vx_ly_d, vy_lx_d, vy_ly_d, vx_rx_d, vx_ry_d, vy_rx_d, vy_ry_d,
            P_lx_d, P_rx_d, P_ly_d, P_ry_d, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();

        periodic_boundary(rho_lx_d);
        periodic_boundary(rho_ly_d);
        periodic_boundary(rho_rx_d);
        periodic_boundary(rho_ry_d);
        periodic_boundary(vx_lx_d);
        periodic_boundary(vx_ly_d);
        periodic_boundary(vx_rx_d);
        periodic_boundary(vx_ry_d);
        periodic_boundary(vy_lx_d);
        periodic_boundary(vy_ly_d);
        periodic_boundary(vy_rx_d);
        periodic_boundary(vy_ry_d);
        periodic_boundary(P_lx_d);
        periodic_boundary(P_ly_d);
        periodic_boundary(P_rx_d);
        periodic_boundary(P_ry_d);
    }

    void compute_fluxes()
    {
        extrapolate_to_interface();

        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        compute_fluxes_kernel<<<grid_dim, block_dim>>>(
            mass_flux_x_d, momx_flux_x_d, momy_flux_x_d, energy_flux_x_d,
            mass_flux_y_d, momx_flux_y_d, momy_flux_y_d, energy_flux_y_d,
            rho_lx_d, rho_ly_d, rho_rx_d, rho_ry_d,
            vx_lx_d, vx_ly_d, vy_lx_d, vy_ly_d,
            vx_rx_d, vx_ry_d, vy_rx_d, vy_ry_d,
            P_lx_d, P_rx_d, P_ly_d, P_ry_d,
            gamma_, nx_, ny_);
        cudaGetLastError();
        cudaDeviceSynchronize();

        periodic_boundary(mass_flux_x_d);
        periodic_boundary(mass_flux_y_d);
        periodic_boundary(momx_flux_x_d);
        periodic_boundary(momx_flux_y_d);
        periodic_boundary(momy_flux_x_d);
        periodic_boundary(momy_flux_y_d);
        periodic_boundary(energy_flux_x_d);
        periodic_boundary(energy_flux_y_d);
    }

    void update_conserved(float dt)
    {
        dim3 block_dim(16, 16);
        dim3 grid_dim((nx_ + block_dim.x - 1) / block_dim.x,
                      (ny_ + block_dim.y - 1) / block_dim.y);
        update_conserved_kernel<<<grid_dim, block_dim>>>(
            mass_d, momx_d, momy_d, energy_d,
            mass_flux_x_d, momx_flux_x_d, momy_flux_x_d, energy_flux_x_d,
            mass_flux_y_d, momx_flux_y_d, momy_flux_y_d, energy_flux_y_d,
            nx_, ny_, dx_, dy_, dt);
        cudaGetLastError();
        cudaDeviceSynchronize();

        periodic_boundary(mass_d);
        periodic_boundary(momx_d);
        periodic_boundary(momy_d);
        periodic_boundary(energy_d);
    }

    void output(int n)
    {
        cudaDeviceSynchronize();
        cudaMemcpy(rho_h.data(), rho_d, bytes_, cudaMemcpyDeviceToHost);
        std::ofstream outfile("../data/output_rho_" + std::to_string(n) + ".csv");
        for (int j = 1; j < ny_ - 1; j++) {
            for (int i = 1; i < nx_ - 1; i++) {
                int idx = i + j * nx_;
                outfile << rho_h[idx];
                // std::cout << idx << ": " << rho_h[idx] << '\n';
                if (i != nx_ - 2)
                    outfile << ", ";
                else
                    outfile << std::endl;
            }
        }
        outfile.close();
    }

private:
    void cudaMallocFluidSolver()
    {
        // primitive variables
        cudaMalloc(&rho_d, bytes_);
        cudaMalloc(&vx_d, bytes_);
        cudaMalloc(&vy_d, bytes_);
        cudaMalloc(&P_d, bytes_);
        // temp array during primitive update
        cudaMalloc(&rho_tmp_d, bytes_);
        cudaMalloc(&vx_tmp_d, bytes_);
        cudaMalloc(&vy_tmp_d, bytes_);
        cudaMalloc(&P_tmp_d, bytes_);

        // conserved variables
        cudaMalloc(&mass_d, bytes_);
        cudaMalloc(&momx_d, bytes_);
        cudaMalloc(&momy_d, bytes_);
        cudaMalloc(&energy_d, bytes_);

        // fluxes for conserved variables
        cudaMalloc(&mass_flux_x_d, bytes_);
        cudaMalloc(&momx_flux_x_d, bytes_);
        cudaMalloc(&momy_flux_x_d, bytes_);
        cudaMalloc(&energy_flux_x_d, bytes_);
        cudaMalloc(&mass_flux_y_d, bytes_);
        cudaMalloc(&momx_flux_y_d, bytes_);
        cudaMalloc(&momy_flux_y_d, bytes_);
        cudaMalloc(&energy_flux_y_d, bytes_);

        // extrapolating to cell interface
        cudaMalloc(&rho_lx_d, bytes_);
        cudaMalloc(&rho_ly_d, bytes_);
        cudaMalloc(&rho_rx_d, bytes_);
        cudaMalloc(&rho_ry_d, bytes_);
        cudaMalloc(&vx_lx_d, bytes_);
        cudaMalloc(&vx_ly_d, bytes_);
        cudaMalloc(&vx_rx_d, bytes_);
        cudaMalloc(&vx_ry_d, bytes_);
        cudaMalloc(&vy_lx_d, bytes_);
        cudaMalloc(&vy_ly_d, bytes_);
        cudaMalloc(&vy_rx_d, bytes_);
        cudaMalloc(&vy_ry_d, bytes_);
        cudaMalloc(&P_lx_d, bytes_);
        cudaMalloc(&P_ly_d, bytes_);
        cudaMalloc(&P_rx_d, bytes_);
        cudaMalloc(&P_ry_d, bytes_);

        cudaMalloc(&speeds_d, (nx_ - 2) * (ny_ - 2) * sizeof(float));
        cudaMalloc(&max_speed_d, sizeof(float));
    }

    void cudaFreeFluidSolver()
    {
        // primitive variables
        cudaFree(rho_d);
        cudaFree(vx_d);
        cudaFree(vy_d);
        cudaFree(P_d);
        // temp array during primitive update
        cudaFree(rho_tmp_d);
        cudaFree(vx_tmp_d);
        cudaFree(vy_tmp_d);
        cudaFree(P_tmp_d);

        // conserved variables
        cudaFree(mass_d);
        cudaFree(momx_d);
        cudaFree(momy_d);
        cudaFree(energy_d);

        // fluxes for conserved variables
        cudaFree(mass_flux_x_d);
        cudaFree(momx_flux_x_d);
        cudaFree(momy_flux_x_d);
        cudaFree(energy_flux_x_d);
        cudaFree(mass_flux_y_d);
        cudaFree(momx_flux_y_d);
        cudaFree(momy_flux_y_d);
        cudaFree(energy_flux_y_d);

        // extrapolating to cell interface
        cudaFree(rho_lx_d);
        cudaFree(rho_ly_d);
        cudaFree(rho_rx_d);
        cudaFree(rho_ry_d);
        cudaFree(vx_lx_d);
        cudaFree(vx_ly_d);
        cudaFree(vx_rx_d);
        cudaFree(vx_ry_d);
        cudaFree(vy_lx_d);
        cudaFree(vy_ly_d);
        cudaFree(vy_rx_d);
        cudaFree(vy_ry_d);
        cudaFree(P_lx_d);
        cudaFree(P_ly_d);
        cudaFree(P_rx_d);
        cudaFree(P_ry_d);

        cudaFree(speeds_d);
        cudaFree(max_speed_d);
    }

    float lx_;
    float ly_;
    int nx_;            // padded dim
    int ny_;
    float dx_;
    float dy_;
    float gamma_;
    float output_dt_;
    size_t n_;          // total number of cells
    size_t bytes_;

    // primitive variables
    std::vector<float> rho_h;
    float* rho_d { nullptr };
    float* vx_d { nullptr };
    float* vy_d { nullptr };
    float* P_d { nullptr };
    // temp array during primitive update
    float* rho_tmp_d { nullptr };
    float* vx_tmp_d { nullptr };
    float* vy_tmp_d { nullptr };
    float* P_tmp_d { nullptr };

    // conserved variables
    float* mass_d { nullptr };
    float* momx_d { nullptr };
    float* momy_d { nullptr };
    float* energy_d { nullptr };

    // fluxes for conserved variables
    float* mass_flux_x_d { nullptr };
    float* momx_flux_x_d { nullptr };
    float* momy_flux_x_d { nullptr };
    float* energy_flux_x_d { nullptr };
    float* mass_flux_y_d { nullptr };
    float* momx_flux_y_d { nullptr };
    float* momy_flux_y_d { nullptr };
    float* energy_flux_y_d { nullptr };

    // extrapolating to cell interface
    float* rho_lx_d { nullptr };
    float* rho_ly_d { nullptr };
    float* rho_rx_d { nullptr };
    float* rho_ry_d { nullptr };
    float* vx_lx_d { nullptr };
    float* vx_ly_d { nullptr };
    float* vx_rx_d { nullptr };
    float* vx_ry_d { nullptr };
    float* vy_lx_d { nullptr };
    float* vy_ly_d { nullptr };
    float* vy_rx_d { nullptr };
    float* vy_ry_d { nullptr };
    float* P_lx_d { nullptr };
    float* P_ly_d { nullptr };
    float* P_rx_d { nullptr };
    float* P_ry_d { nullptr };

    // speeds for dt calculation
    float* speeds_d { nullptr };
    float* max_speed_d { nullptr };
};
