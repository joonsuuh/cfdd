# 2D Computational Fluid Dynamics

# Report

[main.pdf](./latex/report/main.pdf)

# File Structure

```zsh
cfdd
├── README.md
├── data
│   ├── kelvin_helmholtz
│   │   ├── *.png
│   │   ├── *.csv
│   │   └── kelvin_helmholtz.md
│   └── sound_wave
│       ├── *.png
│       ├── *.csv
│       └── sound_wave.md
├── kh.ipynb
├── kh_growth.ipynb
├── kh_k100.mp4
├── kh_k12.mp4
├── kh_slow.ipynb
├── latex
│   └── report
│       ├── custom.sty
│       ├── draculatheme.sty
│       ├── images
│       │   └── *.png
│       ├── main.pdf
│       └── main.tex
├── sound_wave.ipynb
├── sound_wave_c2.mp4
└── src
    ├── fluid_solver_2d.h
    ├── fluid_solver_2d_kh.h
    ├── fluid_solver_2d_kh_MP.h
    ├── fs_2d_cmax.h
    ├── fs_2d_mp_fast.h
    ├── fs_2d_temp.h
    ├── kelvin_helmholtz.cpp
    ├── kelvin_helmholtz_MP.cpp
    ├── kh_mp_rand.cpp
    ├── kh_mp_rand_k12.cpp
    ├── kh_mp_slow.cpp
    ├── kh_mp_vy.cpp
    ├── kh_mp_vy.h
    ├── sound_wave.cpp
    ├── sw_cmax.cpp
    └── sw_cmax.cpp
```

# References 

- Chandrasekhar, S. 1961, International Series of Monographs on Physics, Oxford: Clarendon

- Landau, L., & Lifshitz, E. 2013, Fluid Mechanics: Volume 6 No. v. 6 (Elsevier Science)

- Rusanov, V. V. 1970, Journal of Computational Physics, 5, 507, doi: [10.1016/0021-9991(70)
90077-X](http://doi.org/10.1016/0021-9991(70)90077-X)

- Thorne, K. S., & Blandford, R. D. 2017, Modern classical physics: optics, fluids, plasmas, elasticity,
relativity, and statistical physics (Princeton: Princeton University Press)