# DXMODE-Algorithm
A Dynamic Explorative Multi-Operator Differential Evolution Algorithm for Engineering Optimization Problems

Official MATLAB implementation of the DXMODE algorithm.

## 📄 Paper

Mohamed Reda, Ahmed Onsy, Amira Y. Haikal, and Ali Ghanbari  
**DXMODE: A dynamic explorative multi-operator differential evolution algorithm for engineering optimization problems**  
Information Sciences, 2025  
DOI: https://doi.org/10.1016/j.ins.2025.122271

---

## Overview

DXMODE is a Differential Evolution (DE) algorithm designed to improve:

- Exploration capability
- Population diversity
- Convergence performance

Main contributions:
- Error-based Linear Population Decay (ELPD)
- Multi-operator mutation strategies
- Chaotic exploration (MNCE)
- Gaussian exploration (AGE)
- Adaptive archive
- Hybrid local search

---

## 📂 Repository Structure

```text
DXMODE-Algorithm/
│
├── DXMODE_algorithm.m
├── cost_cec2020.m
├── cost_cec2022.m
│
├── cec_functions/
│   ├── cec20_func.cpp
│   ├── cec22_test_func.cpp
│   └── (compiled mex files)
│
├── input_data/
│   ├── cec2020/
│   └── cec2022/
│
├── examples/
│   └── run_DXMODE_demo.m
│
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore
```


---

## ⚙️ Requirements

- MATLAB (recommended R2023a or later)
- Optimization Toolbox (`fmincon`)
- Statistics Toolbox (`lhsdesign`, `normrnd`)

---

## ▶️ How to Run

```matlab
[goalReached, GlobalBest, countFE] = DXMODE_algorithm();


## 📌 Notes

Benchmark settings are defined inside the main file.

Modify:
- dimension
- function number
- bounds

Requires CEC benchmark files.

## 📚 Citation

If you use this code, please cite:

```bibtex
@article{reda2025dxmode,
  title   = {DXMODE: A dynamic explorative multi-operator differential evolution algorithm for engineering optimization problems},
  author  = {Reda, Mohamed and Onsy, Ahmed and Haikal, Amira Y. and Ghanbari, Ali},
  journal = {Information Sciences},
  volume  = {717},
  pages   = {122271},
  year    = {2025},
  doi     = {10.1016/j.ins.2025.122271}
}
```

## 📜 License

This project is released under the MIT License. See the LICENSE file for details.

## 📧 Contact

Mohamed Reda  
Emails: 
  mohamed.reda.mu@gmail.com
  mohamed.reda@mans.edu.eg


## Acknowledgement

This repository accompanies the published DXMODE paper in Information Sciences.
