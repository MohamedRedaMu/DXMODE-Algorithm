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
```

## 📌 Notes

Benchmark settings are defined inside the main file.

Modify:
- dimension
- function number
- bounds

Requires CEC benchmark files.

---

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
---

## 📜 License

This project is released under the MIT License. See the LICENSE file for details.

---


## 📧 Contact

**Dr. Mohamed Reda**  
University of Central Lancashire, UK  
Mansoura University, EGY

- 📩 Personal: [mohamed.reda.mu@gmail.com](mailto:mohamed.reda.mu@gmail.com)  
- 📩 Academic: [mohamed.reda@mans.edu.eg](mailto:mohamed.reda@mans.edu.eg)  

---

## 🌐 Academic Profiles

- 🧑‍🔬 ORCID: https://orcid.org/0000-0002-6865-1315  
- 🎓 Google Scholar: https://scholar.google.com/citations?user=JmuB2qwAAAAJ  
- 📊 Scopus: https://www.scopus.com/authid/detail.uri?authorId=57220204540  
- 📚 Web of Science: https://www.webofscience.com/wos/author/record/3164983  
- 🧾 SciProfiles: https://sciprofiles.com/profile/Mreda  

---

## 🔗 Professional & Social Links

- 💼 LinkedIn: https://www.linkedin.com/in/mraf  
- 🔬 ResearchGate: https://www.researchgate.net/profile/Mohamed-Reda-8  
- 🎓 Academia: https://mansoura.academia.edu/MohamedRedaAboelfotohMohamed  
- 📘 SciLit: https://www.scilit.net/scholars/12099081  
- 🧮 MATLAB Central: https://uk.mathworks.com/matlabcentral/profile/authors/36082525  
- ▶️ YouTube: https://youtube.com/@mredacs  

---

## Acknowledgement

This repository accompanies the published DXMODE paper in Information Sciences.
