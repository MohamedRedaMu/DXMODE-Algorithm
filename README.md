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
├── run_DXMODE_main.m  
│
├── cec20_func.cpp
├── cec20_func.mex*
├── cec22_test_func.cpp
├── cec22_test_func.mex*
│
├── input_data2020/
├── input_data2022/
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

You can run the algorithm in three ways.

### Option 1: Run with default values
If no input arguments are passed, the algorithm uses its internal default settings.

```matlab
[goalReached, GlobalBest, countFE] = DXMODE_algorithm();
```

Default values inside `DXMODE_algorithm.m`:
```matlab
CECyear = 2020;
fNo     = 3;
nd      = 20;
lb      = -100;
ub      = 100;
```

### Option 2: Run by passing parameters directly
You can directly pass the benchmark settings to `DXMODE_algorithm.m`:

```matlab
[goalReached, GlobalBest, countFE] = DXMODE_algorithm(CECyear, fNo, nd, lb, ub);
```

Example for CEC2020:
```matlab
[goalReached, GlobalBest, countFE] = DXMODE_algorithm(2020, 3, 20, -100, 100);
```

Example for CEC2022:
```matlab
[goalReached, GlobalBest, countFE] = DXMODE_algorithm(2022, 5, 10, -100, 100);
```

### Option 3: Run using the main configuration file
You can also edit the settings in `run_DXMODE_main.m`, then run:

```matlab
run_DXMODE_main
```

---

## 🔧 Input Parameters

The function format is:

```matlab
DXMODE_algorithm(CECyear, fNo, nd, lb, ub)
```

### Parameters Description

- `CECyear` : Benchmark year  
  Supported values:
  - `2020`
  - `2022`

---

- `fNo` : Benchmark function number  

  - For **CEC2020**:
    - valid range: `1 – 10`

  - For **CEC2022**:
    - valid range: `1 – 12`

---

- `nd` : Problem dimension  

  Supported values in this implementation:
  - `10`
  - `20`

---

- `lb` : Lower bound  
  Typically:
  ```matlab
  -100
  ```

- `ub` : Upper bound  
  Typically:
  ```matlab
  100
  ```

---

## 📌 Important Notes

- Make sure that:
  - The correct CEC input data is available
  - The corresponding benchmark functions are compiled (MEX files)

---

## 🧪 Example Runs

### CEC2020 (Function 3, 20D)
```matlab
DXMODE_algorithm(2020, 3, 20, -100, 100);
```

### CEC2020 (Function 10, 10D)
```matlab
DXMODE_algorithm(2020, 10, 10, -100, 100);
```

### CEC2022 (Function 5, 20D)
```matlab
DXMODE_algorithm(2022, 5, 20, -100, 100);
```

### CEC2022 (Function 12, 10D)
```matlab
DXMODE_algorithm(2022, 12, 10, -100, 100);
```


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
