# Algorithms for SEM-EDS Mineral Dust Classification
[![status](https://joss.theoj.org/papers/c2564d4c44b4ee77c24ac32f7431a6b2/status.svg)](https://joss.theoj.org/papers/c2564d4c44b4ee77c24ac32f7431a6b2)

## About

This repository hosts the Julia-equivalent functions from the original [`eds-classification`](https://github.com/weber1158/eds-classification) repository. The Julia functions have been migrated here so that Julia users do not need to sift through the 80+ files in the original repository (which is mostly for MATLAB users).

If you want the MATLAB version of this repository, please see the File Exchange submission here: [![View my project on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/170771)


## Documentation
See the Jupyter Notebook or HTML files in the `docs` folder for tutorials and documentation.

## Installation
1. Download the GitHub repository and save the algorithms somewhere on your PC.
2. Use the `include` function to import the algorithms within your Julia scripts. For example: `include("donarummo_classification.jl")`. Note that the algorithm files will need to be on the search path for this to work.


## Disclaimer
Please note that `eds-classication.jl` is in development as an alternative to the original MATLAB repository, and as such it is not as robust as its counterpart. There are several functions that still need to be translated from MATLAB into Julia code. In addition, the machine learning algorithm in the Julia repository should not be considered homologous to the MATLAB version, but it is similar.


## How to cite
<a href="https://joss.theoj.org/papers/c2564d4c44b4ee77c24ac32f7431a6b2">
  <img src="https://joss.theoj.org/papers/c2564d4c44b4ee77c24ac32f7431a6b2/status.svg" width="200" height="26" alt="status">
</a>

Please use the information below for citing the software:

#### APA-like
Weber, Austin M., (2025). Algorithms for SEM-EDS mineral dust classification. _Journal of Open Source Software_, *10*(107), 7533, https://doi.org/10.21105/joss.07533

#### `BibTeX`:
```tex
@article{weber2025,
    author = {Weber, Austin M.},
    title = {Algorithms for {SEM-EDS} mineral dust classification},
    journal = {Journal of Open Source Software},
    volume = {10},
    number = {107},
    pages = {7533},
    year = {2025},
    DOI = {10.21105/joss.07533}
}
```
