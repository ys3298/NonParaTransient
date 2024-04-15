# Non-Parametric Analysis of Transient Data: a Pseudo-Competing Event Approach

This repository contains simulation code used in the paper titled "Non-Parametric Analysis of Transient Data: a Pseudo-Competing Event Approach".

## Overview

The code provided here accompanies the research paper and facilitates the reproduction of the results. It comprises simulation scripts written in R, designed to emulate the scenarios described in the paper.


## Usage

The algorithm are in the 'algorithm' folder:

- calibrateAUC.R: calculate AUC

- CI_functions.R: calculate bootstrap CI of AUC

- generate_data.R: generate data for simulations

The main simulations (section 3.2.2) are conducted with simulation_para.R.

The MB test simulations (section 3.2.1) are conducted with test_MB.R (MB test) and sim_MB_test.R (proposed method).

The sensitivity analysis (section 4) are conducted with simulation_sensitivity.R.

Tables and figures are generated using functions in 'visualize results' folder.

Note: we scaled the time from t=0.15 year to month by multiplying the AUC by 30 and devided $\lambda$ by 30. 