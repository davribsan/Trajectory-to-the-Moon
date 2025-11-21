# Interplanetary Trajectory to the Moon

This repository contains Python scripts for designing and analyzing interplanetary trajectories from a Low Earth Orbit (LEO) to lunar orbits, minimizing the Δv budget. 

Specifically, it includes:

- **Izzo’s Lambert solver** — used to solve Lambert’s problem for transfers between two position vectors over a specified time of flight.
- **Orbital mechanics utilities** — functions for orbital parameters, anomalies, and trajectory computations.
- **Trajectory propagation** — numerical integration of the two-body problem using `solve_ivp`.
- **Visualization tools** — scripts for plotting trajectories and orbits.

## Features
- Modular Python scripts organized in folders:  
  - `solvers/` — Lambert solver implementation 
  - `utils/` — orbital and geometric functions  
  - `data/` — orbital constants and mission-specific inputs 

## Requirements

### Core Software
- Python 3.12.4 (recommended)

### Mandatory Python Packages
- numpy
- scipy
- matplotlib

### Internal Modules
- `solvers/izzo.py`
- `solvers/linalg.py`  
- `utils/orbital.py`  
- `utils/geometry.py`  
- `data/orbital_mechanics_constants.py`  
- `data/inputs.py`  

## Usage
1. Adjust mission parameters in `data/inputs.py`.
2. Run the main script:  
   - `main.py`
3. The results will not be shown on screen, however they will be saved locally in the same folder under the name **Trajectory_delta_v.txt**.
