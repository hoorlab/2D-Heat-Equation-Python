# 2D Heat Equation Solver in Python

Solves the 2D heat equation ∂u/∂t = α(∂²u/∂x² + ∂²u/∂y²) using the Finite Difference Method.

## Features
- Crank-Nicolson implicit scheme for stability
- Dirichlet and Neumann boundary conditions  
- NumPy vectorized implementation
- Animated visualization of temperature evolution

## Run it
```bash
pip install -r requirements.txt
python heat_solver.py
```

## Method
Uses central difference in space, Crank-Nicolson in time. Stable for all time steps.

## Author
Malika Hurain – MSc student applying to UoA PhD in Optimization/ML
