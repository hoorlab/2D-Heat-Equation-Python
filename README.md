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
