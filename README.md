# Simplex Method Solvers

This repository contains implementations of various simplex method algorithms for solving linear programming problems as well as integer programming problems. The following algorithms are included:

- **Dual Simplex Method**
- **Two-Phase Simplex Method**
- **Big M Simplex Method**
- **Gomory's Cutting Plane Method**

## Overview

Linear programming is a mathematical method for determining a way to achieve the best outcome in a given mathematical model. The simplex method is one of the most popular algorithms for solving linear programming problems.

### Algorithms Implemented

1. **Dual Simplex Method**
   - This method is used to solve linear programming problems where the primal is infeasible but the dual is feasible. The algorithm iteratively improves the dual solution while maintaining feasibility in the primal problem.

2. **Two-Phase Simplex Method**
   - This approach is used when the initial basic feasible solution is not readily available. The algorithm runs in two phases:
     - **Phase I:** Finds a feasible solution by minimizing the sum of artificial variables.
     - **Phase II:** Optimizes the objective function based on the feasible solution found in Phase I.

3. **Big M Simplex Method**
   - This method handles artificial variables in the simplex method by adding a large penalty (Big M) to the objective function. This ensures that the artificial variables are driven out of the solution as the algorithm progresses.

4. **Gomory's Cutting Plane Method**
  - This method is used to solve integer programming problem by first applying the usual simplex algorithms (Simplex and BigM method) and then using gomory's cutting plane method and dual simplex to get an integer solution.

## Getting Started

### Prerequisites

- C++ compiler (e.g., g++, clang++)
- Standard Template Library (STL)

### Compilation

To compile the individual files, you can use the following command in your terminal:

```bash
g++ -o <output_file_name> <source_file_name>.cpp
