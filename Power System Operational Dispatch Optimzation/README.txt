This project is a power system operational dispatch optimalization program using 2-step Robust Optimization using YALMIP and Gurobi optimization module. The motive is to optimize the operations of a modern power system network consisting of generation and demand uncertainties coming from variable loads and PV power plant.

Code explanations:
Primal_Problem : 1st stage robust optimization to find dispatch conditions/numbers with the minimum operational expenditure/cost.
Dual_Problem : 2nd stage robust optimization to include operational uncertainties (PV power plant and variable loads) and thus finding the operational dispatch for minimum operational cost in the worst case scenario of the uncertainties (lowest PV generation with highest variable load condition)