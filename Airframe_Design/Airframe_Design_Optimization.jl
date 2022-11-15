using SNOW
using VortexLattice
using LinearAlgebra

include("Optimization_Functions.jl")

#Initial conditions and bounds for the optimizer
x_0 = [0.41, 0.5, 0.15, 0.2, 0.31, 0.001, 0.3, 0.3, 1.8]

lx = [0.1, 0.5, 0.01, 0.2, 0.3, 0.0, 0.1, 0.1, 1.6]

ux = [0.5, 0.75, 0.2, 0.5, 0.5, 0.1, 0.5, 0.5, 2.0]

lg = [-Inf, -Inf, -Inf, -Inf]

ug = [0.0, 0.0, 0.0, 0.0]

ng = 4

options = Options(solver=IPOPT())


xopt, fopt, info = minimize(objective, x_0, ng, lx, ux, lg, ug)

println(xopt)

println(-1*fopt)

"""
[0.5000000099996722, 0.7500000099999041, 0.20000000999938128, 0.500000009999548, 0.2999999900001041, 0.000000099991631, 0.30000381469726345, 0.3000085830687774, 1.7]
124.249
"""