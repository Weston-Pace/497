using SNOW
using VortexLattice
using LinearAlgebra

include("Optimization_Functions.jl")

#Initial conditions and bounds for the optimizer
x_0 = [0.41, 0.5, 0.15, 0.2, 0.31, 0.001, 0.3, 0.3, 2.0]

lx = [0.1, 0.5, 0.01, 0.2, 0.3, 0.0, 0.1, 0.1, 1.7]

ux = [0.5, 0.75, 0.2, 0.5, 0.5, 0.05, 0.5, 0.5, 2.5]

lg = [-Inf, -Inf, -Inf, -Inf]

ug = [0.0, 0.0, 0.0, 0.0]

ng = 4

options = Options(solver=IPOPT())


xopt, fopt, info = minimize(objective, x_0, ng, lx, ux, lg, ug)

println(xopt)

println(-1*fopt)

"""
[0.5000000099706225, 0.7500000099904578, 0.20000000996657682, 0.5000000099605839, 0.2999999900122995, -9.929404636797622e-9, 0.29999427309849536, 0.30000476836969897, 1.6999999894261193]
124.24983405856212
"""