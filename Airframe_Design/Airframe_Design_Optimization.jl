using SNOW
using VortexLattice
using LinearAlgebra

include("Optimization_Functions.jl")

#Initial conditions and bounds for the optimizer
x_0 = [0.41, 0.5, 0.15, 0.2, 0.31, 0.001, 0.3, 0.3, 1.8]

lx = [0.4, 0.5, 0.1, 0.2, 0.3, 0.0, 0.2, 0.2, 1.6]

ux = [0.5, 0.75, 0.2, 0.4, 0.5, 0.1, 0.4, 0.4, 2.0]

lg = [-Inf, -Inf, -Inf, -Inf]

ug = [0.0, 0.0, 0.0, 0.0]

ng = 4

options = Options(solver=IPOPT())


xopt, fopt, info = minimize(objective, x_0, ng, lx, ux, lg, ug)

println(xopt)

println(-1*fopt)

"""
[0.5000000099993505, 0.7500000099998233, 0.5000000099996648, 0.4000000099991854, 0.29999999000018374, -9.998559696625576e-9, 0.3, 0.29999809265135313, 1.5999999840682646]
97.82047470461352
"""