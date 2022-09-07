using Xfoil
using Printf
 
"""
some important functions from the example include:
- Xfoil.set_coordinates(x,y) = Creates the airfoil from the data that is given
- Xfoil.pane() = Breaks up the airofil into panels
- Xfoil.solve_alpha(......) = Solves for the coefficients that we desire: such as lift, drag, moment
- adding _cs to the end of the functions lets the solve be over complex steps
- Xfoil.alpha_sweep(.....) = simplifies the process of collecting the coefficient data by going through the process automatically

Some pseudo code:
- solve for coefficients normally:

Get the airfoil dimensions into lists that are seperated by x and y coordinates, and then use Xfoil.set_coordinates to have the program create the airfoil.

Then use Xfoil.pane to panelize the airfoil

Then define the values of alpha, Reynold's number, and mach. Alpha should be a set of numbers with equals steps between them.

Then initialize the outputs. The outputs need to have the same length as the list of alpha values.

Then run the Xfoil.solve_alpha function to get the values, which are then printed

- solve for the derivatives of the coefficients with respect to alpha:

Get the airfoil dimensions into lists that are seperated by x and y coordinates, and then use Xfoil.set_coordinates to have the program create the airfoil.

Then use Xfoil.pane to panelize the airfoil

Then define the values of alpha, Reynold's number, and mach. Alpha should be a set of numbers with equals steps between them.

Choose a very small step in alpha that will be used to find the change in the coefficients

Then initialize the outputs. The outputs need to have the same length as the list of alpha values.

Solve for two sets of coefficients using Xfoil.solve_alpha; one will be at alpha, and the other will be at alpha + the small step.

Using the two different sets of data, find the change in the them by finding the difference, and then divide this by the small increment that was chosen.
This value will need to be put into radians by multiply by 180/pi.

Print the results

-Using the complex step method:

The steps are similar to the previous two examples, but the addition of _cs needs to be added to the end of the functions

*note: I need to ask Adam what the complex step method really does because I kind of have an idea about it, but It isn't a strong one

-Using alpha_sweep:

The process is very similar; however, the set_coordinate and pane function are not needed. Also, the outputs do not need to be initialized since it is 
automatically done through alpha_sweep.
"""

#Get the airfoil data

#Get delimited files
x = [1,
0.996032,
0.983139,
0.970241,
0.957338,
0.944435,
0.931531,
0.918623,
0.905712,
0.892802,
0.879895,
0.866986,
0.854082,
0.841185,
0.828289,
0.815397,
0.802513,
0.789659,
0.776802,
0.76394,
0.751071,
0.738193,
0.725313,
0.712429,
0.699539,
0.68665,
0.67376,
0.660873,
0.647983,
0.637643,
0.632642,
0.630142,
0.627483,
0.622164,
0.609518,
0.59663,
0.58374,
0.570852,
0.557963,
0.545083,
0.532213,
0.519369,
0.506538,
0.493693,
0.480826,
0.467897,
0.454937,
0.441952,
0.428961,
0.415972,
0.402988,
0.390006,
0.377036,
0.36408,
0.351151,
0.338231,
0.325325,
0.312435,
0.299551,
0.286678,
0.273814,
0.260959,
0.248102,
0.235245,
0.222389,
0.209522,
0.196646,
0.183777,
0.170915,
0.158072,
0.14526,
0.132466,
0.1197,
0.107001,
0.094402,
0.082002,
0.069834,
0.057913,
0.046208,
0.03492,
0.024439,
0.0156,
0.009506,
0.005818,
0.003521,
0.002004,
0.000981,
0.000336,
0.000034,
0,
0.000034,
0.000336,
0.00098,
0.002004,
0.003521,
0.005818,
0.009505,
0.0156,
0.024439,
0.03492,
0.046207,
0.057912,
0.069834,
0.082002,
0.094401,
0.107001,
0.1197,
0.132465,
0.14526,
0.158071,
0.170915,
0.183776,
0.196646,
0.209522,
0.222389,
0.235244,
0.248101,
0.260958,
0.273814,
0.286677,
0.299551,
0.312434,
0.325324,
0.33823,
0.35115,
0.36408,
0.377036,
0.390005,
0.402987,
0.415972,
0.428961,
0.441952,
0.454937,
0.467897,
0.480825,
0.493692,
0.506538,
0.519368,
0.532213,
0.545083,
0.557963,
0.570852,
0.58374,
0.59663,
0.609518,
0.622164,
0.627483,
0.630142,
0.632642,
0.637643,
0.647983,
0.660873,
0.67376,
0.68665,
0.699539,
0.712429,
0.725313,
0.738193,
0.751071,
0.76394,
0.776802,
0.789659,
0.802513,
0.815397,
0.828289,
0.841185,
0.854082,
0.866985,
0.879894,
0.892802,
0.905712,
0.918623,
0.931529,
0.944435,
0.957338,
0.97024,
0.983139,
0.996032,
1]
y = [0.001801,
0.0021,
0.003111,
0.004063,
0.004962,
0.005811,
0.006614,
0.007376,
0.008099,
0.008784,
0.009434,
0.010052,
0.010638,
0.011195,
0.011711,
0.012164,
0.012552,
0.01287,
0.013115,
0.013276,
0.013358,
0.01336,
0.013282,
0.013127,
0.012901,
0.012602,
0.012238,
0.01181,
0.01132,
0.010884,
0.010659,
0.010545,
0.010715,
0.011053,
0.011814,
0.01254,
0.013219,
0.013855,
0.014451,
0.01502,
0.015559,
0.016076,
0.016567,
0.017041,
0.017496,
0.017936,
0.018362,
0.018773,
0.01917,
0.019553,
0.019921,
0.020273,
0.020613,
0.020941,
0.021253,
0.021556,
0.021848,
0.022127,
0.022396,
0.022654,
0.022897,
0.023125,
0.023337,
0.023531,
0.023703,
0.023851,
0.023965,
0.02404,
0.024068,
0.024037,
0.023935,
0.023751,
0.023467,
0.023063,
0.022517,
0.021807,
0.020906,
0.019777,
0.018367,
0.016611,
0.014453,
0.01199,
0.009653,
0.007724,
0.006112,
0.004665,
0.003303,
0.001952,
0.000626,
0.000012,
-0.000602,
-0.001928,
-0.003278,
-0.004641,
-0.006088,
-0.0077,
-0.009629,
-0.011967,
-0.01443,
-0.016588,
-0.018343,
-0.019754,
-0.020883,
-0.021785,
-0.022495,
-0.023042,
-0.023445,
-0.02373,
-0.023915,
-0.024016,
-0.024047,
-0.02402,
-0.023946,
-0.023831,
-0.023684,
-0.023512,
-0.023318,
-0.023107,
-0.022879,
-0.022636,
-0.022379,
-0.02211,
-0.021831,
-0.021539,
-0.021237,
-0.020925,
-0.020598,
-0.020258,
-0.019906,
-0.019538,
-0.019156,
-0.018758,
-0.018348,
-0.017923,
-0.017482,
-0.017028,
-0.016554,
-0.016063,
-0.015546,
-0.015008,
-0.014439,
-0.013843,
-0.013208,
-0.012529,
-0.011803,
-0.011043,
-0.010705,
-0.010535,
-0.010648,
-0.010874,
-0.01131,
-0.0118,
-0.012228,
-0.012593,
-0.012892,
-0.013118,
-0.013274,
-0.013352,
-0.013351,
-0.013269,
-0.013108,
-0.012864,
-0.012545,
-0.012158,
-0.011705,
-0.011189,
-0.010633,
-0.010047,
-0.00943,
-0.008779,
-0.008095,
-0.007372,
-0.006611,
-0.005808,
-0.004959,
-0.00406,
-0.003108,
-0.002098,
-0.001799]

Xfoil.set_coordinates(x,y)

#repanelize the coordinates
Xfoil.pane()

#set values
alpha = -10:1:10
re = 9.5e5
mach = 0.0

#initialize values
n_a = length(alpha)
c_l = zeros(n_a)
c_d = zeros(n_a)
c_dp = zeros(n_a)
c_m = zeros(n_a)
converged = zeros(Bool, n_a)
for i = 1:n_a
    c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter = 100)
end

#display results
println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
for i = 1:n_a
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end




