using Xfoil
using Printf
using DelimitedFiles
using Plots
 
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

x = []
y = []

#This omitted section is used to find the airfoil data using the solve_alpha method
"""
Xfoil.set_coordinates(x,y)
Xfoil.pane()

#set values
alpha = -5:0.1:20

mach = 0.0
n_a = length(alpha)
#initialize values
re = 5e5
c_l = zeros(n_a)
c_d = zeros(n_a)
c_dp = zeros(n_a)
c_m = zeros(n_a)
converged = zeros(Bool, n_a)
for i = 1:n_a
    c_l[i], c_d[i], c_dp[i], c_m[i], converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter = 500)
end

#display results

    
println("Angle\t\tCl\t\tCd\t\tCm\t\tConverged")
for i = 1:n_a
    @printf("%8f\t%8f\t%8f\t%8f\t%d\n",alpha[i],c_l[i],c_d[i],c_m[i],converged[i])
end
"""

#This region is used to find the airfoil data using the alpha_sweep method
alpha = -5:0.1:20
re = 5e5
c_l, c_d, c_dp, c_m, converged = Xfoil.alpha_sweep(x, y, alpha, re, iter=500, zeroinit=false, printdata=true, npan = 200)

#finds lift to drag ratio from solved data
ratio = []
for i in 1:length(alpha)
    push!(ratio, c_l[i]/c_d[i])
end


#This omitted section is used to find the airfoil geometry for a NACA airfoil
"""
function NACA(m,p,t)
    global y_u = []
    global x_u = []
    global y_l = []
    global x_l = []
    for i in 0:.01:1
        yt = 5*t*(0.2969*sqrt(i)-0.1260*i-0.3516*(i)^(2)+0.2843*(i)^(3)-0.1015*(i)^(4))
        if 0<=i<=p
            yc=(m/(p^2))*(2*p*i-i^2)
        else p<i<=1
            yc=(m/(1-p)^2)*((1-2*p)+2*p*i-i^2)
        end
        if 0<=i<=p
            derivative = (2*m/(p^2))*(p-i)
        else p<i<=1
            derivative = (2*m/(1-p)^2)*(p-i)
        end
        theta = atan(derivative)
        push!(x_u, i-yt*sin(theta))
        push!(y_u, yc+yt*cos(theta))
        push!(x_l, i+yt*sin(theta))
        push!(y_l, yc-yt*cos(theta))
    end
end

"""
#derivative finder for airfoil data
Xfoil.set_coordinates(x,y)
Xfoil.pane(npan = 200)

#set values
alpha = -5:0.1:20

mach = 0.0
n_a = length(alpha)
#initialize values
re = 5e5
h = 1e-6


# initialize outputs
n_a = length(alpha)
c_l_a = zeros(n_a)
c_d_a = zeros(n_a)
c_dp_a = zeros(n_a)
c_m_a = zeros(n_a)
converged = zeros(Bool, n_a)

for i = 1:n_a
    c_l1, c_d1, c_dp1, c_m1, converged[i] = Xfoil.solve_alpha(alpha[i], re; mach, iter=500)
    c_l2, c_d2, c_dp2, c_m2, converged[i] = Xfoil.solve_alpha(alpha[i]+h, re; mach, iter=500)
    c_l_a[i] = (c_l2 - c_l1)/h * 180/pi
    c_d_a[i] = (c_d2 - c_d1)/h * 180/pi
    c_m_a[i] = (c_m2 - c_m1)/h * 180/pi
end
        



