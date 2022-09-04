using plots
using Xfoil
 
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









