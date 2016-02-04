The code works as follows:

In the folder start are the necessary files needed to generate the samples

First, open up opt2.py.  Change the value for n and l depending on the number of
wavefunctions for a power of n and angular momentum of at max l.  the number
of functions it generates is n*(l+1)**2 since to count each m state for each possible
l state, we find we need 1 + 3 + 5 ... + 2l+1 = (l+1)**2 and we have n of each of these
functions.

We can also change the resolution and the accuracy of our groundstate.  This program looks
at the ground state of our helium with our wavefunctions using different nu values and
powers of nu (called space).  We alter these and see where the minimum is.  The error value
put in only finds the difference from the last minimum value, so the error represents
how much different our new min is.  Smaller error means we have the minimum of our
function, but it does not be we have the minimum energy.  Certain wf give better
minima.

Run opt2.py and it will make the .dat files needed to run the code with the minimized
wavefunctions

next, open up script2.py.  We can change the range of our intensities, the wavelength, and
the width of our laser.  After, we just run the script and it finds the groundstate prob
as time propogates.  After a long time, the laser is gone and the probability of being in the
groundstate is roughly equal to 1-prob(ionization).  The program finds this value by averaging
the probability over the last few timesteps and puts it into a file called knee_

All of the generated groundstates are kept in a folder with the knee graph.  This can be done
for different widths and frequencies for a single wf set.



There is a script for making graphs, but its cumbersome.  Will rewrite soon.



