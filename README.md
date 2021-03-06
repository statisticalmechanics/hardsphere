# Hard Sphere

HSMD.c is the main code,
compiled by typing "make", using the "Makefile".
Parameters are read from the file "input".
The code is for molecular dynamics simulation of tenary (with three components A B C) hard spheres.
It uses collision dynamics, i.e. the system is updated till the next collision event, instead of with uniform time step dt.

In addition to constant volume V simulation, the code can also do volume compression. 
The compression consists of sequential constant V simulations, each for a time interval T. 
At the end of each T, the system is compressed affinely till the first collision event. 
Overall, this leads to an exponential volume decrease, and the compression rate is controlled by T.
The compression terminates untill the collision frequency, measured during each T, reaches a certain threshold, which implies reaching a certain packing fraction. 

Periodic boundary conditions are used.
Cell list is used to speed up the simulation and the cell size is scaled as volume shrinks.
The code calculate g(r) and output coordinate xyz file simutaneously 

Reference:

K. Zhang, W. W. Smith, M. Wang, Y. Liu, J. Schroers, M. D. Shattuck, C. S. O'Hern, Phys. Rev. E 90, 032311 (2014) 
