The file hfk_class_graph2.cpp contains the modified Baldwin-Gillam code for grid homology for links in S^3. 
It adds a second edge set representing the map tau on an equivariant grid diagram 
for a freely 2-periodic knot, and then applies Gaussian elimination on the second edge set to compute
the action of the tau on grid homology. 

The file hfkdemo_graph2.cpp allows the user to input different equivariant grid diagrams
and prints the tau edges on grid homology along with the ranks of the quotient of grid homology
by 1+tau.

12n403allgradings.pdf shows the output of this code on the 10 by 10 grid diagram for 
the freely 2-periodic knot 12n403.

rp3_12n403_quotient_ALL_A.cpp is based on Daniele Celoria's grid homology for lens
spaces code, which is originally written in C++ and computes grid homoogy 
with Z coefficients. This file is written in C++ and only computes with F2 coefficients
which allows for computations on much larger grids.

rp3_12n403_quotient_ALL_A.pdf contains the output of this code.
