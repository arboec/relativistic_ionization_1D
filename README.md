# A numerical solution of the 1D time-dependent Dirac equation
Supplimentary meterials for the article

This is a MATLAB implementation of the code, written by Heiko Bauke in Python in "Computational Strong-Field Quantum Dynamics", edited by Dieter Bauer, chapter 3 "Time-dependent relativistic wave equations: Numerics of the Dirac and the Klein-Gordon equation", listing 3, pages 106-107.

The code implements a basic split-operator approach for the solution of TDDE.

Structure:

======

run.m -- a main scrupt, lauches the computations.

=====

calcA.m -- function, calculates an electro-magnetic vector potential A.
calcV.m -- function, calculates a central potential for a given grid.
calcVAbsorb.m -- function, calculates an absorbing potential on the boundaries of the mesh.

=====
