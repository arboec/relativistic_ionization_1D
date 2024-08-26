# relativistic_ionization_1D
Supplimentary meterials for the article

This code solves a modified FW equation in scaling coordinates. The details of the algorithm can be found in a corresponding paper (submitted to computer physics communications). The code was run for a particular set of parameters using MATLAB R2022b version. The code is provided  as a supplimentary material to the article and is not intended to be a stand-alone program. You are free to use or modify the code, however, the author does not guarantee the correctness of the results.

Structure:
==========
run.m – the main script, which runs the computations

post_analysis.m – a routine to plot some graphs after run.m finishes.

==========

initializeGrid.m – initializes a uniform grid. Input parameters: (N = number of nods, L = length of the grid)

initializePGrid_pr.m – initializes a grid in momentum (Fourier) space, which corresponds to the grid in a coordinate space, provided by initializeGrid.m. Input parameters: (N = numer of nods, L = length of the corresponding grid in coordinate space)

initializeScaling.m – initializes scaling parameters, such as w and v. Input parameters: (SCALING = boolean: 0 == no scaling, 1 == scaling is ON).

InitializeSystemConstants.m – initializes initial conditions, such as number of nodes of the grid, size of the simulation box, delta t and so on.

InitializeWaveFunction.m – creates an initial wave function. Input parameters:(FILENAME_RE = input file for Re() part of wavefunction,FILENAME_IM = input file for Im() part of wavefunction , N = number of grid knots)

==========

m.m, mInv.m – a map from computational (uniform) coordinates to physical (non-uniform) ones.  Minv denotes an inverse transform == m^{-1}.

mder1.m, mder2.m – calculates the derivatives d m(x) / dx and d^2 m(x) / dx^2

m*L.m – corresponding versions of functions described above, but for an extended grid.
 
==========

R.m – scaling factor. Input parameters: (t = time, v = group velocity).

Rd.m – time derivative of R: d R(t) / dt.

==========

evolvePotential.m – calculates exp(-i*dt*V*0.5), i.e. propagates a wavefunction in time under the influence of potential.

calcV.m – calculates a central potential.

calcVAbsorb.m – calculates an absorbing potential on the boundaries.

calcVder2.m – calculates a second derivative of the central potential, d^2 V(x) / dx^2.

==========

evolveKinetic.m – calculates exp(-i*dt*K), i.e. propagates a wavefunction in time under the influence of the kinetic operator. Based on Krylov-Arnoldi scheme.

expmArnoldi.m – calculates a matrix exponent using Krylov-Arnoldi algorithm.

applyOperator.m – applays a kinetic operator K to the wave function, I.e. K*wave(x,t).

calculateCoefficients.m – calculates coefficients, i.e. values of T_0(x,t), T_1(x,t) T_2(x,t) and so on. Used for wavefunction propagation by applyOperator.m

darwinTerm_formula – calcualtes a contribution from Darwin-Spin-Orbital term. Used by applyOperator.m In practices, does not influence the results due to its small magnitude for the tested paramteres.

==========

calcA.m – calculates an electro-magnetic vector potential A.
calcExactSolutionA.m – calculates a phase of an exact free (without central potential) solution in form of \int_0^t sqrt(1+(p-A/c)^2/c^2)).

calcPhaseL.m – calculates a phase to be cancelled using the results, produced by calcExactSolutionA.m.

==========

calcProjections.m – calculates projections of the wave function to the spinors, which correspond to either positive or negative energies.

derivative_pr.m – calculates a derivative of a function provided on  a non-uniform mesh.

derivative_uniform_pr.m – calculates a derivative of a function provided on a uniform mesh.

expand.m – expands the grid by a few additional nodes on both ends of the mesh.

trunk.m – truncate the extended grid by removing a few nodes on both ends of the mesh.

timeOfNextRefining.m – determines a moment of time, when the mesh must be refined.

==========
