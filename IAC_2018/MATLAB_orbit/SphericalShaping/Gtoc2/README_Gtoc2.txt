%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File written for the Spherical Shaping applied to the GTOC2 competition (to the first 12 combinations)

% Guidelines: The Matlab files are divided in three folders and 4 files out in the main foder. This 4 files are the general functions, described below.
The 3 folders are: plot of the trajectory, functions that implement the spherical shaping algorithm and functions that implements the optimal trajectory computation.



%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Matlab Files
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% General functions (Data Loading, Main File, Time Evolution Computation and Optimality Verification and Computation)

-AstComb: file the combinations of asteroids and the departure and arrival date vectors for each transfer
 (the file ast_ephem.txt contains the ephemeris of all the asteroids of GTOC2) 
-readFileImportdata: function for reading external files
-Gtoc2_test: MAIN FILE that runs the spherical shaping on a combination of GTOC2
-timeEvolution: function that computes the time evolution of the trajectory
-verifyOptimalityGTOC2PhasesSeparated: function that first verifies the optimality of the obtained inverse solution and from that computes the optimal one, for each phase separately

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% 'OptimalityFunctions' FOLDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential equations used:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-adjointEq: adjoints differential equation for rk4 integrator (from minimum energy optimal control problem)
-adjointEqOde: adjoints differential equation for ode45 integrator


-diffEq: differential equations of dynamics of 2BP and adoints (minimum energy opt control prob)
-diffEqTangU: differential equation with a tangential constant trhust
-diffEqWithControl: differential equation with a predetermined control history along all the integration time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Others:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-comparisonLRVU: comparsion of input R,V,U from 2 different methods (and eventually of Lagrange Multipliers)
-euler: euler integrator
-rk4: rk4 integrator
-scalingFunction: scaling function
-TrajectoryFromControl: obtain the trajectory (r,v) from a known control profile in cartesian coordinates and inital conditions (r0,v0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions that prepare the conversion from the spherical shaped trajectory to the optimal one.
The selection of which of these function to use is done by a menu in verifyOptimality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Computation of optimal trajectory with direct transcription method. Matrix version (the variables are stored in a matrix)
-OptimalTrajectoryDirectMat


%%%% Computation of optimal trajectory with direct transcription method, solved with OPTI Toolbox. Matrix version (the variables are stored in a matrix)
-OptimalTrajectoryDirectOPTI


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum energy opt cont problem. Version with FSOLVE
-OptimalTrajectoryIndirectFSOLVEMultiple


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint and objective functions for the different schemes of direct transcription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Constraints and objective for Gauss Lobatto Transcription:
-DirTranscFunctionFMINCONGLConst
-DirTranscFunctionFMINCONGLObj


%%%% Constraints and objective for Collocation Transcription, implemented with different schemes:
euler,runge-kutta,hermite-simpson (recommended) and crank-nicholson. Matrix form for the variables
-DirTranscFunctionFMINCONMatConst
-DirTranscFunctionFMINCONMatObj


%%%% Constraints and objective for Collocation Transcription, implemented with different schemes:
euler,runge-kutta,hermite-simpson (recommended) and crank-nicholson. Vector form for the variables
-DirTranscFunctionFMINCONVecConst
-DirTranscFunctionFMINCONVecObj


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint and objective functions for the different schemes of indirect methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% Multiple Shooting function files for fmincon (standard version)
-ShootingFunctionFMINCONConst
-ShootingFunctionFMINCONObj


%%%% Multiple Shooting function files for fsolve (standard version)
-ShootingFunctionFSOLVEMultiple


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%'SphericalShapingFunctions' FOLDER. Containts the essential functions needed for the spherical shaping.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions called by the main functions of the shaping:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-addExcessVelocity: function that adds an excess velocity to the spacecraft at departure
-baseFuncDeriv: compute derivatives of base functions of Spherical Shaping (needed for Spherical Shaping)
-cartesian2Spherical: convert cartesian positon into spherical coordinates (R,theta,Phi)
-KeplerElem2rv: convert Kepler elements to geocentric r,v 
-KeplerElem2rvMatrix: convert Kepler elements to geocentric r,v (matrix version)
-M2theta: convert M to theta
-timeFunc: compute derivatives of time (needed for Spherical Shaping)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main functions that realize the shaping and compute the trajectory and the control:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% The following 4 do the same thing but in different ways:
-SphericalShaping: SPHERICAL SHAPING FUNCTION WITH STEP ACCELERATION. Version with Newton Raphson cycle to find the unknown coefficient to solve satifying the boundary
conditions and the TOF constraint. Version with step acceleration algorithm (pay attention, is dependent upon the transfer); this was added because is some cases
the solver was stuck, decreasing with a very short step but not arriving to convergence
-SphericalShapingWithFSolve: SPHERICAL SHAPING FUNCTION, USE FSOLVE TO SOLVE THE CYCLE. Version that uses FSOLVE to solve the nonlinear problem. It is the best performing.


%%%% The following is called by FSOLVE
-SSFunctionFSOLVE: nonlinear function solved by FSOLVE.



%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%'TrajectoryPlotting' FOLDER. Contains the different functions that plot the trajectory obtained in the 3D space and the control.

-M2theta: convert M to theta
-plotTrajectoryC: plot the trajectory obtained from the C code
-plotTrajectoryMultBodiesCartContMultTraj: plot the trajectory and the control, using as input the cartersian coordinates, considering that the total trajectory is stored
in a multidimensional matrix
-plotTrajectoryMultBodiesSpherNoCont: plot the trajectory without control, using as input the cartersian coordinates, considering that the total trajectory is stored
in a multidimensional matrix
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------