%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File written for the Spherical Shaping applied to single trajectories (Earth -> Mars/1989-ML/Tempel-1/Neptune)

% Guidelines: The Matlab files are divided in three folders and 4 files out in the main foder. This 4 files are the general functions, described below.
The 3 folders are: plot of the trajectory, functions that implement the spherical shaping algorithm and functions that implements the optimal trajectory computation.


%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Maple Files
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-dPrime: verification of D derivative



%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Matlab Files
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% General functions (Data Loading, Main File, Time Evolution Computation and Optimality Verification and Computation)

-LoadData: file that loads data on planets, measurements units. It is possible to add other bodies in the structure 'body'
-SpShTest: MAIN FILE that runs the spherical shaping on a combination Earth-(Mars,Tempel-1,1989ML,Neptune)
 (can generate the output file 'traj_rvut.mat' containing info on a trajectory if specified inside the file)
-timeEvolution: function that computes the time evolution of the trajectory
-verifyOptimality: function that first verifies the optimality of the obtained inverse solution and from that computes the optimal one

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% 'OptimalityFunctions' FOLDER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential equations used:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-adjointEq: adjoints differential equation for rk4 integrator (from minimum energy optimal control problem)
-adjointEqOde: adjoints differential equation for ode45 integrator


-diffEq: differential equations of dynamics of 2BP and adoints (minimum energy opt control prob)
-diffEqLinU: differential equations of dynamics of 2BP and adoints (minimum fuel opt control prob)
-diffEqLinUScaled: differential equations of dynamics of 2BP and adoints (minimum fuel opt control prob). Scaled version
-diffEqScaled: differential equations of dynamics of 2BP and adoints (minimum energy opt control prob). Scaled version
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
-transfer_bcs: boundary conditions for the transfer, to be used by Matlab bvp4c function that solves bvp (in the file OptimalTrajectoryIndirectTPBVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions that prepare the conversion from the spherical shaped trajectory to the optimal one.
The selection of which of these function to use is done by a menu in verifyOptimality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Computation of optimal trajectory with direct transcription method. Matrix version (the variables are stored in a matrix)
-OptimalTrajectoryDirectMat


%%%% Computation of optimal trajectory with direct transcription method, solved with OPTI Toolbox. Matrix version (the variables are stored in a matrix)
-OptimalTrajectoryDirectOPTI


%%%% Computation of optimal trajectory with direct transcription method. Vector version
-OptimalTrajectoryDirectVec


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum energy opt cont problem. Version with FSOLVE
-OptimalTrajectoryIndirectFSOLVEMultiple


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum fuel opt cont problem. Version with FSOLVE
-OptimalTrajectoryIndirectFSOLVEMultipleLinU


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum fuel opt cont problem. Version with FSOLVE, with scaled equations
-OptimalTrajectoryIndirectFSOLVEMultipleLinUScaled


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum energy opt cont problem. Version with FSOLVE, with scaled equations
-OptimalTrajectoryIndirectFSOLVEMultipleScaled


%%%% Computation of optimal trajectory with single shooting indirect method, for minimum energy opt cont problem. Version with FSOLVE
-OptimalTrajectoryIndirectFSOLVESingle


%%%% Computation of optimal trajectory with multiple shooting indirect method, for minimum energy opt cont problem. Version with OPTI
-OptimalTrajectoryIndirectOPTIMultiple


%%%% Computation of optimal trajectory with single shooting indirect method, for minimum energy opt cont problem. Version with Newton Raphson cycle
-OptimalTrajectoryIndirectSingle


%%%% Computation of optimal trajectory for minimum energy opt cont problem, version that uses Matlab bvp4c for solving the bvp
-OptimalTrajectoryIndirectTPBVP


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


%%%% Multiple Shooting function files for fsolve (version with linear control)
-ShootingFunctionFSOLVEMultipleLinU


%%%% Multiple Shooting function files for fsolve (version with linear control, scaled eq)
-ShootingFunctionFSOLVEMultipleLinUScaled


%%%% Multiple Shooting function files for fsolve (version with scaled eq)
-ShootingFunctionFSOLVEMultipleScaled


%%%% Single Shooting function files for fsolve (standard version)
-ShootingFunctionFSOLVESingle


%%%% Multiple Shooting function files for OPTI (standard version)
-ShootingFunctionOPTIMultipleConst
-ShootingFunctionOPTIMultipleObj





%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%'SphericalShapingFunctions' FOLDER. Containts the essential functions needed for the spherical shaping.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions called by the main functions of the shaping:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-baseFuncDeriv: compute derivatives of base functions of Spherical Shaping (needed for Spherical Shaping)
-cartesian2Spherical: convert cartesian positon into spherical coordinates (R,theta,Phi)
-KeplerElem2rv: convert Kepler elements to geocentric r,v 
-KeplerElem2rvMatrix: convert Kepler elements to geocentric r,v (matrix version)
-M2theta: convert M to theta
-rk4: runge-kutta 4 integrator
-timeFunc: compute derivatives of time (needed for Spherical Shaping)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main functions that realize the shaping and compute the trajectory and the control:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% The following 4 do the same thing but in different ways:
-SphericalShaping: SPHERICAL SHAPING FUNCTION WITH STEP ACCELERATION. Version with Newton Raphson cycle to find the unknown coefficient to solve satifying the boundary
conditions and the TOF constraint. Version with step acceleration algorithm (pay attention, is dependent upon the transfer); this was added because is some cases
the solver was stuck, decreasing with a very short step but not arriving to convergence
-SphericalShapingNoStepAcc: SPHERICAL SHAPING FUNCTION ORIGINAL VERSION. Version with Newton Raphson cycle to find the unknown coefficient to solve satifying the boundary
conditions and the TOF constraint. Original version.
-SphericalShapingWithFMincon: SPHERICAL SHAPING FUNCTION, USE FMINCON TO SOLVE THE CYCLE. Version that uses FMINCON to solve the nonlinear problem.
-SphericalShapingWithFSolve: SPHERICAL SHAPING FUNCTION, USE FSOLVE TO SOLVE THE CYCLE. Version that uses FSOLVE to solve the nonlinear problem. It is the best performing.



%%%% The following 2 are called by FMINCON or FSOLVE
-SSFunctionFMINCON: nonlinear function solved by FMINCON.
-SSFunctionFSOLVE: nonlinear function solved by FSOLVE.





%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%'TrajectoryPlotting' FOLDER. Contains the different functions that plot the trajectory obtained in the 3D space and the control.

-plotTrajectoryC: plot the trajectory obtained from the C code
-plotTrajectory: plot the trajectory and the control, using as input the spherical coordinates (R,theta,Phi) as given by the spherical shaping
-plotTrajectoryCartesianState: plot the trajectory and the control, using as input the cartersian coordinates
-plotTrajectoryMovie: plot a movie of the trajectory

%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------