%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File written for the Hill Shaping (uncomplete) and Pseudo Shaping

% NOTE: the Pseudo Shaping works properly, the trajectory recovered is a good one, but cannot be completely verified because the analytic functions for v,a are not
working fine and the numerical differentiation is not accurate (this is very strange because for other shapes it works).
% The Hill Shaping does't work because due to the imposed shape there can be H/G > 1 which causes a complex value of i (inclination).
% More time is required to find a good shaping approach (functions that don't cause this behavior)


%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Maple Files
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

VelAcc_Pseudo_Exp_Shape
VelAcc_Pseudo_Lin_Shape

Files that calculate velocity and acceleration taking into account the imposed shape of the equinoctial elements (exponential and linear).
They were checked and no mistakes were found (apparently). 


%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Matlab Files
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Main scripts

-Test_Shape_1: test on the shape. Integrate the trajectory with a known control profile and obtain the correspondent Hill (and Pseudo) elements.
-test_ver: test the verification algorithm. Generate a shape (a spiral) obtain its velocity,acc and control and then integrate dynamics with that control and check the shape obtained.

-Test_Hill_Shaping_1: test on the Hill Shaping. Doesn't work completely because due to the imposed shape there can be H/G > 1 which causes a complex value of i. More time needed to solve this.
-Test_Pseudo_Shaping_1: test on the Pseudo Shaping. Works but when using finite differences to obain the vel and acc the results are not very accurate.
                        Analytic formulas have been computed but don't produce required results (probably there is a mistake in computation somewhere)

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Auxiliary functions

-Cart2Kep: cartesian to Keplerian
-diffEqTangU: diff equation with tangential control acceleration of maximum value
-diffEqTangUDec: diff equation with tangential control acceleration of value decreasing with square of distance from Sun
-diffEqWithControl: diff equation with external control acceleration
-Hill2Kep, Pseudo2Kep: conversion from Hill (Pseudo) to Kepler elem
-Kep2Hill,Kep2Pseudo,Kep2rv: conversion from Keplerian to Hill,Pseudo,rv
-LoadData,M2theta,plotTrajectoryCartesianState: load the data, conver M to theta, plot trajectory
-rk4: rk4 integrator
-timeEvolution: compute the time evolution of the trajectory
-TrajectoryFromControl: verify the trajectory obtained by integration of control in the dynamics

%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Main functions

-Pseudo_Shaping_FMINCON: function that solves for the Pseudo Shaping as an optimization problem
-Pseudo_Shaping_Obj: objective function of above (Delta V)
-Pseudo_Shaping_Const: constraint function of above (TOF)
-Pseudo_Shaped: function that calculates the evolution of the elements given the angle vector and the coefficients

% Same thing of above but for Hill elements

-Hill_Shaping_FMINCON
-Hill_Shaping_Obj
-Hill_Shaping_Const
-Hill_Shaped
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------