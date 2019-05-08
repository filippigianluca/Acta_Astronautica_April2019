function [FunVal] = ShootingFunctionFSOLVESingle(lambda,inCond,fCond,mu,timeVec)

%ShootingFunctionFSOLVESingle
% Function that calculates difference between given and computed final
% state, which is the function to drive to zero and find the correspondent
% solution vector lambda
%
% INPUT
% lambda: vector of initial adjoints to be solved for from fsolve
% inCond: struct with initial conditions
% fCond: struct with final conditions
% mu: gravitational parameter
% timeVec: time vector (can be used for integration instead of t0,tf)
%
%
%OUTPUT
% FunVal: value of the function r(end)-rf and v(end)-vf to drive to zero in
% order to find the root (lambda, the initial value of the lagrange multipliers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the function of which to find the root

% Write down the initial and final conditions
r0 = inCond.r0;
v0 = inCond.v0;
t0 = inCond.t0;
rf = fCond.rf;
vf = fCond.vf;
tf = fCond.tf;

% Execute the integration in order to find the final state
%options_ODE_NUM = odeset('AbsTol',1e-10,'RelTol',1e-10);
%[tVec,Y] = ode45(@(t,x) diffEq(t,x,mu),[t0 tf],[r0 v0 lambda]);
[tVec,Y] = rk4(@diffEq,timeVec,[r0 v0 lambda],1,mu);
 Y = Y';


% Extract the state vectors
r = Y(:,1:3);
v = Y(:,4:6);
lambdaRes = Y(:,7:12);

% Evaluate the difference between the obtained cartesian state and the
% desired cartesian state in order to drive the function to zero and find
% the root (lambda)
FunVal(1:3) = abs((r(end,:) - rf)).^2;%./(rf - r0);
FunVal(4:6) = abs((v(end,:) - vf)).^2;%./(vf - v0);

end

