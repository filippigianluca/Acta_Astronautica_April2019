function [r,v,lambda,tVec] = OptimalTrajectoryIndirectFSOLVESingle(timeVec,t0,tf,r0,rf,v0,vf,lambda0,mu)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method

% Method: indirect single shooting with FSOLVE

% INPUT
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% lambda0: initial lagrange multipliers
% M: number of grid points
% nT: number of integration points in a phase
% mu: gravitational parameter
%
% OUTPUT
% r: matrix with position
% v: matrix with velocity
% lambda: matrix with lagrange multipliers
% tVec: vector with time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define errors and tolerances for the shooting method

% Write initial and final conditions
inCond.r0 = r0;
inCond.v0 = v0;
inCond.t0 = t0;
fCond.rf = rf;
fCond.vf = vf;
fCond.tf = tf;

% Find the correct initial conditions for the lagrange multipliers from
% fsolve
optionsFS = optimoptions('fsolve','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',100000,'MaxIter',100000);
[lambdaRes,fval] = fsolve(@(lambda) ShootingFunctionFSOLVESingle(lambda,inCond,fCond,mu,timeVec),lambda0,optionsFS);

% Integrate to obtain the solution
%options_ODE_NUM = odeset('AbsTol',1e-10,'RelTol',1e-10);
%[tVec,Y] = ode45(@(t,x) diffEq(t,x,mu),[t0 tf],[r0 v0 lambdaRes]);
[tVec,Y] = rk4(@diffEq,timeVec,[r0 v0 lambdaRes],1,mu);
 Y = Y';

% Extract the state vectors
r = Y(:,1:3);
v = Y(:,4:6);
lambda = Y(:,7:12);
   
% Print the difference between obtained final values and desired final
% values
fprintf('The obtained final position is: %f  \n',r(end,:)); 
fprintf('The desired final position is: %f  \n',rf); 
fprintf('The obtained final velocity is: %f  \n',v(end,:)); 
fprintf('The desired final velocity is: %f  \n',vf); 

  
end

