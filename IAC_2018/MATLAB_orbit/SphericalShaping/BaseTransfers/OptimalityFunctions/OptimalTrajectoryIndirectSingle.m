function [r,v,lambda,tVec] = OptimalTrajectoryIndirectSingle(timeVec,t0,tf,r0,rf,v0,vf,lambda0,mu)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method

% Method: indirect multiple shooting with manual Newton Raphson cycle

% INPUT
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% lambdAd: lagrange multipliers
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

% Relative tolerance for error on final position and velocity (norms)
toll = 1e-2;
err = toll + 1;
cont = 0;
maxIter = 100;

% Define lambda initial and fact (percentage modification factor)
lambdaIn = lambda0;
fact = 0.001;
h = 1e-14;

% Start the cycle
disp('Entering the optimal control cycle');
while ((err > toll) && (cont < maxIter)) 
   
   % Integrate
   options_ODE_NUM = odeset('AbsTol',1e-10,'RelTol',1e-10); 
   [tVec,Y] = ode45(@(t,x) diffEq(t,x,mu),[t0 tf],[r0 v0 lambdaIn],options_ODE_NUM);
   %[tVec,Y] = rk4(@diffEq,timeVec,[rf vf lambdaIn],0,mu);
   %Y = Y';
   % Extract the state vectors
   r = Y(:,1:3);
   v = Y(:,4:6);
   lambda = Y(:,7:12);
   
   % Calculate the errors on final position and velocity
   err_rf = (rf-r(end,:))./rf;
   err_vf = (vf-v(end,:))./vf;
   err = max([norm(err_rf) norm(err_vf)]);
%    err_r0 = (r0-r(1,:))./r0;
%    err_v0 = (v0-v(1,:))./v0;
%    err = max([norm(err_r0) norm(err_v0)]);
   
   % Modify the lambdaIn
   % Modification based on random values
   %lambdaIn = lambdaIn + rand(1,6).*lambdaIn*fact;
   
   % Modification based on Newton Raphson cycle
   [dum,Yp] = ode45(@(t,x) diffEq(t,x,mu),[t0 tf],[r0 v0 lambdaIn+h],options_ODE_NUM);
   [dum,Ym] = ode45(@(t,x) diffEq(t,x,mu),[t0 tf],[r0 v0 lambdaIn-h],options_ODE_NUM);
   rp = Yp(:,1:3);
   vp = Yp(:,4:6);
   rm = Ym(:,1:3);
   vm = Ym(:,4:6);
   err_rfp = (rf-rp(end,:))./rf;
   err_vfp = (vf-vp(end,:))./vf;
   err_rfm = (rf-rm(end,:))./rf;
   err_vfm = (vf-vm(end,:))./vf;
   derr_r = (err_rfp - err_rfm)/(2*h);
   derr_v = (err_vfp - err_vfm)/(2*h);
   lambdaIn = lambdaIn - [err_rf err_vf]./[derr_r derr_v];
   
%    [dum,Yp] = rk4(@diffEq,timeVec,[rf vf lambdaIn+h],0,mu);
%    [dum,Ym] = rk4(@diffEq,timeVec,[rf vf lambdaIn-h],0,mu);
%    Yp = Yp';
%    Ym = Ym';
%    rp = Yp(:,1:3);
%    vp = Yp(:,4:6);
%    rm = Ym(:,1:3);
%    vm = Ym(:,4:6);
%    err_r0p = (r0-rp(1,:))./r0;
%    err_v0p = (v0-vp(1,:))./v0;
%    err_r0m = (r0-rm(1,:))./r0;
%    err_v0m = (v0-vm(1,:))./v0;
%    derr_r = (err_r0p - err_r0m)/(2*h);
%    derr_v = (err_v0p - err_v0m)/(2*h);
%    lambdaIn = lambdaIn + [err_r0 err_v0]./[derr_r derr_v];


   % Update cont
   cont = cont+1;
   fprintf('Cont %f  \n', cont);
end

if (cont == maxIter)
 disp('ERROR - Cycle for optimal control, direct shooting, not converged');
else
 disp('OK - Cycle for optimal control, direct shooting, has converged');
end

end

