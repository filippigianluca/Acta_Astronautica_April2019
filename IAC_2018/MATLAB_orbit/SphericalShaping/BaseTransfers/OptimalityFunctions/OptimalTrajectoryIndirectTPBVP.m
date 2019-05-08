function [r,v,lambda,tVec] = OptimalTrajectoryIndirectTPBVP(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,n,mu)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method


% Method: tpbvp solver of Matlab

% INPUT
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% lambdAd: lagrange multipliers
% n: number of grid points
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
%% Code

%--------------------------------------------------------------------------
% Define parameters
%--------------------------------------------------------------------------

% Time vector of the mesh
tVec = linspace(t0,tf,n);

% An interpolation is needed in order to have the values of the variables
% at each moment in time corresponding to the mesh vector
 for i=1:3
  rInt(i,:) = spline(timeVec,R(i,:),tVec);
  vInt(i,:) = spline(timeVec,Vcart(i,:),tVec);
  lambdarInt(i,:) = spline(timeVec,lambdaAd(i,:),tVec);
  lambdavInt(i,:) = spline(timeVec,lambdaAd(i+3,:),tVec);
 end

%  %Add some 'noise' to the initial guess data (a randomatic variation given
% %as a percentual of the original value)
% a = -0.001;
% b = 0.001;
%  for i=1:3
%   rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,n)).*rInt(i,:);
%   vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,n)).*vInt(i,:);
%   lambdarInt(i,:) = lambdarInt(i,:) + (a + (b-a).*rand(1,n)).*lambdarInt(i,:);
%   lambdavInt(i,:) = lambdavInt(i,:) + (a + (b-a).*rand(1,n)).*lambdavInt(i,:);
%  end
 
% Write initial guess for the solution (at all the mesh points)
solinit.x = tVec;
solinit.y(1:3,:) = rInt;
solinit.y(4:6,:) = vInt;
solinit.y(7:9,:) = lambdarInt;
solinit.y(10:12,:) = lambdavInt;


%--------------------------------------------------------------------------
% Call bvp4c
%--------------------------------------------------------------------------

%Compute the solution with bvp4c
options = bvpset('RelTol',1e-5,'AbsTol',1e-7,'NMax',1000);
sol = bvp4c(@(t,x) diffEq(t,x,mu),@(X0,Xf) transfer_bcs(X0,Xf,r0,rf,v0,vf),solinit,options);

%%Compute the solution with bvp4c using the scaled version
%options = bvpset('RelTol',1e-5,'AbsTol',1e-7,'NMax',50000);
%sol = bvp4c(@(t,x) diffEqScaled(t,x,mu,norm(r0),norm(v0),norm(lambdaAd(1:3,1)),norm(lambdaAd(4:6,1)),tf),@(X0,Xf) transfer_bcs(X0,Xf,r0,rf,v0,vf),solinit,options);
   
% %--------------------------------------------------------------------------
% %Cycle to refine the mesh
% %--------------------------------------------------------------------------
% 
%Compute the solution with bvp4c
% No = 500;
% for i=1:10
%  solinit.x = sol.x;
%  solinit.y = sol.y;
%  options = bvpset('RelTol',1e-4,'AbsTol',1e-6,'NMax',No*i);
%  sol = bvp4c(@(t,x) diffEq(t,x,mu),@(X0,Xf) transfer_bcs(X0,Xf,r0,rf,v0,vf),solinit,options);
% end

%--------------------------------------------------------------------------
% Write the output evaluated on a time vector of length equal to the mesh
% chosen by bvp4c
%--------------------------------------------------------------------------
% Save the original solution as output
tVec = linspace(t0,tf,length(sol.x));

% Interpolated it to get an equispaced time vector
solEval = deval(sol,tVec);

% Outputs
r = (solEval(1:3,:))';
v = (solEval(4:6,:))';
lambda = (solEval(7:12,:))';

% %Interpolate to get vdot
% for i=1:3
%    vdot(i,:) = spline(sol.x,sol.yp(i+3,:),tVec); 
% end
% vdot = vdot';

%--------------------------------------------------------------------------
% Write the output evaluated on the original time vector (size n)
%--------------------------------------------------------------------------
% Evaluate the solution in the original interpolated time vector
% solEval = deval(sol,tVec);
% 
% % Outputs
% r = (solEval(1:3,:))';
% v = (solEval(4:6,:))';
% lambda = (solEval(7:12,:))';

end

