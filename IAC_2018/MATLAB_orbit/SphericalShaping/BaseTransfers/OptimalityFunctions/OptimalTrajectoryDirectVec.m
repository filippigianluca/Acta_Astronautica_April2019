function [r,v,u,tVec] = OptimalTrajectoryDirectVec(N,timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,Ucart,mu)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method
% Method: direct transcription with variables in vector form

% INPUT
% N: number of grid points
% timeVec: vector with the time of the trajectory
% t0,tf,r0,rf,v0,vf: initial and final time,position,velocity
% R: radius vector
% Vcart: velocity
% Ucart: control acceleration
% mu: gravitational parameter
% uMax: maximum control
%
% OUTPUT
% r: matrix with position
% v: matrix with velocity
% u: matrix with control acceleration
% tVec: vector with time
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters and variables for direct transcription   

% Find time step from number of nodes and time vector
%h = (tf - t0)/(N-1);
tVec = linspace(t0,tf,N);
h = tVec(2)-tVec(1);

% Write initial and final conditions
inCond.r0 = r0;
inCond.v0 = v0;
inCond.t0 = t0;
fCond.rf = rf;
fCond.vf = vf;
fCond.tf = tf;

% Write initial guess vector: it must be similar in form to y that is found
% by fmincon, so comprising the state and the control at each moment of time
% y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
% a 3x1 vector (velocity) and u is a 3x1 vector (control)
% It is a matrix [3x3n x 1]

% An interpolation is needed in order to have the values
 for i=1:3
  rInt(i,:) = spline(timeVec,R(i,:),tVec);
  vInt(i,:) = spline(timeVec,Vcart(i,:),tVec);
  uInt(i,:) = spline(timeVec,Ucart(i,:),tVec);
 end

for i=1:N
  rInVec(3*(i-1)+1:3*(i-1)+3) =  rInt(:,i);
  vInVec(3*(i-1)+1:3*(i-1)+3) =  vInt(:,i);
  uInVec(3*(i-1)+1:3*(i-1)+3) =  uInt(:,i);
end

x0 = [rInVec vInVec uInVec];

% Find the correct initial conditions for the lagrange multipliers from
% fsolve
options = optimoptions('fmincon','Algorithm','sqp','TolCon',1e-6,'TolFun',1e-6,'MaxFunEvals',800);
xSol = fmincon(@(y) DirTranscFunctionFMINCONVecObj(y,N),x0,[],[],[],[],[],[],@(y) DirTranscFunctionFMINCONVecConst(y,inCond,fCond,N,h,mu),options);

% Extract the obtained state and save it in matrices
rAft = xSol(1:3*N);
vAft = xSol(3*N+1:6*N);
uAft = xSol(6*N+1:9*N);

for i=1:N
  r(:,i) = rAft(3*(i-1)+1:3*(i-1)+3);
  v(:,i) = vAft(3*(i-1)+1:3*(i-1)+3);
  u(:,i) = uAft(3*(i-1)+1:3*(i-1)+3);
end

% % Print the difference between obtained final values and desired final
% % values
% format long;
% fprintf('The obtained final position is: %f  \n',r(:,end)); 
% fprintf('The desired final position is: %f  \n',rf); 
% fprintf('The obtained final velocity is: %f  \n',v(:,end)); 
% fprintf('The desired final velocity is: %f  \n',vf); 


  
end