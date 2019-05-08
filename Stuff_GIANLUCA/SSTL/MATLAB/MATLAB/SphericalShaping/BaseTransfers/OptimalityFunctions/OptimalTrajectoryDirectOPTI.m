function [r,v,u,tVec] = OptimalTrajectoryDirectOPTI(N,timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,Ucart,mu,uMax)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method
% Method: direct transcription solved with OPTI

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
tVec = linspace(t0,tf,N);
h = tVec(2)-tVec(1);

% Write initial and final conditions
inCond.r0 = r0;
inCond.v0 = v0;
inCond.t0 = t0;
fCond.rf = rf;
fCond.vf = vf;
fCond.tf = tf;

% An interpolation is needed in order to have the initial values of the var
 for i=1:3
  rInt(i,:) = spline(timeVec,R(i,:),tVec);
  vInt(i,:) = spline(timeVec,Vcart(i,:),tVec);
  uInt(i,:) = spline(timeVec,Ucart(i,:),tVec);
 end

% %Add some 'noise' to the initial guess data (a randomatic variation given
% %as a percentual of the original value)
% a = -0.1;
% b = 0.1;
%  for i=1:3
%   rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,N)).*rInt(i,:);
%   vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,N)).*vInt(i,:);
%   uInt(i,:) = uInt(i,:) + (a + (b-a).*rand(1,N)).*uInt(i,:);
%  end

% Collect in initial guess x0 in a vector form (requirement of OPTI)
for i=1:N
  rInVec(3*(i-1)+1:3*(i-1)+3) =  rInt(:,i);
  vInVec(3*(i-1)+1:3*(i-1)+3) =  vInt(:,i);
  uInVec(3*(i-1)+1:3*(i-1)+3) =  uInt(:,i);
end

x0 = [rInVec vInVec uInVec];


%--------------------------------------------------------------------------
% Adapt to formalism of OPTI
%--------------------------------------------------------------------------

% Write limits of nonlinear constraints
cl = zeros((12+6*(N-1)),1);
cu = zeros((12+6*(N-1)),1);

% Options for OPTI NLP
%opts = optiset('solver','filtersd','display','iter','maxiter',5000,'maxfeval',10000,'tolrfun',1e-10,'tolafun',1e-10,'tolint',1e-7);
opts = optiset('solver','ipopt','display','iter');
%opts = optiset('solver','nlopt','display','iter');

% Options for OPTI GNLP (Global NLP)
%opts = optiset('solver','nomad','display','iter'); % Only works for problems with inequality constraints
%opts = optiset('solver','pswarm','display','iter');
%opts = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'),'display','iter');
%opts = optiset('solver','scip','display','iter');

% Build OPTI Problem
Opt = opti('fun',@(y) DirTranscFunctionFMINCONVecObj(y,N,tVec),'nl',@(y) DirTranscFunctionFMINCONVecConst(y,inCond,fCond,N,h,mu,uMax),cl,cu,'x0',x0,'options',opts);

% Solve with OPTI
[xSol,fval,exitflag,info] = solve(Opt);

% Extract the obtained state and save it in matrices
rAft = xSol(1:3*N);
vAft = xSol(3*N+1:6*N);
uAft = xSol(6*N+1:9*N);

for i=1:N
  r(:,i) = rAft(3*(i-1)+1:3*(i-1)+3);
  v(:,i) = vAft(3*(i-1)+1:3*(i-1)+3);
  u(:,i) = uAft(3*(i-1)+1:3*(i-1)+3);
end



end