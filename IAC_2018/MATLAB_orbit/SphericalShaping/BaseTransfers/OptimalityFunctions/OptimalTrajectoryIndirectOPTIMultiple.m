function [r,v,lambda,tVec] = OptimalTrajectoryIndirectOPTIMultiple(timeVec,t0,tf,r0,rf,v0,vf,R,Vcart,lambdaAd,M,nT,mu)
% Function that calculates the optimal trajectory of the base transfer to
% check the results of the spherical shaping method

% Method: indirect multiple shooting with OPTI

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

% Write initial and final conditions
inCond.r0 = r0;
inCond.v0 = v0;
inCond.t0 = t0;
fCond.rf = rf;
fCond.vf = vf;
fCond.tf = tf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write initial guess of solution vector, which is made by M*4 vector
% variables (M*4*3 scalar variables)

% First interpolate to obtain the values at such time instants (with
% initial and final instants)
timeInt = linspace(t0,tf,M+1);

% An interpolation is needed in order to have the values
 for i=1:3
  rInt(i,:) = spline(timeVec,R(i,:),timeInt);
  vInt(i,:) = spline(timeVec,Vcart(i,:),timeInt);
  lambdarInt(i,:) = spline(timeVec,lambdaAd(i,:),timeInt);
  lambdavInt(i,:) = spline(timeVec,lambdaAd(i+3,:),timeInt);
 end

% %  Add some 'noise' to the initial guess data (a randomatic variation given
% % as a percentual of the original value)
% a = -0.001;
% b = 0.001;
%  for i=1:3
%   rInt(i,:) = rInt(i,:) + (a + (b-a).*rand(1,M+1)).*rInt(i,:);
%   vInt(i,:) = vInt(i,:) + (a + (b-a).*rand(1,M+1)).*vInt(i,:);
%   lambdarInt(i,:) = lambdarInt(i,:) + (a + (b-a).*rand(1,M+1)).*lambdarInt(i,:);
%   lambdavInt(i,:) = lambdavInt(i,:) + (a + (b-a).*rand(1,M+1)).*lambdavInt(i,:);
%  end
 

% Collect in initial guess x0 in a vector form (requirement of OPTI)
for i=1:M
  rInVec(3*(i-1)+1:3*(i-1)+3) =  rInt(:,i);
  vInVec(3*(i-1)+1:3*(i-1)+3) =  vInt(:,i);
  lambdarInVec(3*(i-1)+1:3*(i-1)+3) =  lambdarInt(:,i);
  lambdavInVec(3*(i-1)+1:3*(i-1)+3) =  lambdavInt(:,i);
end

x0 = [rInVec vInVec lambdarInVec lambdavInVec];


%--------------------------------------------------------------------------
% Adapt to formalism of OPTI
%--------------------------------------------------------------------------

%%%% OPTI NLEQ %%%%
% Options for OPTI NLEQ
% opts = optiset('solver','filtersd','display','iter','maxiter',5000,'maxfeval',10000,'tolrfun',1e-10,'tolafun',1e-10,'tolint',1e-7);
% opts = optiset('solver','ipopt','display','iter');
% 
% Create OPTI Object
% Opt = opti('nleq',@(y) ShootingFunctionOPTIMultipleConst(y,inCond,fCond,mu,timeInt,M,nT),'x0',x0);
% 
% Solve the SNLE problem
% [xSol,fval,exitflag,info] = solve(Opt);


%%%% OPTI NLP %%%%
% Options for OPTI NLP
%opts = optiset('solver','filtersd','display','iter','maxiter',5000,'maxfeval',10000,'tolrfun',1e-10,'tolafun',1e-10,'tolint',1e-7);
opts = optiset('solver','ipopt','display','iter','maxiter',5000,'maxfeval',10000,'tolrfun',1e-8,'tolafun',1e-8,'tolint',1e-7);
%opts = optiset('solver','nlopt','display','iter');

% Write limits of nonlinear constraints
cl = zeros((12+12*(M-1)),1);
cu = zeros((12+12*(M-1)),1);

% Build OPTI Problem
Opt = opti('fun',@(y) ShootingFunctionOPTIMultipleObj(y,M,timeInt,nT,mu),'nl',@(y) ShootingFunctionOPTIMultipleConst(y,inCond,fCond,mu,timeInt,M,nT),cl,cu,'x0',x0,'options',opts);

% Solve with OPTI
[xSol,fval,exitflag,info] = solve(Opt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the obtained state and save it in matrices
rAft = xSol(1:3*M);
vAft = xSol(3*M+1:6*M);
lambdarAft = xSol(6*M+1:9*M);
lambdavAft = xSol(9*M+1:12*M);

for i=1:M
  rSol(:,i) = rAft(3*(i-1)+1:3*(i-1)+3);
  vSol(:,i) = vAft(3*(i-1)+1:3*(i-1)+3);
  lambdarSol(:,i) = lambdarAft(3*(i-1)+1:3*(i-1)+3);
  lambdavSol(:,i) = lambdavAft(3*(i-1)+1:3*(i-1)+3);
end

% Version modified without r(1:3,1) as optimization variable (modify also
% the starting i in the cycle)
% i=1;
% [tOde,Y] = ode45(@(t,x) diffEq(t,x,mu),linspace(timeInt(i),timeInt(i+1),nT+1),[r0 v0 lambda0]);
% tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(1:end-1);
% r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,1:3);
% v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,4:6);
% lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(1:end-1,7:12);
      
% Execute M integrations on the different subintervals
for i=1:M 
    % Initial conditions for each sub-interval
    r_in = rSol(:,i);
    v_in = vSol(:,i);
    lambdar_in = lambdarSol(:,i);
    lambdav_in = lambdavSol(:,i);
    
    % Integration on each sub-interval
    %[tOde,Y] = ode45(@(t,x) diffEq(t,x,mu),linspace(timeInt(i),timeInt(i+1),nT+1),[r_in;v_in;lambda_in]);
    timeInteg = linspace(timeInt(i),timeInt(i+1),nT+1);
    [tOde,Y] = rk4(@diffEq,timeInteg,[r_in ;v_in;lambdar_in;lambdav_in],1,mu);
    Y = Y';
     
    % Saving of elements for each sub-interval (exclude the last point
    % because will be the first point of the next sub-interval integration) 
    tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(1:end-1);
    r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,1:3);
    v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,4:6);
    lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(1:end-1,7:12);
   
%     % Different way: saving of elements for each sub-interval (exclude the first point
%     % because will be the last point of the previous sub-interval integration) 
%     tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(2:end);
%     r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(2:end,1:3);
%     v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(2:end,4:6);
%     lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(2:end,7:12);
   
end

% Add the last element of the integration in the last sub-interval (because
% is not saved in the cycle), which corresponds to the last time point
tVec(end+1) = tOde(end);
r(end+1,1:3) = Y(end,1:3);
v(end+1,1:3) = Y(end,4:6);
lambda(end+1,1:6) = Y(end,7:12);

% %Add the first element of the integration in the first sub-interval (because
% %is not saved in the cycle), which corresponds to the first time point
% tVec = [t0 tVec];
% r = [rSol(1:3,1)';r];
% v = [vSol(1:3,1)';v];
% lambda = [lambdaSol(1:6,1)';lambda];


end

