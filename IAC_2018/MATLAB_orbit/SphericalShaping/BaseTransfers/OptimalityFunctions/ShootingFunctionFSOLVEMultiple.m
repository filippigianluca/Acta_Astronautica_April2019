function [FunVal] = ShootingFunctionFSOLVEMultiple(y,inCond,fCond,mu,timeInt,M,nT)

%ShootingFunctionFSOLVEMultiple (Optimal Control Function for fsolve)
% Function that calculates difference between given and computed final
% state and patching conditions at the intervals boundaries,
% which is the function to drive to zero and find the correspondent
% solution vector lambda
%
% INPUT
% y: matrix with grid points (position,velocity,adjoints) corresponding to
% the solution (dimension (12 x n) )
% inCond: struct with initial conditions
% fCond: struct with final conditions
% mu: gravitational parameter
% timeVec: time vector (can be used for integration instead of t0,tf)
% M: number of grid points
% nT: number of integration points within a phase
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

% Extract from the matrix variable y the state matrix with r and v and the
% lambda adjoints matrix lambdar and lambdav
r = y(1:3,:);
v = y(4:6,:);
lambda(1:3,:) = y(7:9,:);
lambda(4:6,:) = y(10:12,:);

% Create matrices to keep the results of integration (only the last
% elements that are needed in order to write the continuity constraints)
rOde = [];
vOde = [];
lambdaOde = []; 

% Version modified without r(1:3,1) as optimization variable (modify also
% the starting i in the cycle)
% i=1;
% [tOde,Y] = ode45(@(t,x) diffEq(t,x,mu),linspace(timeInt(i),timeInt(i+1),nT+1),[r0 v0 lambda0]);
% rOde(1:3,i) = Y(end,1:3)';
% vOde(1:3,i) = Y(end,4:6)';
% lambdaOde(1:6,i) = Y(end,7:12)';
      
for i=1:M 
    
    % Initial conditions for each sub-interval
    r_in = r(:,i);
    v_in = v(:,i);
    lambda_in = lambda(:,i);
    
    % Integration on each sub-interval
    %[tVec,Y] = ode45(@(t,x) diffEq(t,x,mu),[timeInt(i) timeInt(i+1)],[r_in;v_in;lambda_in]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Change integration procedure respect to ode45
    
    %%% Euler, n steps %%%
    %timeInteg = linspace(timeInt(i),timeInt(i+1),nT+1);
    %[tt,Y] = euler(@diffEq,timeInteg,y(:,i),mu);
    %Y = Y';
    
     %%% Rk 4, n steps %%%
     timeInteg = linspace(timeInt(i),timeInt(i+1),nT+1);
     [tt,Y] = rk4(@diffEq,timeInteg,y(:,i),1,mu);
     Y = Y';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Saving of last elements for each sub-interval
    rOde(1:3,i) = Y(end,1:3)';
    vOde(1:3,i) = Y(end,4:6)';
    lambdaOde(1:6,i) = Y(end,7:12)';
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write the constraints now that the integration in each sub-interval has 
% been done. The constraints are written as the function to drive to zero

% Initial boundary constraints (r,v)
FunVal(1:3) = abs((r(:,1) - r0'));
FunVal(4:6) = abs((v(:,1) - v0'));

% Final boundary constraints (r,v)
FunVal(7:9) = abs((rOde(:,end) - rf'));
FunVal(10:12) = abs((vOde(:,end) - vf'));
%FunVal(1:3) = abs((rOde(:,end) - rf'));
%FunVal(4:6) = abs((vOde(:,end) - vf'));


% Add other 6 constraints, to link the last element of the input vectors
% (only in the case M+1 is passed, otherwise r(:,end) corresponds to
% rOde(:,end-1))
 %FunVal(13:15) = abs((rOde(:,end) - r(:,end)));
 %FunVal(16:18) = abs((vOde(:,end) - v(:,end)));
% 
% % Further linking (to obtain a square system when M+1 is passed)
%FunVal(7:9) = abs((r(:,end) - rf'));
%FunVal(10:12) = abs((v(:,end) - vf'));

% Continuity boundary constraints(M-1)
for i=1:M-1
  % r 
  FunValr(3*(i-1)+1:3*(i-1)+3)  = abs((rOde(:,i) - r(:,i+1))); 
  
  % v 
  FunValv(3*(i-1)+1:3*(i-1)+3)  = abs((vOde(:,i) - v(:,i+1)));
    
  % lambda 
  FunVall(6*(i-1)+1:6*(i-1)+6) = abs((lambdaOde(:,i) - lambda(:,i+1)));
  %FunVall(3*(i-1)+1:3*(i-1)+3) = abs((lambdaOde(4:6,i) - lambda(4:6,i+1)));
end

% Add everything in a single function vector
FunVal = [FunVal FunValr FunValv FunVall];



end

