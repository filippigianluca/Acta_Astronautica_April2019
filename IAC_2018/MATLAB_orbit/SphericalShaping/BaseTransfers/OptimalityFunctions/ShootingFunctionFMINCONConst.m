function [c,ceq] = ShootingFunctionFMINCONConst(y,inCond,fCond,mu,timeInt,M,nT)

%ShootingFunctionFMINCONConst (Optimal Control Function for fsolve)
% Function that calculates the equality constraints for the shooting
% function with fmincon
%
% INPUT
% y: matrix with grid points (position,velocity,adjoints) corresponding to
% the solution (dimension (12 x n) )
% inCond: struct with initial conditions
% fCond: struct with final conditions
% timeInt: time vector
% M: number of grid points
% nT: number of integration points within a phase
%
%OUTPUT
% c: [] no inequalities constraints
% ceq: equalities constraints. In this case they are the inital and final
% values of the state and the patching at the intervals boundaries

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

% Execute M integrations on the different subintervals
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
     [tOde,Y] = rk4(@diffEq,timeInteg,y(:,i),1,mu);
     Y = Y';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Saving of last elements for each sub-interval
    rOde(1:3,i) = Y(end,1:3)';
    vOde(1:3,i) = Y(end,4:6)';
    lambdaOde(1:6,i) = Y(end,7:12)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the constraints now that the integration in each sub-interval has
% been done. The constraints are written as the function to drive to zero

% Initial boundary constraints (r,v)
FunVal(1:3) = abs((r(:,1) - r0'));
FunVal(4:6) = abs((v(:,1) - v0'));

% Final boundary constraints (r,v)
FunVal(7:9) = abs((rOde(:,end) - rf'));
FunVal(10:12) = abs((vOde(:,end) - vf'));

% Add other 6 constraints, to link the last element of the input vectors
% (only in the case M+1 is passed, otherwise r(:,end) corresponds to
% rOde(:,end-1))
% FunVal(13:15) = abs((rOde(:,end) - r(:,end)));
% FunVal(16:18) = abs((vOde(:,end) - v(:,end)));
% 
% % Further linking (to obtain a square system when M+1 is passed)
%  FunVal(19:21) = abs((rf' - r(:,end)));
%  FunVal(22:24) = abs((vf' - v(:,end)));

% Add the constraints on the Lagrange multipliers (? nonsense?0
%FunVal(13:18) = abs((lambda(:,1) - lambda0'));

% Continuity boundary constraints(M-1)
for i=1:M-1
  % r 
  FunValr(3*(i-1)+1:3*(i-1)+3)  = abs((rOde(:,i) - r(:,i+1)));
  
  % v 
  FunValv(3*(i-1)+1:3*(i-1)+3)  = abs((vOde(:,i) - v(:,i+1)));
    
  % lambda 
  FunVall(6*(i-1)+1:6*(i-1)+6) = abs((lambdaOde(:,i) - lambda(:,i+1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%% Constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear inequalities
c = [];

% Nonlinear equalities
% Add everything in a single function vector
ceq = [FunVal FunValr FunValv FunVall];



end

