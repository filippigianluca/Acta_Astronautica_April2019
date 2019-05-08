function [FunVal] = ShootingFunctionFMINCONObj(y,M,timeInt,nT,mu)

%ShootingFunctionFMINCONObj
% Function that calculates discretized performance index for the direct
% transcription problem of the optimal control problem applied to low
% thrust trajectory optimization
%
% INPUT
% y: matrix with grid points (position,velocity,adjoints) corresponding to
% the solution (dimension (12 x n) )
% timeInt: time vector
% M: number of grid points
% nT: number of integration points within a phase
%
%OUTPUT
% FunVal: value of the function to minimize, in this case sum of the square
% of the control vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the function of which to find the root

% Obtain the solution (as matrix with initial values of r,v,lambda in each
% sub-interval)
rSol = y(1:3,:);
vSol = y(4:6,:);
lambdaSol(1:3,:) = y(7:9,:);
lambdaSol(4:6,:) = y(10:12,:);

% Integrate to obtain
for i=1:M 
    % Initial conditions for each sub-interval
    r_in = rSol(:,i);
    v_in = vSol(:,i);
    lambda_in = lambdaSol(:,i);

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
    
    % Saving of elements for each sub-interval (exclude the last point
    % because will be the first point of the next sub-interval integration) 
    tVec(nT*(i-1)+1:nT*(i-1)+nT) = tOde(1:end-1);
    %r(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,1:3);
    %v(nT*(i-1)+1:nT*(i-1)+nT,1:3) = Y(1:end-1,4:6);
    lambda(nT*(i-1)+1:nT*(i-1)+nT,1:6) = Y(1:end-1,7:12);
   
end

% Add the last element of the integration in the last sub-interval (because
% is not saved in the cycle), which corresponds to the last time point
tVec(end+1) = tOde(end);
%r(end+1,1:3) = Y(end,1:3);
%v(end+1,1:3) = Y(end,4:6);
lambda(end+1,1:6) = Y(end,7:12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract control history from integration results 
 u = -lambda(:,4:6)';

%Calculate the value of the function to minimize with CS quadrature
tVec = scalingFunction(tVec,1/(tVec(end)-tVec(1)));
tStep = tVec(3) - tVec(1);
uSquared = u(1,:).^2 + u(2,:).^2 + u(3,:).^2;
DV = (tStep/6)*(uSquared(1) + 2*sum( uSquared(3:2:end-1) ) + ...
       + 4*sum(uSquared(2:2:end-1)) + uSquared(end));
FunVal = 0.5*DV;

% Calculate the value of the function to minimize with just a summation
% N = length(tVec);
% uSquared = zeros(N,1);
% for i=1:N
%   uSquared(i) = norm(u(i,:))^2;
% end
% FunVal = 0.5*sum(uSquared);


%FunVal = 1;
end

