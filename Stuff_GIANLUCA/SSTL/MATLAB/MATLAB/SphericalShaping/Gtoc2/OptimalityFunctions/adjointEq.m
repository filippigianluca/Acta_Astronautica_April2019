function lambdadot = adjointEq(t,lambda,R,timeVec,mu)

%adjointEq: function that writes the adjoints differential equations for
%the optimal control problem (minimum Delta V)
%
% INPUT
% t: time
% lambda: Lagrange multipliers (r,v). Row vector
% R: radius vector in cartesian coordinates [x,y,z]. Is a matrix of size
% (3 x timeStep)
% timeVec: time vector of size (1 x timeStep)
% mu: gravitational parameter
% OUTPUT
% lambdadot: time derivative of Lagrange multipliers (r,v). Row vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, January 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write the time derivatives

% Find the index for selecting the vector r corresponding to that time step
index = find(timeVec == t);

if ( mod(index,1) == 0 )
    
  %disp('OK - Time in input present in timeVec'); 
  
  % Select the correspondent radius vector and transform into row vector
  r = R(:,index)';
  
else
    
  %warning('Error - Index doesnt exist! An interpolation is done in order to obtain the values for r');  
  %disp('Error - Index doesnt exist! An interpolation is done in order to obtain the values for r'); 
   
%   for i=2:length(timeVec)
%     if ( ( t < timeVec(i)) && ( t > timeVec(i-1)) )
%      indexSup = i;
%      break;
%     end
%   end
  
  % Interpolate to obtain r
  r(1) = interp1(timeVec,R(1,:),t);
  r(2) = interp1(timeVec,R(2,:),t);
  r(3) = interp1(timeVec,R(3,:),t);
  
end

% Write the lambda multipliers
lambdar = lambda(1:3);
lambdav = lambda(4:6);

% % Write the differential equations (from Vasile paper)
% lambdardot = -lambdav.*(3*mu/norm(r)^5*(r.^2) - mu/norm(r)^3*[1 1 1]);
% lambdavdot = -lambdar;

% Write the differential equations (from Homotopic approach and
% pseudospectral method...)
lambdardot = lambdav*mu/norm(r)^3 - 3*mu/norm(r)^5*(dot(lambdav,r))*r;
lambdavdot = -lambdar;

% Wrtie the time derivatives as outputs
lambdadot(1:3) = lambdardot;
lambdadot(4:6) = lambdavdot;

end

