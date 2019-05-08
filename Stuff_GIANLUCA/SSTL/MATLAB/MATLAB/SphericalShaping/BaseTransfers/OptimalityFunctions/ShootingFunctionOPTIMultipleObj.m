function [FunVal] = ShootingFunctionOPTIMultipleObj(y,M,timeInt,nT,mu)

%ShootingFunctionOPTIMultipleObj (Optimal Control Function for OPTI NLP)
% Fake objective function for NLP solver of OPTI
% INPUT
% y: array comprising the state and the control at each moment of time
% y = [r1..rn, v1..vn, u1..un] where r is a 3x1 vector (position), v is also
% a 3x1 vector (velocity) and u is a 3x1 vector (control)
% It is a matrix [3x3n x 1]
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

% Extract from the matrix variable y the state matrix with r and v and the
% lambda adjoints matrix lambdar and lambdav
r_V = y(1:3*M);
v_V = y(3*M+1:6*M);
lambdar_V = y(6*M+1:9*M);
lambdav_V = y(9*M+1:12*M);

for i=1:M
  rSol(:,i) = r_V(3*(i-1)+1:3*(i-1)+3);
  vSol(:,i) = v_V(3*(i-1)+1:3*(i-1)+3);
  lambdarSol(:,i) = lambdar_V(3*(i-1)+1:3*(i-1)+3);
  lambdavSol(:,i) = lambdav_V(3*(i-1)+1:3*(i-1)+3);
end

lambdaSol = [lambdarSol; lambdavSol];

% Integrate to obtain
for i=1:M 
    % Initial conditions for each sub-interval
    r_in = rSol(:,i);
    v_in = vSol(:,i);
    lambda_in = lambdaSol(:,i);

    % Integration on each sub-interval
     %%% Rk 4, n steps %%%
     timeInteg = linspace(timeInt(i),timeInt(i+1),nT+1);
     [tOde,Y] = rk4(@diffEq,timeInteg,[r_in;v_in;lambda_in],1,mu);
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

